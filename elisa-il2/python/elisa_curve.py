"""IL-2 sandwich ELISA: standard-curve fitting and back-calculation.

Fits three calibration models to the IL-2 standards (0-600 pg/ml):
  * linear           OD = m*c + b
  * semi-log linear  OD = m*log10(c) + b          (standards > 0 only)
  * 4PL              OD = d + (a - d) / (1 + (c/C)**B)

then back-calculates the unknown samples, checks replicate CVs and
standard recovery, and estimates the limit of detection (LOD) and
quantification (LOQ).

Usage:
    python elisa_curve.py                       # uses ../data, writes ../results
    python elisa_curve.py --dilution 101
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

STD_LOW, STD_HIGH = 3.1, 600.0   # lowest / highest non-zero standard (pg/ml)


# ----------------------------------------------------------------------------
# Models
# ----------------------------------------------------------------------------
def four_pl(x, a, b, c, d):
    """4-parameter logistic. a = OD at zero, d = plateau, c = EC50, b = slope."""
    x = np.asarray(x, dtype=float)
    return d + (a - d) / (1.0 + np.power(x / c, b))


def four_pl_inverse(y, a, b, c, d):
    """Concentration for a given OD. NaN if OD lies outside the asymptotes."""
    y = np.asarray(y, dtype=float)
    lo, hi = min(a, d), max(a, d)
    with np.errstate(invalid="ignore", divide="ignore"):
        x = c * np.power((a - d) / (y - d) - 1.0, 1.0 / b)
    return np.where((y > lo) & (y < hi), x, np.nan)


def r_squared(y, y_hat):
    ss_res = np.sum((y - y_hat) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    return 1.0 - ss_res / ss_tot


def fit_models(conc: np.ndarray, od: np.ndarray) -> dict:
    """Fit all three models on individual replicates (not means)."""
    models = {}

    # Linear
    m, b = np.polyfit(conc, od, 1)
    pred = m * conc + b
    models["linear"] = dict(
        params=dict(m=m, b=b),
        fwd=lambda c, m=m, b=b: m * np.asarray(c, float) + b,
        inv=lambda y, m=m, b=b: (np.asarray(y, float) - b) / m,
        r2=r_squared(od, pred),
        rmse=np.sqrt(np.sum((od - pred) ** 2) / (len(od) - 2)),
    )

    # Semi-log linear (log10 of concentration; zero standard excluded)
    pos = conc > 0
    ms, bs = np.polyfit(np.log10(conc[pos]), od[pos], 1)
    pred = ms * np.log10(conc[pos]) + bs
    models["semilog"] = dict(
        params=dict(m=ms, b=bs),
        fwd=lambda c, m=ms, b=bs: m * np.log10(np.asarray(c, float)) + b,
        inv=lambda y, m=ms, b=bs: 10 ** ((np.asarray(y, float) - b) / m),
        r2=r_squared(od[pos], pred),
        rmse=np.sqrt(np.sum((od[pos] - pred) ** 2) / (pos.sum() - 2)),
    )

    # 4PL (Levenberg-Marquardt / trust-region via scipy)
    p0 = [od.min(), 1.0, np.exp(np.mean(np.log(conc[pos]))), od.max() * 1.15]
    bounds = ([-np.inf, 0.05, 1e-3, -np.inf], [np.inf, 10.0, 1e5, np.inf])
    popt, _ = curve_fit(four_pl, conc, od, p0=p0, bounds=bounds, maxfev=20000)
    pred = four_pl(conc, *popt)
    models["4pl"] = dict(
        params=dict(zip("aBCd", popt)),
        fwd=lambda c, p=popt: four_pl(c, *p),
        inv=lambda y, p=popt: four_pl_inverse(y, *p),
        r2=r_squared(od, pred),
        rmse=np.sqrt(np.sum((od - pred) ** 2) / (len(od) - 4)),
        popt=popt,
    )
    return models


# ----------------------------------------------------------------------------
# Detection limit
# ----------------------------------------------------------------------------
def detection_limits(std: pd.DataFrame, models: dict) -> dict:
    """Blank + 3/10 SD if blank wells exist, otherwise ICH Q2 on low standards."""
    fpl = models["4pl"]
    blank = std.loc[std.conc_pgml == 0, ["od1", "od2"]].to_numpy().ravel()
    blank = blank[~np.isnan(blank)]

    if blank.size >= 2:
        mu, s = blank.mean(), blank.std(ddof=1)
        lod = float(fpl["inv"](mu + 3 * s))
        loq = float(fpl["inv"](mu + 10 * s))
        method = "blank + 3 SD / 10 SD (4PL)"
    else:
        low = std[(std.conc_pgml > 0) & (std.conc_pgml <= 50)]
        x = np.repeat(low.conc_pgml.to_numpy(), 2)
        y = low[["od1", "od2"]].to_numpy().ravel()
        m, b = np.polyfit(x, y, 1)
        s_res = np.sqrt(np.sum((y - (m * x + b)) ** 2) / (len(x) - 2))
        lod, loq = 3.3 * s_res / m, 10 * s_res / m
        method = "ICH Q2: 3.3*s/slope, 10*s/slope (standards <= 50 pg/ml)"

    a = fpl["popt"][0]
    lod_4pl = float(fpl["inv"](a + 3 * fpl["rmse"]))
    return dict(method=method, lod=lod, loq=loq, lod_4pl=lod_4pl)


# ----------------------------------------------------------------------------
# Analysis per person
# ----------------------------------------------------------------------------
def cv_percent(a, b):
    return np.std([a, b], ddof=1) / np.mean([a, b]) * 100


def analyse(person, std, smp, dilution, out_dir):
    have_blank = std.loc[std.conc_pgml == 0, ["od1", "od2"]].notna().all(axis=None)
    use = std if have_blank else std[std.conc_pgml > 0]
    conc = np.repeat(use.conc_pgml.to_numpy(float), 2)
    od = use[["od1", "od2"]].to_numpy(float).ravel()

    models = fit_models(conc, od)
    lim = detection_limits(std, models)
    fpl = models["4pl"]

    # --- samples -----------------------------------------------------------
    rows = []
    for _, s in smp.iterrows():
        mean_od = (s.od1 + s.od2) / 2
        c4 = float(fpl["inv"](mean_od))
        if np.isnan(c4):
            status = "outside asymptotes"
        elif c4 > STD_HIGH:
            status = "above top standard"
        elif c4 < lim["lod"]:
            status = "below LOD"
        elif c4 < lim["loq"]:
            status = "below LOQ"
        else:
            status = "in range"
        rows.append(dict(
            person=person, sample=s["sample"], mean_od=mean_od,
            od_cv_pct=cv_percent(s.od1, s.od2),
            pgml_4pl=c4,
            pgml_4pl_rep1=float(fpl["inv"](s.od1)),
            pgml_4pl_rep2=float(fpl["inv"](s.od2)),
            pgml_semilog=float(models["semilog"]["inv"](mean_od)),
            pgml_linear=float(models["linear"]["inv"](mean_od)),
            supernatant_ngml_4pl=c4 * dilution / 1000,
            software_pgml_mean=(s.software_pgml_1 + s.software_pgml_2) / 2,
            status=status,
        ))
    results = pd.DataFrame(rows)

    # --- standard recovery -------------------------------------------------
    rec = []
    for _, r in std[std.conc_pgml > 0].iterrows():
        mean_od = (r.od1 + r.od2) / 2
        rec.append(dict(
            person=person, conc_pgml=r.conc_pgml, mean_od=mean_od,
            od_cv_pct=cv_percent(r.od1, r.od2),
            **{f"recovery_{k}_pct": float(models[k]["inv"](mean_od)) / r.conc_pgml * 100
               for k in ("4pl", "semilog", "linear")},
        ))
    recovery = pd.DataFrame(rec)

    # --- fit summary -------------------------------------------------------
    fit = pd.DataFrame([
        dict(person=person, model=k, r2=v["r2"], rmse_od=v["rmse"],
             **{f"p_{n}": val for n, val in v["params"].items()})
        for k, v in models.items()
    ])
    fit["lod_pgml"], fit["loq_pgml"], fit["lod_4pl_pgml"] = lim["lod"], lim["loq"], lim["lod_4pl"]
    fit["lod_method"] = lim["method"]

    plot_curve(person, conc, od, models, smp, lim, out_dir)
    return results, recovery, fit


def in_notebook() -> bool:
    """True inside Jupyter, Google Colab, VS Code notebooks, etc."""
    try:
        from IPython import get_ipython
        shell = get_ipython()
        return shell is not None and shell.__class__.__name__ != "TerminalInteractiveShell"
    except ImportError:
        return False


def show_table(df: pd.DataFrame, title: str):
    """Rich table in notebooks, plain text in a terminal."""
    df = df.round(3)
    if in_notebook():
        from IPython.display import display, Markdown
        display(Markdown(f"### {title}"))
        display(df.reset_index(drop=True))
    else:
        print(f"\n{title}\n{df.to_string(index=False)}")


def project_root() -> Path:
    """Folder that contains data/. Works as a script and in notebooks (no __file__)."""
    try:
        return Path(__file__).resolve().parent.parent
    except NameError:
        cwd = Path.cwd()
        for cand in (cwd, cwd / "elisa-il2", cwd.parent):
            if (cand / "data" / "standards.csv").exists():
                return cand
        return cwd


def plot_curve(person, conc, od, models, smp, lim, out_dir):
    fig, ax = plt.subplots(figsize=(7, 4.6))
    xx = np.logspace(np.log10(2.2), np.log10(800), 300)
    ax.plot(xx, models["4pl"]["fwd"](xx), lw=2.2, color="#1f5fbf", label="4PL")
    ax.plot(xx, models["semilog"]["fwd"](xx), "--", lw=1.6, color="#c2410c", label="Semi-log linear")
    ax.plot(xx, models["linear"]["fwd"](xx), ":", lw=1.6, color="#6b7a8c", label="Linear")
    pos = conc > 0
    ax.plot(conc[pos], od[pos], "s", mfc="white", mec="black", ms=7, label="Standards")

    for _, s in smp.iterrows():
        y = (s.od1 + s.od2) / 2
        x = float(models["4pl"]["inv"](y))
        if np.isfinite(x):
            ax.plot([2.2, x, x], [y, y, 0], color="#8a6a00", lw=1, ls=(0, (3, 3)))
            ax.plot(x, y, "o", color="#f2c200", mec="#8a6a00", ms=8)
            ax.annotate(s["sample"], (x, y), xytext=(6, 6), textcoords="offset points")

    ax.axvspan(2.2, lim["lod"], color="grey", alpha=0.12, label=f"< LOD ({lim['lod']:.1f} pg/ml)")
    ax.set_xscale("log")
    ax.set_xlim(2.2, 800)
    ax.set_ylim(0, None)
    ax.set_xticks([3.1, 12.5, 50, 300, 600])
    ax.set_xticklabels(["3.1", "12.5", "50", "300", "600"])
    ax.set_xlabel("IL-2 (pg/ml, log scale)")
    ax.set_ylabel("OD 450 nm")
    ax.set_title(f"IL-2 standard curve, Person {person}  "
                 f"(4PL R² = {models['4pl']['r2']:.3f})")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(out_dir / f"curve_person_{person}.png", dpi=200)
    # figure is left open; plt.show() in main() displays all of them


def main():
    root = project_root()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", type=Path, default=root / "data")
    ap.add_argument("--out-dir", type=Path, default=root / "results" / "python")
    ap.add_argument("--dilution", type=float, default=100.0,
                    help="sample dilution factor before the assay (protocol: 10 µl + 1000 µl)")
    args, _ = ap.parse_known_args()   # ignores the extra arguments notebooks pass in
    args.out_dir.mkdir(parents=True, exist_ok=True)

    std_all = pd.read_csv(args.data_dir / "standards.csv")
    smp_all = pd.read_csv(args.data_dir / "samples.csv")

    res, rec, fit = [], [], []
    for person in sorted(std_all.person.unique()):
        r, c, f = analyse(person,
                          std_all[std_all.person == person].reset_index(drop=True),
                          smp_all[smp_all.person == person].reset_index(drop=True),
                          args.dilution, args.out_dir)
        res.append(r); rec.append(c); fit.append(f)

    res, rec, fit = pd.concat(res), pd.concat(rec), pd.concat(fit)
    res.to_csv(args.out_dir / "sample_results.csv", index=False, float_format="%.3f")
    rec.to_csv(args.out_dir / "standard_recovery.csv", index=False, float_format="%.3f")
    fit.to_csv(args.out_dir / "fit_summary.csv", index=False, float_format="%.5g")

    pd.set_option("display.width", 160)
    show_table(res[["person", "sample", "mean_od", "od_cv_pct", "pgml_4pl", "pgml_semilog",
                    "pgml_linear", "supernatant_ngml_4pl", "software_pgml_mean", "status"]],
               "IL-2 sample results (in-well pg/ml, supernatant ng/ml)")
    show_table(rec, "Standard recovery (%)")
    show_table(fit[["person", "model", "r2", "rmse_od", "lod_pgml", "loq_pgml", "lod_4pl_pgml"]],
               "Fit quality and detection limits")
    print(f"\nResults written to {args.out_dir}")
    plt.show()   # inline in notebooks, separate windows when run as a script


if __name__ == "__main__":
    main()
