"""lac operon induction kinetics from the shared course table.

Steps
  1. Read every shared table in data/shared/ (xlsx or csv, one per group, or
     one pooled file) and clean it: column-name aliases, decimal commas,
     spelling variants of carbon source / inducer / strain, empty and
     duplicate rows.
  2. Where a group's raw readings are available (data/lab_notebook/*.csv),
     recompute beta-galactosidase activity in Miller-type units
         U = A420 * 1000 / (t * V * A600)
     and cross-check the shared table: wrong condition labels, swapped rows
     and calculation errors are reported and corrected from the raw data.
  3. Flag readings outside the valid A420 range (> 1.0 per the script,
     "> 3" = photometer overflow -> lower bound only, < 0.02 = near the
     detection limit).
  4. Pool groups (mean, SD, n per strain x carbon source x inducer x time),
     compute fold changes and the glycerol/glucose and IPTG/no-IPTG ratios,
     and draw a multi-curve kinetics plot.

Usage
    python lac_kinetics.py
    %run python/lac_kinetics.py          (Jupyter / Spyder / Colab)
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ----------------------------------------------------------------------------
# Assay parameters (Praktikum Regulationsbiologie 2025)
# ----------------------------------------------------------------------------
T_REACTION_MIN = 10.0                         # ONPG reaction time
A600_DILUTION = 5.0                           # 200 µl suspension + 800 µl Z buffer
SUSPENSION_ML = {"glucose": 0.020, "glycerol": 0.010}   # volume in the enzyme assay
A420_MAX = 1.0                                # upper limit of the linear range (script)
A420_OVERFLOW = 3.0                           # photometer limit, "> 3"
A420_MIN = 0.020                              # below this: near the detection limit

COLUMN_ALIASES = {
    "gruppe": "group", "grp": "group",
    "nr": "no", "nr.": "no", "probe": "no", "sample": "no", "no.": "no",
    "stamm": "strain",
    "c-source": "csource", "c_source": "csource", "c-quelle": "csource",
    "kohlenstoffquelle": "csource", "medium": "csource",
    "zeit": "time", "t": "time", "t (min)": "time", "time (min)": "time", "time_min": "time",
    "induktor": "inducer",
    "activity": "betagal", "aktivität": "betagal", "aktivitaet": "betagal",
    "miller": "betagal", "units": "betagal", "beta-gal": "betagal",
}
CSOURCE_MAP = {"glu": "glucose", "glc": "glucose", "glucose": "glucose", "glukose": "glucose",
               "gly": "glycerol", "glyc": "glycerol", "glycerol": "glycerol", "glycerin": "glycerol",
               "lac": "lactose", "lactose": "lactose", "laktose": "lactose"}
NO_INDUCER = {"", "-", "–", "ohne", "none", "no", "kein", "0", "nan", "w/o"}


# ----------------------------------------------------------------------------
# Cleaning helpers
# ----------------------------------------------------------------------------
def to_number(value):
    """'1,347.83', '1347,83', ' 0.5 ' -> float. Returns (value, is_lower_bound)."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return np.nan, False
    s = str(value).strip().replace(" ", "")
    lower = s.startswith(">")
    s = s.lstrip("><=")
    if "," in s and ("." not in s or s.rfind(",") > s.rfind(".")):
        s = s.replace(".", "").replace(",", ".")   # German: 1.347,83 or 1347,83
    s = s.replace(",", "")                         # English thousands separator
    try:
        return float(s), lower
    except ValueError:
        return np.nan, False


def norm_csource(v):
    key = str(v).strip().lower()
    return CSOURCE_MAP.get(key, key)


def norm_inducer(v):
    key = str(v).strip().lower()
    if key in NO_INDUCER:
        return "none"
    return "IPTG" if "iptg" in key else key


def norm_strain(v):
    key = re.sub(r"[\s_]", "", str(v).strip().lower())
    if key.startswith("laci") and ("neg" in key or key.endswith("-")):
        return "lacI-"
    if key in {"wt", "wildtyp", "wildtype", "wildtypus"}:
        return "WT"
    if key.startswith("crp") and ("neg" in key or key.endswith("-")):
        return "crp-"
    return str(v).strip()


def read_shared(folder: Path) -> pd.DataFrame:
    frames = []
    for f in sorted(folder.glob("*")):
        if f.suffix.lower() in {".xlsx", ".xls"}:
            df = pd.read_excel(f, dtype=str)
        elif f.suffix.lower() == ".csv":
            df = pd.read_csv(f, dtype=str, sep=None, engine="python")
        else:
            continue
        df["source_file"] = f.name
        frames.append(df)
    if not frames:
        raise FileNotFoundError(f"No shared tables found in {folder}")
    raw = pd.concat(frames, ignore_index=True)

    raw.columns = [COLUMN_ALIASES.get(c.strip().lower(), c.strip().lower()) for c in raw.columns]
    raw = raw.dropna(how="all", subset=[c for c in raw.columns if c != "source_file"])

    out = pd.DataFrame({
        "group": raw["group"].map(lambda v: to_number(v)[0]).astype("Int64"),
        "no": raw["no"].map(lambda v: to_number(v)[0]).astype("Int64"),
        "strain": raw["strain"].map(norm_strain),
        "csource": raw["csource"].map(norm_csource),
        "inducer": raw["inducer"].map(norm_inducer),
        "time_min": raw["time"].map(lambda v: to_number(v)[0]),
        "betagal_shared": raw["betagal"].map(lambda v: to_number(v)[0]),
        "source_file": raw["source_file"],
    })
    before = len(out)
    out = out.drop_duplicates(subset=["group", "no", "strain", "csource", "inducer", "time_min"])
    out.attrs["duplicates_removed"] = before - len(out)
    return out.reset_index(drop=True)


def read_notebooks(folder: Path) -> pd.DataFrame:
    frames = [pd.read_csv(f, dtype=str) for f in sorted(folder.glob("*.csv"))]
    if not frames:
        return pd.DataFrame()
    nb = pd.concat(frames, ignore_index=True)
    a420 = nb["a420"].map(to_number)
    nb = pd.DataFrame({
        "group": nb["group"].astype(int), "no": nb["no"].astype(int),
        "strain": nb["strain"].map(norm_strain), "csource": nb["csource"].map(norm_csource),
        "inducer": nb["inducer"].map(norm_inducer), "time_min": nb["time_min"].astype(float),
        "a600_diluted": nb["a600_diluted"].map(lambda v: to_number(v)[0]),
        "a420": [v for v, _ in a420], "a420_overflow": [lb for _, lb in a420],
    })
    vol = nb["csource"].map(SUSPENSION_ML)
    nb["a600"] = nb["a600_diluted"] * A600_DILUTION
    nb["activity_U"] = nb["a420"] * 1000 / (T_REACTION_MIN * vol * nb["a600"])
    return nb


def quality_flag(row) -> str:
    if row.get("a420_overflow") is True or row.get("a420_overflow") is np.True_:
        return "lower bound (A420 > 3)"
    a = row.get("a420", np.nan)
    if pd.isna(a):
        return "not verified (no raw data)"
    if a > A420_MAX:
        return "A420 > 1.0 (outside linear range)"
    if a < A420_MIN:
        return "A420 < 0.02 (near detection limit)"
    return "ok"


# ----------------------------------------------------------------------------
# Cross-check shared table against raw readings
# ----------------------------------------------------------------------------
def reconcile(shared: pd.DataFrame, nb: pd.DataFrame):
    labels = ["strain", "csource", "inducer", "time_min"]
    issues = []
    if nb.empty:
        data = shared.rename(columns={"betagal_shared": "activity_U"}).copy()
        data["quality"] = "not verified (no raw data)"
        return data, pd.DataFrame(issues)

    m = shared.merge(nb, on=["group", "no"], how="outer", suffixes=("_shared", ""), indicator=True)
    for _, r in m.iterrows():
        g, n = r["group"], r["no"]
        if r["_merge"] == "left_only":
            continue
        if r["_merge"] == "right_only":
            issues.append(dict(group=g, no=n, issue="missing from shared table", detail=""))
            continue
        wrong = [c for c in labels if r[f"{c}_shared"] != r[c]]
        if wrong:
            shared_lbl = tuple(r[f"{c}_shared"] for c in labels)
            partner = m[(m.group == g) & (m.no != n) &
                        m.apply(lambda q: tuple(q[c] for c in labels) == shared_lbl, axis=1)]
            swap = ""
            if len(partner):
                p = partner.iloc[0]
                if tuple(p[f"{c}_shared"] for c in labels) == tuple(r[c] for c in labels):
                    swap = f"; swapped with no. {int(p['no'])}"
            issues.append(dict(
                group=g, no=n, issue="condition label wrong in shared table",
                detail=(f"shared says {r['csource_shared']}/{r['inducer_shared']}/{r['time_min_shared']:g} min, "
                        f"raw data says {r['csource']}/{r['inducer']}/{r['time_min']:g} min{swap}")))
        ref = r["activity_U"]
        if np.isfinite(ref) and abs(r["betagal_shared"] - ref) > 0.01 * abs(ref):
            issues.append(dict(group=g, no=n, issue="activity differs from recalculation",
                               detail=f"shared {r['betagal_shared']:.2f} U, recalculated {ref:.2f} U"))

    # final data: raw-data labels + recalculated activity where available
    verified = nb.copy()
    verified["source"] = "raw readings"
    have_raw = set(zip(nb.group.astype(int), nb.no.astype(int)))
    keep = [(int(g), int(n)) not in have_raw for g, n in zip(shared.group, shared.no)]
    unverified = shared[keep].rename(columns={"betagal_shared": "activity_U"})
    unverified["source"] = "shared table only"
    data = pd.concat([verified, unverified], ignore_index=True)
    data["quality"] = data.apply(quality_flag, axis=1)
    return data, pd.DataFrame(issues)


# ----------------------------------------------------------------------------
# Kinetics summary
# ----------------------------------------------------------------------------
def add_preinduction_baseline(data: pd.DataFrame) -> pd.DataFrame:
    """IPTG is added right after the t = 0 sample, so the uninduced t0 of the
    same strain / carbon source / group is the t0 of the IPTG curve too."""
    t0 = data[(data.time_min == 0) & (data.inducer == "none")].copy()
    have = data[data.inducer == "IPTG"][["group", "strain", "csource"]].drop_duplicates()
    keys = ["group", "strain", "csource"]
    base = t0.merge(have, on=keys)
    already = set(map(tuple, data.loc[(data.inducer == "IPTG") & (data.time_min == 0), keys].values))
    base = base[[tuple(k) not in already for k in base[keys].values]]
    base["inducer"] = "IPTG"
    base["source"] = "t0 before IPTG (shared)"
    return pd.concat([data, base], ignore_index=True)


def pooled_kinetics(data: pd.DataFrame) -> pd.DataFrame:
    usable = data[data.quality == "ok"]
    g = usable.groupby(["strain", "csource", "inducer", "time_min"])["activity_U"]
    out = g.agg(n="count", mean_U="mean", sd_U="std").reset_index()
    allpts = data.groupby(["strain", "csource", "inducer", "time_min"]).size().rename("n_all").reset_index()
    return allpts.merge(out, how="left").fillna({"n": 0})


def ratio(num, den):
    """Ratio with a text qualifier when one side is not a clean value."""
    if num is None or den is None:
        return np.nan, "missing"
    v = num["activity_U"] / den["activity_U"]
    qn, qd = num["quality"], den["quality"]
    if qn == "ok" and qd == "ok":
        return v, "reliable"
    if qn.startswith("lower bound") and qd == "ok":
        return v, "lower limit (≥)"
    return v, f"uncertain ({'numerator' if qn != 'ok' else 'denominator'}: {qn if qn != 'ok' else qd})"


def key_ratios(data: pd.DataFrame) -> pd.DataFrame:
    rows = []
    pick = lambda d, **k: (lambda s: None if s.empty else s.iloc[0])(
        d.loc[np.logical_and.reduce([d[c] == v for c, v in k.items()])])
    for (g, strain), d in data.groupby(["group", "strain"]):
        d = d[d.source != "t0 before IPTG (shared)"]
        v, q = ratio(pick(d, csource="glycerol", inducer="none", time_min=0),
                     pick(d, csource="glucose", inducer="none", time_min=0))
        rows.append(dict(group=g, strain=strain, comparison="glycerol / glucose (t = 0, catabolite repression)",
                         time_min=0, ratio=v, reliability=q))
        for cs in ("glucose", "glycerol"):
            for t in sorted(d.time_min[d.time_min > 0].unique()):
                v, q = ratio(pick(d, csource=cs, inducer="IPTG", time_min=t),
                             pick(d, csource=cs, inducer="none", time_min=t))
                rows.append(dict(group=g, strain=strain, comparison=f"IPTG / no IPTG on {cs}",
                                 time_min=t, ratio=v, reliability=q))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------------
COLORS = {"glucose": "#1f5fbf", "glycerol": "#c2410c", "lactose": "#2f855a"}
MARKERS = {"ok": ("o", True), "A420 < 0.02 (near detection limit)": ("o", False),
           "A420 > 1.0 (outside linear range)": ("s", False), "lower bound (A420 > 3)": ("^", True),
           "not verified (no raw data)": ("D", False)}


def plot_kinetics(data: pd.DataFrame, out_dir: Path):
    strains = sorted(data.strain.unique())
    fig, axes = plt.subplots(1, len(strains), figsize=(7.5 * len(strains), 5), squeeze=False)
    for ax, strain in zip(axes[0], strains):
        d = data[data.strain == strain]
        for (cs, ind), c in d.groupby(["csource", "inducer"]):
            pts = c.groupby("time_min").agg(U=("activity_U", "mean"), q=("quality", "first")).reset_index()
            color = COLORS.get(cs, "grey")
            ls = "-" if ind == "IPTG" else "--"
            ax.plot(pts.time_min, pts.U, ls, color=color, lw=1.8, label=f"{cs}, {'+ IPTG' if ind == 'IPTG' else 'no IPTG'}")
            for _, p in pts.iterrows():
                mk, filled = MARKERS.get(p.q, ("x", True))
                ax.plot(p.time_min, p.U, mk, ms=8, mec=color, mfc=color if filled else "white", mew=1.6)
                if p.q.startswith("lower bound"):
                    ax.annotate("", xy=(p.time_min, p.U * 2.2), xytext=(p.time_min, p.U),
                                arrowprops=dict(arrowstyle="->", color=color))
        ax.set_yscale("log")
        ax.set_xlabel("time after IPTG addition (min)")
        ax.set_ylabel("β-galactosidase activity (U = A420·1000 / (t·V·A600))")
        ax.set_title(f"{strain} strain (groups: {', '.join(map(str, sorted(d.group.unique())))})")
        ax.set_xticks(sorted(d.time_min.unique()))
        ax.grid(alpha=0.3, which="both")
        # legend: curves + marker meaning
        handles, labels = ax.get_legend_handles_labels()
        for q, (mk, filled) in MARKERS.items():
            if q in set(d.quality):
                handles.append(plt.Line2D([], [], ls="", marker=mk, mec="black",
                                          mfc="black" if filled else "white", ms=7))
                labels.append(q)
        ax.legend(handles, labels, fontsize=8, loc="best")
    fig.tight_layout()
    fig.savefig(out_dir / "lac_kinetics.png", dpi=200)
    return fig


# ----------------------------------------------------------------------------
# Display helpers (same as the ELISA project)
# ----------------------------------------------------------------------------
def in_notebook() -> bool:
    try:
        from IPython import get_ipython
        shell = get_ipython()
        return shell is not None and shell.__class__.__name__ != "TerminalInteractiveShell"
    except ImportError:
        return False


def show_table(df: pd.DataFrame, title: str):
    df = df.round(2)
    if in_notebook():
        from IPython.display import display, Markdown
        display(Markdown(f"### {title}"))
        display(df.reset_index(drop=True))
    else:
        print(f"\n{title}\n{df.to_string(index=False)}")


def project_root() -> Path:
    try:
        return Path(__file__).resolve().parent.parent
    except NameError:
        cwd = Path.cwd()
        for cand in (cwd, cwd / "lac-operon", cwd.parent):
            if (cand / "data" / "shared").exists():
                return cand
        return cwd


def main():
    root = project_root()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", type=Path, default=root / "data")
    ap.add_argument("--out-dir", type=Path, default=root / "results" / "python")
    args, _ = ap.parse_known_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    shared = read_shared(args.data_dir / "shared")
    nb = read_notebooks(args.data_dir / "lab_notebook")
    data, issues = reconcile(shared, nb)
    data = add_preinduction_baseline(data)
    data = data.sort_values(["group", "strain", "csource", "inducer", "time_min"]).reset_index(drop=True)
    pooled = pooled_kinetics(data)
    ratios = key_ratios(data)

    cols = ["group", "no", "strain", "csource", "inducer", "time_min", "a600", "a420", "activity_U", "quality", "source"]
    cols = [c for c in cols if c in data.columns]
    data[cols].to_csv(args.out_dir / "clean_data.csv", index=False)
    issues.to_csv(args.out_dir / "data_issues.csv", index=False)
    pooled.to_csv(args.out_dir / "pooled_kinetics.csv", index=False)
    ratios.to_csv(args.out_dir / "key_ratios.csv", index=False)
    plot_kinetics(data, args.out_dir)

    print(f"Read {len(shared)} rows from data/shared "
          f"({shared.attrs.get('duplicates_removed', 0)} duplicates removed), "
          f"groups: {sorted(shared.group.dropna().unique().tolist())}")
    show_table(issues if len(issues) else pd.DataFrame({"issue": ["none found"]}),
               "Problems found in the shared table")
    show_table(data[cols], "Cleaned data (labels and activity from raw readings where available)")
    show_table(ratios, "Key ratios")
    show_table(pooled, "Pooled kinetics (reliable points only)")
    print(f"\nResults written to {args.out_dir}")
    plt.show()


if __name__ == "__main__":
    main()
