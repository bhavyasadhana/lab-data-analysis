# Lab data analysis in Python and MATLAB

Bhavya Sadhana, B.Sc. Biosystems Engineering (Biosystemtechnik), Otto von Guericke University Magdeburg

Two small analysis projects I wrote in November 2026 to evaluate data from my own university laboratory courses. Each analysis exists twice, once in Python and once in MATLAB, and both versions implement the same steps.

| Project | Lab course | What the code does |
|---|---|---|
| [IL-2 ELISA analysis](elisa-il2/) | Immunology block course, August 2025 | Fits linear, semi-log and 4PL standard curves, back-calculates IL-2 in the samples, checks replicate CVs and standard recovery, estimates LOD and LOQ |
| [lac operon kinetics](lac-operon/) | Regulatory biology, June 2025 | Cleans the results table, recalculates beta-galactosidase activity from raw photometer readings, flags readings outside the valid range, plots the induction kinetics |

## Repository layout

```
elisa-il2/
  python/elisa_curve.py      Python version
  matlab/elisa_curve.m       MATLAB version
  data/                      plate reader values (standards and samples)
  results/                   tables and plots written by the scripts
lac-operon/
  python/lac_kinetics.py
  matlab/lac_kinetics.m
  data/shared/               results table as reported in the course
  data/lab_notebook/         raw photometer readings
  results/
```

## How to run

Python 3.9 or newer:

```
pip install -r requirements.txt
python elisa-il2/python/elisa_curve.py
python lac-operon/python/lac_kinetics.py
```

MATLAB: open the `matlab` folder of a project and run `elisa_curve` (R2016b or newer, no toolboxes needed) or `lac_kinetics` (R2020b or newer).

Each script writes CSV tables and PNG plots to `results/python` or `results/matlab` inside its project folder.

## Data

All data are my own measurements from laboratory courses at OVGU Magdeburg. The ELISA data contain the plates of my lab partner and me (Person A and Person B). The lac operon data are the readings of my own group only.

## Contact

[LinkedIn](https://www.linkedin.com/in/bhavya-sadhana/)
