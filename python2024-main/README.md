
---

# Algorithmic Bioinformatics – Python Coursework (2024)

This repository contains my solutions to **Assignments 01–04** and a **multi-part programming project** on **Atomic Force Microscopy (AFM)** data analysis.

Grading uses **automated private tests**. Each task includes at least one public example test you can run locally (see *Testing*).

---

## Environment

* Python 3.10+ (uses `match/case` in A03/T3)
* NumPy (A04 & project; `mamba install numpy` or `pip install numpy`)
* Matplotlib (project plotting)
* Streamlit + Plotly (project Part III app)

  ```bash
  pip install numpy matplotlib streamlit plotly biopython pandas logomaker
  ```

Repo layout (key parts):

```
python2024-main/
├── assignments/
│   ├── assignments_01.py ... assignments_04.py
│   ├── test_assignment_01.py ... test_assignment_04.py
│   └── sample inputs for tests (e.g., ass_02.3.fasta)
└── project-afm/
    ├── sample.txt              # provided 12-measurement sample
    ├── plotafm.py              # CLI for Parts I & II
    └── afm.py                  # Streamlit app for Part III
```

---

## Testing & Grading

* Private tests grade your submission. Public example tests are in `assignments/`.
* Run with `pytest -q` (recommended) or `python test_assignment_XX.py`.
* **Project Part II formatting rule:** the program must print **only** the per-spectrum lines (and any diagnostic lines must start with `#`). Otherwise, auto-tests will fail.

---

## Project: AFM (Atomic Force Microscopy)

### Overview

We analyze force–distance curves measured by an AFM over a **128×128** grid (two series per pixel: push **s=0** and lift **s=1**). Each measurement has \~700 points of **distance (m)** and **force (N)**.

* **Part I:** Parse the text data; plot each measurement.
* **Part II:** Estimate the **slope** of the linear-looking left region (esp. series **s=0**); add the fitted line to plots; print `s i j slope`.
* **Part III:** Build a **Streamlit** app to visualize **heights** and **slopes** as heatmaps and inspect any curve.

### Data

* **Sample:** `project-afm/sample.txt` (12 measurements).
* **Full dataset (provided by course):**

  * `afm.heights.npy` (height per pixel)
  * `afm.data.pickled` (all curves for 2×128×128 positions)
  * Place them next to `afm.py` and pass the **prefix** (e.g., `afm`).

---

### Part I — Parse & Plot

**Goal:** Extract distance `d` (m) and force `f` (N) for each `(s,i,j)`; plot **f vs d**.

**CLI (implemented in `project-afm/plotafm.py`):**

```bash
python project-afm/plotafm.py --textfile project-afm/sample.txt --show --plotprefix curve
```

* Creates 12 PNGs: `curve-<s>-<iii>-<jjj>.png` (e.g., `curve-0-004-000.png`)
* `--show` displays plots interactively (omit for batch run)
* X-axis: **distance (m)**, Y-axis: **force (N)**, labeled with meaningful titles

**Parsing notes implemented:**

* Lines starting with `#` are metadata.
* New measurement begins at `index/iIndex/jIndex` block.
* Two series per `(i,j)`: first `s=0` (push), second `s=1` (lift).
* Measurement data begins after `# units`; the first two columns are used (d,f).
* Number of rows from `# recorded-num-points`.

---

### Part II — Slope Estimation & Reporting

**Goal:** Estimate the slope of the approximately linear region on the **left part** of each curve (esp. `s=0`), overlay the fitted line on the plot, and print one summary line per spectrum.

**CLI:**

```bash
# pure text output with required formatting
python project-afm/plotafm.py -t project-afm/sample.txt

# also write plots
python project-afm/plotafm.py -t project-afm/sample.txt --plotprefix curve
```

**Output format (required):**

```
# parsing sample.txt...
# processing 12 spectra...
0 000 000 -0.01065
1 000 000 -0.01108
0 001 000 -0.01079
...
```

* Each non-comment line is: `series i j slope`
* Any extra info must start with `#`

**Slope method (implemented, reproducible):**

* Focus on the **lower-distance** section (left side).
* Smooth (optional) and evaluate a sliding window; pick the window with **highest linearity** (max |corr| / R²).
* Fit `force = a*distance + b` in that window via least squares; report **slope = a**.
* Works for both `s=0` and `s=1`.

*(Your exact values may differ slightly from the example; that’s acceptable.)*

---

### Part III — Streamlit App (Heatmaps + Curves)

**Goal:** Load preprocessed data, estimate slopes for all pixels (for a chosen series), visualize:

* **Upper heatmap:** heights (from `afm.heights.npy`)
* **Lower heatmap:** slopes (computed per pixel for selected **s**)
* **Curve viewer:** the f–d curve at `(s,i,j)` with the fitted line

**Run:**

```bash
streamlit run project-afm/afm.py -- afm
# where 'afm' is the dataset prefix; the app loads:
#   afm.heights.npy
#   afm.data.pickled
```

**App features:**

* Sidebar controls: **series s (0/1)**, **i** (0–127), **j** (0–127)
* Heatmaps via **Plotly** (or Matplotlib fallback)
* Curve plot with fitted line via Matplotlib
* (Bonus, optional): click a heatmap to set `(i,j)` using `st.session_state` + `streamlit_plotly_events`

---

## Assignments (Brief)

> Full task specs are in the `assignments/` folder; highlights below.

**Assignment 01**

* `process_list` (odd/even transform + 7-multiple position offset)
* `DNA_complement` (reverse complement over A,C,G,T)
* `list_kmers`, `number_of_unique` (k-mer utilities)
* `kmer_code`, `kmer_decode` (base-4 bitwise encoding)

**Assignment 02**

* `anagrams` (largest anagram group)
* `hamming_distance` (min/max across file, RNA→DNA, validate input)
* `common_kmers` (FASTA, set intersection)
* `is_vampire`, `vampire_generator` (vampire numbers)

**Assignment 03**

* Parcel stacks: `read_config`, `rearrange_parcels`
* PWM & consensus: `pwm`, `consensus_sequence`
* Mini interpreter: `solve(n, "ops")` using `match/case`

**Assignment 04** (uses NumPy)

* `CountingBloomFilter` class (uint8 counters, insert/delete/contains/count)
* `get_treasure` (circular list reordering puzzle)
* `quantile_normalize` (rank-based normalization)

---

## License & Author

* Author: **Abdel (Helala)** — Algorithmic Bioinformatics, 2024

---
