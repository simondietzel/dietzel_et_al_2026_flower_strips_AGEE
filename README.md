# Data and Code for Dietzel et al., 2026

## Publication
Title: Synergistic effects of flower strips and landscape complexity buffer arthropods against warming and drought – A meta-analysis.
Agriculture, Ecosystems & Environment.

> Simon Dietzel, Péter Batáry, Alina Twerski, Anita Kirmer, Johannes Kollmann & Christina Fischer
> [https://doi.org/10.xxxx/xxxxx](https://doi.org/10.xxxx/xxxxx)

---

## Overview
This repository contains the R code to reproduce all analyses and figures of the above manuscript. The study is a meta-analysis examining how flower strips, landscape complexity, and climatic conditions (temperature, precipitation) jointly affect beetle and spider abundance and species richness in agricultural landscapes.

---

## Citation
**Paper**
Dietzel S, Batáry P, Twerski A, Kirmer A, Kollmann J & Fischer C (2026). Synergistic effects of flower strips and landscape complexity buffer arthropods against warming and drought – A meta-analysis. Agriculture, Ecosystems & Environment (xx). https://doi.org/xxx.yyy

**Dataset**
Dietzel S (2026). Synergistic effects of flower strips and landscape complexity buffer arthropods against warming and drought – Raw data. Zenodo. https://doi.org/10.5281/zenodo.18016074.

---

## Usage
### Analysis
The script will automatically download the data from Zenodo into `data/`:

```r
source("R/Rcode_meta_analysis_dietzel_et_al.R")
```

All figures are saved to `output/figures/` and model result tables to `output/tables/`.

---

### github repository content
.
├── data/
    └── Zenodo: 10.5281/zenodo.18016074
├── output/
│   ├── figures/
│   └── tables/
├── R/
│   └── Rcode_meta_analysis_dietzel_et_al.R
├── dietzel_et_al_2026_flower_strips_AGEE.Rproj
├── LICENSE
├── LICENSE-DATA
└── README.md