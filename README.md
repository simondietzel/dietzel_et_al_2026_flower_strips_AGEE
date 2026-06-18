# Data and Code for Dietzel et al., 2026

## Publication
Synergistic effects of flower strips and landscape complexity buffer arthropods against warming and drought – A meta-analysis. Agriculture, Ecosystems & Environment. https://doi.org/10.xxxx/xxxxx](https://doi.org/10.xxxx/xxxxx.

Simon Dietzel <a href="https://orcid.org/0000-0003-4319-8195"><img src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16"/></a>
Péter Batáry <a href="https://orcid.org/0000-0002-1017-6996"><img src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16"/></a>
Alina Twerski <a href="https://orcid.org/0000-0001-7966-1335"><img src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16"/></a>
Anita Kirmer <a href="https://orcid.org/0000-0002-2396-713X"><img src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16"/></a>
Johannes Kollmann <a href="https://orcid.org/0000-0002-4990-3636"><img src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16"/></a>
Christina Fischer <a href="https://orcid.org/0000-0001-7784-1105"><img src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16"/></a>

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
```r
├── data
│   └── Zenodo: https://doi.org/10.5281/zenodo.18016074
├── img
│   ├── beetle_img.png
│   └── spider_img.png
├── output
│   ├── figures
│   └── tables
├── R
│   └── Rcode_meta_analysis_dietzel_et_al.R
├── renv
│   ├── activate.R
│   ├── library
├── dietzel_et_al_2026_flower_strips_AGEE.Rproj
├── README.md
├── LICENSE
├── LICENSE-DATA
└── renv.lock
```

---

### Copyright
**Code**: Published under MIT License (LICENSE)

**Data**: Published under CC-BY 4.0 (LICENSE-DATA)