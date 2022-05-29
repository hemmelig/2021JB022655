# Supplementary analysis/visualisation code
Data and code required to reproduce the visualisations presented in:

Bacon, C.A., Johnson, J.H., White, R.S. and Rawlinson, N., 2022. On the origin of seismic anisotropy in the shallow crust of the Northern Volcanic Zone, Iceland. Journal of Geophysical Research: Solid Earth, 127(1), p.e2021JB022655.

Article: [![DOI](https://img.shields.io/badge/DOI-10.1029/2021JB022655-blue)](https://doi.org/10.1029/2021JB022655)

Supplementary datafiles: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5007022.svg)](https://doi.org/10.5281/zenodo.5007022)

Analysis and visualisation: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5636924.svg)](https://doi.org/10.5281/zenodo.5636924)

## Steps to reproduce figures
1. Clone this repository and navigate to it, e.g.:

```
git clone https://github.com/hemmelig/2021JB022655
cd 2021JB022655
```

2. Install the packages listed in the environment.yml file, either manually, or using (for example) conda:

```
conda env create
```

3. Add `mfast_summary.mplstyle` (a matplotlib stylesheet) to your `mpl_configdir` (usually found at `~/.config/matplotlib`)

4. Optional: Install Helvetica font for Matplotlib

5. Navigate to each figure directory and run the `.gmt` (as `bash <script>.gmt`) or `.py` (as `python <script>.py`) scripts.

## Notes
These figures were prepared using Linux 20.04.
