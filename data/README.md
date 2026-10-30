# data/

This folder holds the scripts that download the raw data. The raw files and cleaned files are not committed. See `../DATA.md` for the full data record.

| Item | What it is |
|---|---|
| `download_lalonde.R` | Downloads the NSW (Dehejia-Wahba), PSID and CPS files from https://users.nber.org/~rdehejia/nswdata2.html into `data/raw/`, checks row counts, and writes `data/raw/manifest.csv`. |
| `raw/` | Created by the download script. Listed in `.gitignore`. |
| `processed/` | Created by `../src/00_prep_data.R`. Listed in `.gitignore`. |

Terms: "distributed for attributable non-commercial use (CC by NC)", as stated on the source page. Citations are listed in `../DATA.md`, Section 2.

The Hillstrom e-mail file is not downloaded. Its terms could not be confirmed (see `../DESIGN.md`, Section 9).

Run from the repository root:

```
Rscript data/download_lalonde.R
```
