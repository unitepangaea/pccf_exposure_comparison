# Data

Nothing in this folder is tracked by git.

- `derived/` is written by the scripts of this repository (`DERIVED` in `config.R`).
- Inputs are read from the pccf_nar_validation checkout (`VALIDATION` in `config.R`):
  the public raw data (`RAW`: dissemination area boundaries, MSDI 2021, MODIS EVI), the
  licensed PCCF+ 8D files (`PCCF_DIR`) and the outputs of its pipeline (`VAL_DERIVED`:
  `linked/`, `p2/`, `p4/`). See the README of
  https://github.com/unitepangaea/pccf_nar_validation for how to obtain them.

PCCF+ is licensed by Statistics Canada. Neither the PCCF+ files nor any address-,
postal-code- or DA-level result derived from them may be placed in a public repository.
