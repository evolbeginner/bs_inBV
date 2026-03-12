### v1.4.3 — 2026-03-12
- **New:** Improved `--param_bs` (or `--bs pbs`), `--best_fit [y]`

### v1.4.2 — 2026-03-12
- **New:** Improved `--param_bs`. Now parametric bs is better performed and particularly the script `lib/pmsf_sitewise_alisim.rb` helps simulating alignment under pmsf.

### v1.4.1 — 2026-03-10
- **New:** Added `--full_pmsf`.

### v1.4.0 — 2026-03-10
- **New:** Added full support for `-b` with `--pmsf` by the script `lib/bs_phylip_noallgap.rb`.
- **Improved:** PMSF under `--best_fit`

### v1.3.4 — 2026-03-07
- **Fixed:** `--best_fit` for pmsf by `-fs` and `-keep-ident`

### v1.3.3 — 2026-03-06
- **Fixed:** `--best_fit`

### v1.3.2 — 2026-03-06
- **Improved:** progressbar improved
- **Fixed:** `--best_fit` for IQ-Tree v2

### v1.3 — 2026-01-29
- **New:** Added support for bootstrapping with the best-fitting subs. model if `--best_fit` is specified.
- **Improved:** The organization of `create_hessian_by_bootstrapping.rb`.
- **Fixed:** N/A.
