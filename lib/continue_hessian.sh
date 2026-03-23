#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# Defaults
indir="best_fit"
rscript_prog="$HOME/lab-tools/bs_inBV/latest/from_bs_to_hessian.R"

# Parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    --indir)
      indir="${2:-}"
      shift 2
      ;;
    --rscript)
      rscript_prog="${2:-}"
      shift 2
      ;;
    -h|--help)
      cat <<EOF
Usage: $0 [--indir DIR] [--rscript FILE]

Options:
  --indir DIR      Top directory containing job folders (default: best_fit)
  --rscript FILE   Path to from_bs_to_hessian.R
  -h, --help       Show this help
EOF
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

[[ -d "$indir" ]] || { echo "ERROR: --indir not found: $indir" >&2; exit 1; }
[[ -f "$rscript_prog" ]] || { echo "ERROR: R script not found: $rscript_prog" >&2; exit 1; }

for i in "$indir"/*; do
  [[ -d "$i" ]] || continue

  echo "Processing: $i"

  split_dirs=( "$i"/split/* )
  if (( ${#split_dirs[@]} == 0 )); then
    echo "  No split dirs under $i/split, skipping."
    continue
  fi

  for s in "${split_dirs[@]}"; do
    [[ -d "$s" ]] || continue

    hessian="$s/iqtree/hessian"
    bootbls="$s/iqtree/boot.bls"
    mcmc_inbv="$s/mcmctree/in.BV"

    if [[ ! -f "$bootbls" || ! -f "$mcmc_inbv" ]]; then
      echo "  Missing boot.bls or in.BV in $s, skipping this split."
      continue
    fi

    Rscript "$rscript_prog" "$bootbls" "$hessian"

    # Replace existing Hessian block
    sed -i '/Hessian/,$d' "$mcmc_inbv"
    printf "Hessian\n\n" >> "$mcmc_inbv"
    cat "$hessian" >> "$mcmc_inbv"
  done

  # Merge split/*/mcmctree/in.BV -> i/mcmctree/in.BV (natural sort)
  parts=( "$i"/split/*/mcmctree/in.BV )
  if (( ${#parts[@]} > 0 )); then
    mkdir -p "$i/mcmctree"
    printf '%s\n' "${parts[@]}" | sort -V | xargs cat > "$i/mcmctree/in.BV"
  else
    echo "  No split in.BV files found for $i"
  fi
done

