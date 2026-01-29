# bs_inBV

A tool with which the branch lengths, gradient, and Hessian used by `MCMCtree` to approximate the likelihood calculation are estimated with [complex amino acid substitution models](http://www.iqtree.org/doc/Complex-Models).  
In particular, `bs_inBV` returns the `in.BV` file containing inferred branch lengths, a gradient vector, and a Hessian matrix, where the Hessian is approximated using a bootstrap approach.

---

## Idea

The use of complex mixture substitution models has become increasingly important in phylogenetic analyses, especially for deep evolutionary timescales. While `MCMCtree` supports many amino acid substitution models, many widely used complex models, such as +I, +R, EX2, C10–C60, their combinations, and branch-length mixture models like GHOST, are not available. By contrast, most if not all of such models are implemented in `IQ-TREE`, among other phylogenetic software.

When approximate likelihood calculation is enabled in `MCMCtree` (i.e. `usedata = 2` followed by `usedata = 3`; see the [PAML documentation](https://github.com/abacus-gene/paml/blob/master/doc/pamlDOC.pdf)), `CODEML` (or `BASEML` for nucleotide data) is first invoked to estimate branch lengths, the gradient vector, and the Hessian matrix under a maximum-likelihood framework. These quantities are stored in the `in.BV` file.

`bs_inBV` replaces this step by:
- estimating branch lengths under complex substitution models using `IQ-TREE`,
- approximating the Hessian matrix as the **negative inverse of the bootstrap covariance matrix of branch-length estimates**, and
- setting all elements of the gradient vector to zero.

This enables the use of complex substitution models during divergence-time estimation in `MCMCtree`.

![image](https://github.com/evolbeginner/bs_inBV/assets/8715751/6b7ae95a-f018-4331-8812-720601f637ed)

---

## Installation

### Ruby gems

Please ensure that [Ruby](https://www.ruby-lang.org/en/) is installed. Ruby **v3.x** is recommended, but earlier versions should work if the required gems are available.

Install the following gems:

```sh
gem install parallel colorize find getoptlong time fileutils tmpdir bio bio-nwk
```

---

### Third-party tools

The following external programs must be installed and accessible:

- [`PAML`](https://github.com/abacus-gene/paml)  
  See the [PAML Wiki](https://github.com/abacus-gene/paml/wiki#installation).

- [`IQ-TREE`](http://www.iqtree.org/)  
  Installation instructions are available on the [IQ-TREE website](http://www.iqtree.org/#download).

- [`newick-utils`](https://github.com/tjunier/newick_utils)  
  If compiling from source leads to errors, download the pre-compiled binaries of Newick Utilities v1.6 from the [link] (https://web.archive.org/web/20210409163921/http://cegg.unige.ch/pub/newick-utils-1.6-Linux-x86_64-disabled-extra.tar.gz) or use bioconda.

- [`R`](https://cran.r-project.org/)  
  Ensure that the `MASS` package is installed. This is required by `reorder_node.rb`.

---

### Required code adjustments

Before running `bs_inBV`, you may need to edit paths in the following scripts:

- `lib/do_mcmctree.rb`  
  Update the path to `MCMCtree` (line ~40) if it is not available in your `PATH`.

- `create_hessian_by_bootstrapping.rb`  
  Update paths or aliases for `IQ-TREE`, `newick-utils`, and `MCMCtree` if necessary (lines ~26–29).

---

## Usage

### Basic command

```sh
ruby create_hessian_by_bootstrapping.rb \
  --ali <alignment.phy> \
  --calibrated_tree <species.trees> \
  --ref <ref.tre> \
  --outdir <output_directory> \
  [OPTIONS]
```

---

### Required arguments

| Option | Description |
|------|-------------|
| `--ali` | Path to the input alignment file (PHYLIP format). |
| `--calibrated_tree` | Path to the calibrated species tree used by `MCMCtree`. |
| `--ref`, `--ref_tree_file` | Path to the reference tree (`ref.tre`) extracted from an existing `in.BV` file. |
| `--outdir` | Output directory where all results will be written. |

---

### Optional arguments

| Option | Default | Description |
|------|--------|-------------|
| `-m` | `LG+G` | Amino acid substitution model used by `IQ-TREE`. |
| `-b` | `1000` | Number of standard bootstrap replicates. |
| `--cpu` | `1` | Number of CPU cores used by `IQ-TREE`. |
| `--pmsf` | off | Use PMSF approximation (recommended for complex mixture models). |
| `--best_fit`, `--best-fit` | off | Select the best-fitting substitution model before bootstrapping. |
| `--force` | off | Overwrite the output directory if it already exists. |
| `--mcmctree_ctl` | Path to the `MCMCtree` control file (`.ctl`). |
| `--run_mcmctree` | off | Run `MCMCtree` immediately after generating `in.BV`. |
| `--no_mcmctree` | off | Do not run `MCMCtree`. |
| `--te` | none | Path to a fixed tree topology used by `IQ-TREE`. |
| `--add_cmd`, `--add_argu` | `-mwopt` | Additional arguments passed to `IQ-TREE`. |
| `--no_mwopt` | off | Disable `-mwopt` in `IQ-TREE` arguments. |

---

### Notes on key options

- **Bootstrap (`-b`)**  
  Only standard bootstrap is supported. Ultrafast bootstrap (UFBoot) cannot be used because `MCMCtree` requires a fixed topology.

- **`--pmsf`**  
  Strongly recommended for mixture models (e.g. C10–C60) to improve computational efficiency.

- **`--best_fit`**  
  Ensures that the bootstrap-based Hessian is consistent with the approximate-likelihood setting used by `MCMCtree`.

- **`--run_mcmctree`**  
  If enabled, `MCMCtree` is executed automatically after the `in.BV` file is generated.  
  Omit this option if you only want the `in.BV` file.

---

### Example

```sh
ruby create_hessian_by_bootstrapping.rb \
  --ali combined.phy \
  --calibrated_tree species.trees \
  --ref ref.tre \
  --outdir C20 \
  --mcmctree_ctl mcmctree.ctl \
  -m LG+G+C20 \
  -b 100 \
  --cpu 8 \
  --pmsf \
  --run_mcmctree \
  --force
```

---

## Notes

1. When `--run_mcmctree` is enabled, `MCMCtree` is executed immediately after `in.BV` is generated.
2. Without `--pmsf`, `IQ-TREE` may be substantially slower for complex mixture models.
3. Potential issues with singular covariance matrices are handled internally by adding a small diagonal regularization term.
4. The tool is applicable in principle to any substitution model supported by `IQ-TREE`. For nucleotide data, set `seqtype = 1` in `mcmctree.ctl`.
5. For methodological details, see **Note S4 and Figs. S9–S10** in the reference below.

---

## Acknowledgement

I am particularly grateful to Sandra Álvarez-Carretero and Edmund Moody (University of Bristol) for testing the scripts and providing valuable feedback. Sandra also contributed to code development and improvements to the documentation.

---

## How to cite

Wang, Sishuo, and Haiwei Luo. "Dating the bacterial tree of life based on ancient symbiosis." Systematic Biology (2025): syae071.

Please also cite the relevant publications for `PAML`, `IQ-TREE`, `Newick Utilities`, and `BioRuby`.
```

