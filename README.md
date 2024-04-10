# bs_inBV

A tool with which the branch lengths, gradient, and Hessian used by `MCMCtree` to approximate the likelihood calculation are estimated with [complex amino acid substitution models](http://www.iqtree.org/doc/Complex-Models). Particularly, `bs_inBV` will return the `in.BV` file with the inferred branch lengths, gradient, and Hessian; the latter is approximated using a bootstrap approach.

## Idea

People have increasingly acknowledged the importance of using more complex amino acid substitution models in modeling sequence evolution, particularly when analysing deep phylogenies. `MCMCtree` supports various amino acid substitution models, being LG+F the most complex one. Nevertheless, more complex models are available in other phylogenetic software such as `IQ-TREE` and `PhyML` among others.

When enabling the approximate likelihood calculation with `MCMCtree` (i.e., `usedata = 2` followed by `usedata = 3`; please consult the [PAML documentation for details](https://github.com/abacus-gene/paml/blob/master/doc/pamlDOC.pdf)), `CODEML` (`BASEML` if nucleotide data) will be first "called" by `MCMCtree` when using `usedata = 2` to estimate the branch lengths, the gradient, and the Hessian under a maximum-likelihood approach; vectors and matrix that will be output in the so-called `in.BV` file. Subsequently, `MCMCtree` will use the branch lengths, gradient, and Hessian in this `in.BV` file to approximate the likelihood calculation during the MCMC when enabling `usedata = 3` in the control file. Particularly, `bs_inBV` has been written to use complex amino acid substitution models available in `IQ-TREE` to estimate the branch lengths and **approximate the Hessian matrix by calculating the negative inverse of the bootstrap covariance matrix of branch length estimates**. In addition, `bs_inBV` will then set all the elements of the gradient vector to zero.

![image](https://github.com/evolbeginner/bs_inBV/assets/8715751/6b7ae95a-f018-4331-8812-720601f637ed)

## Installation

### Ruby gems

First and foremost, please make sure that [Ruby](https://www.ruby-lang.org/en/) is installed. Second, it is recommended to have **Ruby v3.x** installed but earlier version should also work if the following gems are installed.

In addition, you will need to install the following gems before running `bs_inBV` using the command `gem install <gem_name>`. The gems to be installed are the following:

* parallel
* colorize
* find
* getoptlong
* time
* fileutils
* tmpdir
* bio
* bio-nwk

### Third-party tools

In order to run `bs_inBV`, you will need to install the following programs:

* [`PAML`](https://github.com/abacus-gene/paml). Please refer to the [PAML Wiki](https://github.com/abacus-gene/paml/wiki#installation) for more details regarding the installation procedure.
* [`IQ-TREE`](http://www.iqtree.org/). Please refer to the [IQ-TREE website](http://www.iqtree.org/#download) for more details regarding the installation procedure.
* [`newick-utils`](https://github.com/tjunier/newick_utils). If you had any issues when installing `newick-utils` from the source code, please download [the pre-compiled binaries of Newick Utilities v1.6](https://web.archive.org/web/20210409163921/http://cegg.unige.ch/pub/newick-utils-1.6-Linux-x86_64-disabled-extra.tar.gz). These binaries are available on [this website](https://web.archive.org/web/20210409163921/http://cegg.unige.ch/newick_utils).
* [`R`](https://cran.r-project.org/). Please make sure that you follow the installation instructions on the main website. WSL users will need to install R via the terminal too. In addition, please make sure that the `MASS` R package is installed as it is used by script [`reorder_node.rb`](reorder_node.rb).

### Changes to ruby scripts prior to running `bs_inBV`

At the time of writing, `bs_inBV` has pre-defined some paths an aliases to execute third-party tools in some of its ruby scripts. If you are to use this tool under a different environment, please make sure you change the following lines to match your PC settings:

* Users will need to update the path to PAML programs in line 40 in script [`do_mcmctree.rb`](lib/do_mcmctree.rb). If users have an alias to execute `MCMCtree`, they may also have to change line 43 accordingly.
* If users have not exported `IQ-TREE` to their PATH or have an alias different from the one specified in line 26 in script [`create_hessian_by_bootstrapping.rb`](create_hessian_by_bootstrapping.rb), they will need to change that too. The same applies to lines 27-28 for `newick-utils` programs and line 29 for `MCMCtree` in that same script.

## Usage

* If you have already run `MCMCtree` and have the `in.BV` file, please run the following command to obtain a "reference" tree that will be subsequently used by `bs_inBV`:

  ```sh
  sed '4!d' in.BV > ref.tre
  ```

* If you do not have an `in.BV` file, then you can run a quick dummy analysis with `MCMCtree` to obtain one. Note that the main reason for running this program now is just to make sure that the order of the tips that will be output in the `in.BV` is the same as that of the input tree file specified by the user. To save time, you can ask `MCMCtree` to sample from the prior (i.e., no data will be used, set  `usedata = 0` in [the control file](example/mcmctree.ctl)) and specify few iterations (e.g., set `burnin = 1`, `sampfreq = 1`, `nsample = 1` in [the control file](example/mcmctree.ctl)). If you are using the control file given in the [`example`](example/mcmctree.ctl) directory, please also make sure that you set `usedata = 2`. If you have your own control file, please make the changes accordingly. Once you have a "dummy" `in.BV` file, you can extract the "reference" tree:

  ```sh
  sed '4!d' in.BV > ref.tre
  ```

* Now, you are ready to run `bs_inBV`! Please execute this tool by running the following command (make sure that the paths to the third-party tools have been modified to match your settings as [aforementioned](README.md#third-party-tools)):

  ```sh
  # Run from the same directory where you have this script
  # You can use relative/absolute paths to your input files
  ruby create_hessian_by_bootstrapping.rb --ali <path_to_alignment_file> --calibrated_tree <path_to_calibrated_tree_file> --outdir <path_to_output_directory> --ref <path_to_reference_tree>/ref.tre --force -b 1000 --cpu <number_CPUs_for_IQTREE_to_use> -m LG+G+C20 --mcmctree_ctl <path_to_MCMCtree_control_file>/mcmctree.ctl --run_mcmctree --pmsf
  ```

  > Options:
  >
  > * `--ali`: path to alignment file.
  > * `--calibrated_tree`: path to input calibrated tree file that will be used by `MCMCtree`.
  > * `--outdir`: path to the output directory.
  > * `--force`: tells the program to overwrite the output files, if any.
  > * `-b`: number of bootstraps that `IQ-TREE` will run. Please note that only traditional bootstrapping is allowed as the UFB approach implemented in `IQ-TREE` cannot be used when using a fixed tree topology (which is required by `MCMCtree`).
  > * `--cpu`: number of cores that will be used by `IQ-TREE`.
  > * `-m`: amino acid substitution model that `IQ-TREE` will use for branch lengths inference.
  > * `--mcmctree_ctl`: path to the control file that will execute `MCMCtree`. Please make sure that this control file has `usedata = 2 in.BV` if you want to run `MCMCtree` after `bs_inBV` outputs the `in.BV` file (enable this feature by adding option `--run_mcmctree`).
  > * `--pmsf`: use the PMSF approximation

### Example

There is a dataset that you can use in the [`example`](example) directory to test the usage of `bs_inBV`. First, please clone this repository and follow all the steps detailed above with regards to [gems, third-party tools, and code changes during the installation procedure](README.md#installation). Once everything is ready, then go to directory [`example`](example) and execute the following command:

```sh
# Run this command from directory `example`
ruby ../create_hessian_by_bootstrapping.rb --ali combined.phy --calibrated_tree species.trees --outdir C20 --ref ref.tre --force -b 100 --cpu 8 -m LG+G+C20 --mcmctree_ctl mcmctree.ctl --run_mcmctree --pmsf
```

Once the command above finishes, you will see the following output files:

* `C20/split`: this directory will contain all the output files generated by `IQ-TREE`. The name you choose for option `--outdir` will be given to the first directory (e.g., `C20` in this example).
* `C20/mcmctree`:  if option `--run_mcmctree` is not enabled, this directory will contain your input files and the `in.BV` file. If you decide to run `MCMCtree` by enabling option `--run_mcmctree`, then you will also find the output files generated by `MCMCtree` in this directory. Please note that the `in.BV` file will contain the vector of branch length estimated using the amino acid substitution model defined by the user (e.g.,`-m LG+G+C20` in this example), the gradient vector, and the Hessian matrix approximated using a bootstrap approach.

> **NOTE 1**: Please check above to learn how to generate file `ref.tre`. For more details about how to configure the control file for `MCMCtree`, please check the [PAML documentation](https://github.com/abacus-gene/paml/blob/master/doc/pamlDOC.pdf).
> 
> **NOTE 2**: The example data were simulated under model LG+C20+G and a root age of 4.0 Ga. For more details about the simulation, please [read this blog post](https://sishuowang2022.wordpress.com/2023/06/20/substitution-models-lgcxx-and-cxxlg-differ-in-iq-tree/).

## Notes

1. When using `--run_mcmctree`, `MCMCtree` will be execute right after generating the `in.BV` file. If you only want to obtain the `in.BV` file based on your specified amino model (e.g., LG+G+C60), please do not use `--run_mcmctree`.
2. Without `--pmsf`, `IQ-TREE` is much slower because the traditional Cxx model is used by `IQ-TREE`.
3. In case you experience any issues related to inverting the variance-covariance matrix of the branch lengths, please try to (i) increase the number of bootstrap or (ii) avoid including too closely-related species in your phylogeny. However, in the current version (actually since v1.0) the issue of singular covariance matrix is solved by creating a new matrix by adding a small number of diagonal term to the original covariance matrix, so probably no need to worry about it.
4. The tool should **in theory work for any substitution models** but would be **in practice most useful for complex models** as these are not included in `MCMCtree`. It also works for nucleotides in which case however you may need to specify seqtype=1 in the file `mcmctree.ctl`.
5. For more details, please see **Note S4 and Figs. S9-S10** in the original reference (see below).

## How to cite

Dating the bacterial tree of life based on ancient symbiosis Sishuo Wang, Haiwei Luo bioRxiv 2023.06.18.545440; doi: https://doi.org/10.1101/2023.06.18.545440.

You may also need to cite corresponding papers for the use of `PAML`, `IQ-TREE`, and `Newick Utilities`.

