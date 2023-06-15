# bs_inBV
A tool to help generate the file in.BV when using "complex" substitution models for MCMCTree. Particularly, the hessian is approximated by a bootstrap approach.

# Idea #
People have increasingly realized the importance of using more complex substitution models in modeling sequence evolution, particularly for deep-time evolution and in case the sequence divergence is large. MCMCTree supports only a number of substitution models. These models, however, may be available in other software specifically used for tree construction, e.g., IQ-Tree and RAxML.

To work with the approximate likelihood method of MCMCTree (usedata=2), one needs to get the MLE (maximum likelihood estimate) of the branch length, as well as the gradient and hessian evaluated at the MLE. This tool takes the advantages of the abundant substitution models of IQ-Tree. Briefly, the tool uses IQ-Tree's estimate of the branch length, set the gradient to all zeros, and importantly, #approximate the hessian by calculating the negative inverse of the bootstrap covariance matrix of branch length estimates#.

# Installation
Make sure [RUBY] (https://www.ruby-lang.org/en/) is installed.

Many scripts included in this product require the following RUBY packages. If any is not installed, please install it by `gem install package_name`.
* [parallel]
* [colorize]

This product also employs several other computational tools. Please ensure that you have them installed.
* [PAML](https://github.com/abacus-gene/paml)
* [IQ-Tree](http://www.iqtree.org/)

# Usage #
1. Perform a regular MCMCTree analysis. To save time, you can run with only the prior (e.g., set usedata = 0 in mcmctree.ctl) and very few iterations (set burnin=1, sampfreq=1, nsample=1 in mcmctree.ctl).
  The reason for this step is to make sure the order of the tips in the file in.BV generated by MCMCTree is the same as that of the species.tree specified by the user.

2. Extract a "reference" tree from the file in.BV generated by MCMCTree in the previous step.
  
    `sed '4!d' in.BV > ref.tre`

3. Run the following

    `ruby create_hessian_by_bootstrapping.rb --ali alignment --calibrated_tree species_tree --outdir C20 --ref ref.tre --force -b 1000 --cpu 8 -m LG+G+C20 --mcmctree_ctl mcmctree.ctl --run_mcmctree --pmsf`
    
    Arguments:
      * `--ali`: alignment file
      * `--calibrated_tree`: the species tree used in regular MCMCTree
      * `--outdir`: output directory
      * `--force`: tells the program to overwrite the output if it exists
      * `-b`: the no. of bootstraps (note that only traditional bootstrapping is allowed as IQ-TRee's UFB cannot be used for a species tree with a fixed topology)
      * `--cpu`: no. of cores used in IQ-Tree
      * `-m`: the model passed to IQ-Tree
      * `--mcmctree_ctl`: the control file for MCMCTree
      * `--pmsf`: use the PMSF approximation

# Notes #
1. With `--run_mcmctree`, MCMCTree will be run directly after generating in.BV. In case you want only the file in.BV based on your specified model say LG+G+C60, please do not use `--run_mcmctree`.
2. Without `--pmsf` IQ-Tree is much slower because it is the traditional Cxx model that will be applied in IQ-Tree.

# How to cite
You may also need to cite corresponding papers for the use of PAML and IQ-Tree.

