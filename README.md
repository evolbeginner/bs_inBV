# bs_inBV
A tool to help generate the file in.BV when using "complex" substitution models for MCMCTree. Particularly, the hessian is approximated by a bootstrap approach.


# Usage #
1. Perform a regular MCMCTree analysis. To save time, you can run with only the prior (e.g., set usedata = 0 in mcmctree.ctl) and very few iterations (set burnin, ).

2. sed '4!d' in.BV > ref.tre

3. ruby ~/project/Rhizobiales/scripts/dating/hessian/create_hessian_by_bootstrapping.rb --ali combined.phy --ref ref.tre --outdir dating/C60/ -b 1000 --cpu 8 -m LG+G+C60 --mcmctree_ctl mcmctree.ctl --calibrated_tree species.trees --run_mcmctree --pmsf --force 
# the above runs mcmctree (--run_mcmctree) directly after generating a C60-based in.BV. Note that w/o --pmsf it's much slower so PMSF is recommended and they seem to give highly similar time estimates (more work in validation needed though). #
