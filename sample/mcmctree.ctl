          seed = -1
       seqfile = combined.phy
      treefile = species.trees
      mcmcfile = mcmc.txt
       outfile = out.txt

         ndata = 1        * Number of partitions
       seqtype = 2        * 0: nucleotides; 1:codons; 2:AAs
       usedata = 2 in.BV  * 0: no data (prior); 1:exact likelihood;
                          * 2:approximate likelihood; 3:out.BV (in.BV)
         clock = 3        * 1: global clock; 2: independent rates; 3: correlated rates 

         model = 3        * models for AAs or codon-translated AAs:
                          *     0:poisson, 1:proportional,2:Empirical,3:Empirical+F
                          *     6:FromCodon, 8:REVaa_0, 9:REVaa(nr=189)         alpha = 0.5
         ncatG = 4        * No. categories in discrete gamma
         alpha = 0.5      * alpha for gamma rates at sites

       BDparas = 1 1 0.1   * birth, death, sampling
   kappa_gamma = 6 2       * gamma prior for kappa
   alpha_gamma = 1 1       * gamma prior for alpha

   rgene_gamma = 1 50      * gammaDir prior for rate for genes
  sigma2_gamma = 1 10      * gammaDir prior for sigma^2     (for clock=2 or 3)

         print = 1         * 0: no mcmc sample; 1: everything except branch rates 2: everything
        burnin = 100000
      sampfreq = 1000 
       nsample = 20000
