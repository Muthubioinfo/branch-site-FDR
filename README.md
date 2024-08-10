# Statistical properties of branch-site test and FDR

This repository contains the files to recreate the analysis in Chapter 2 - Statistical properties of the branch-site test and false discovery rate to detect episodic positive selection.

To run the analysis, users are required to download and compile ```EVOLVER``` and ```CODEML``` programs available in ```PAML``` package, see [http://abacus.gene.ucl.ac.uk/software/paml.html](http://abacus.gene.ucl.ac.uk/software/paml.html). Use bash/LINUX commands to run PAML programs.

## Simulating codon alignments under different levels of positive selection
To simulate alignments codon alignments, use the ```EVOLVER``` program. Two trees consisting of 8 species and 16 species with uniform codon frequencies (1/61) and uniform branch lengths of 0.3 nucleotide substitution per site.

The simulated trees are constructed under different parameteric setting of the branch-site test that determine the level of positive selection in a codon alignment.  
These are parameters of the branch-site test that determine positive selection are -
| Parameters of branch-site test |Description                    | Values used for different simulation   |
|--------------------------------|--------------------------------|----------------------------------------|
|        $\omega{_2}$            | Strength of positive selection | 1.5, 2, 3, 5, 8, 10 |
| $p_{2}$                        | Proportion of positively selected sites | 0.01, 0.05, 0.15, 0.20, 0.50 |
| $L_{c}$                        | Number of codon sites | 100, 200, 1000, 10000 | 
| %TP                            | Percentage true positives among all simulations | 0.1, 0.2, 0.5, 1, 2, 5, 20, 50, 100 |

Thus, a total of 24 different positive selection setting are achieved.  The positive selection for the two trees (TREE I - 8 species; TREE II - 16 species) are tested by labelling either internal ($\alpha$) and external branch ($\beta$).

Use the ```MCcodonNSbranchsites.dat``` control file from [MCbranchsite_files_simulations](https://github.com/Muthubioinfo/branch-site_FDR/tree/main/MCbranchsite_files_simulations) directory to simulate codon alignments under different positive selection setting. In this directory, the folders ```alpha_branch_...``` and ```beta_branch_...``` has the control files for simulating sequences with positive selection at either alpha (internal) and beta (external) branches. For example, ```alpha_branch_omega``` has the ```MCcodonNSbranchsites.dat``` that has site classes labelled at alpha branch by specifying omega = 1.5, 2, 3, 5, 8 and 10.


In the ```MCcodonNSbranchsites.dat```, change the number of species to 16 to simulate codon alignments under TREE II. When the null hypothesis is true, all the site class of the branch-site test have conserved sites (i.e., site class 0 with proportion of sites, $p_{0}$ having $0 < \omega < 1$) or neutral (Site class 1 with proportion of sites, $p_{1}$ having $\omega = 1$). When the codon alignment is simulated under the alternate hypothesis, the site class 2a and 2b has foreground $\omega_{2} = 4$. 


## Branch-site test for positive selection

### Log-likelihood estimation
For the given sequence alignment, the log-likelihoods are evaluated under the null model ($\ell_{0}$) and alternative model ($\ell_{1}$). Under the null model, the log-likelihood ($\ell_{0}$) is estimated by fixing the omega to 1, whereas the omega and other free parameters are estimated under the alternative model. Use the control files in [codeml_ctl_positive_selection](https://github.com/Muthubioinfo/branch-site_FDR/tree/main/codeml_ctl_positive_selection) to evalute the $\ell_{0}$ and $ell_{1}$ under the null and alternative models respectively.


The log-likelihood under alternative model is estimated for three iterations, each initialised with $\omega = 1, 2,$ and $4$. 
Then, take the largest log-likelihood among the all log-likelihood estimates under both the null and alternative models. This is to ensure that the $\ell_{1} \le \ell_{0}$).

### Likelihood ratio test, P-values and FDR
The ```rst1``` file output from ```CODEML``` are imported in R. Run the ```FDR``` function in ```function_FDR.R``` and calculate the likelihood ratio test or ```LRT``` and ```P-values```, and the ```q-values``` after FDR correction method. In the ```FDR``` function calculates the power of postitive selection using both the BH-FDR (Benjamini-Hochberg, 1995) and ST-FDR (Storey, 2002) methods. 

## Real dataset

## Simulating model misspecification






