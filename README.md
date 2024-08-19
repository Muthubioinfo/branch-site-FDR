# Statistical properties of branch-site test and FDR

This repository contains the files to recreate the analysis in Chapter 2 - Statistical properties of the branch-site test and false discovery rate to detect episodic positive selection.

To run the analysis, users are required to download and compile ```EVOLVER``` and ```CODEML``` programs available in ```PAML``` package, see [http://abacus.gene.ucl.ac.uk/software/paml.html](http://abacus.gene.ucl.ac.uk/software/paml.html). Use bash/LINUX commands to run PAML programs.

## Simulating codon alignments under different levels of positive selection
To simulate alignments codon alignments, use the ```EVOLVER``` program by setting the option to 6 (for codon datasets) and specifying the dat file ```MCcodonNSbranchsites.dat```

```
evolver 6 MCcodonNSbranchsites.dat
```

Phylogenomic dataset is simulated for two trees names as TREE I and TREE II. TREE I consists of 8 species and TREE II has 16 species. Both the trees assume uniform codon frequencies (1/61) and uniform branch lengths of 0.3 nucleotide substitution per site. The total number of simulations are set to 10,000 genes/hypothesis. 

Under this basis of simulation setup, we can different phylogenomic datasets under different levels of positive selection, determined by the parameters of the branch-site test.
These are parameters of the branch-site test that determine positive selection are -
| Parameters of branch-site test |Description                    | Values used for different simulation   |
|--------------------------------|--------------------------------|----------------------------------------|
|        $\omega{_2}$           | Strength of positive selection | 1.5, 2, 3, 5, 8, 10 |
| $p_{2}$                        | Proportion of positively selected sites | 0.01, 0.05, 0.15, 0.20, 0.50 |
| $L_{c}$                     | Number of codon sites | 100, 200, 1000, 10000 | 
| %TP                          | Percentage true positives among all simulations | 0.1, 0.2, 0.5, 1, 2, 5, 20, 50, 100 |


A total of 24 different phylogenomic datasets can be simulated.  The positive selection for the two trees (TREE I has 8 species and TREE II has 16 species) are tested by labelling either internal ($\alpha$) and external branch ($\beta$).

Use the ```MCcodonNSbranchsites.dat``` control files from [MCbranchsite_files_simulations](https://github.com/Muthubioinfo/branch-site_FDR/tree/main/MCbranchsite_files_simulations) directory to simulate codon alignments under different positive selection setting, as shown in the above table. Each subdirectory named ```alpha_branch_...``` and ```beta_branch_...``` contains the control file for each parametric setting of branch-site test (omega, codon sites, ) labelling either alpha (internal) or beta (external) as foreground branches. For example, the ```MCcodonNSbranchsites.dat``` file from ```alpha_branch_omega``` assumes four site classes (0, 1, 2a and 2b). In 2a and 2b, the $\omega{_2} > 1$ at the $\alpha$ branch by specifying omega = 1.5, 2, 3, 5, 8 and 10. 

In the ```MCcodonNSbranchsites.dat```, change the number of species to 16 to simulate codon alignments under TREE II.


## Branch-site test for positive selection

### Log-likelihood estimation
For the given sequence alignment, the log-likelihoods are evaluated under the null model ($\ell_{0}$) and alternative model ($\ell_{1}$) using the ```CODEML``` program. Under the null model, the log-likelihood ($\ell_{0}$) is estimated by fixing the [$\omega$ to 1](https://github.com/Muthubioinfo/branch-site_FDR/tree/main/codeml_ctl_positive_selection/null_model), whereas the omega and other free parameters are estimated under the alternative model. Use the control files from [codeml_ctl_positive_selection](https://github.com/Muthubioinfo/branch-site_FDR/tree/main/codeml_ctl_positive_selection) to evalute the $\ell_{0}$ and $ell_{1}$ under the null and alternative models respectively. The log-likelihood under alternative model is estimated for three iterations each with initial $\omega = [1](https://github.com/Muthubioinfo/branch-site_FDR/tree/main/codeml_ctl_positive_selection/alternate_iteration_1), [2](https://github.com/Muthubioinfo/branch-site_FDR/tree/main/codeml_ctl_positive_selection/alternate_iteration_2),$ and $[4](https://github.com/Muthubioinfo/branch-site_FDR/tree/main/codeml_ctl_positive_selection/alternate_iteration_3)$. 

Running the ```CODEML```  program returns the output in ```rst1``` file containing the estimated parameters and maximum likelihood for each gene. We generate four rst1 files for the four iterations (null model, alternative model 1, 2 and 3). The rst1 files after analysing under null hypothesis is named as ```rst1H0null```, ```rst1H0alt1```,```rst1H0alt2``` and ```rst1H0alt3```. Similarly, the rst1 files for simulations under the alternative hypothesis are named as ```rst1H1null```, ```rst1H1alt1```,```rst1H1alt2``` and ```rst1H1alt3```. 


### Likelihood ratio test, P-values and FDR
Import all the ```rst1``` files in R/Rstudio. Run the ```simFDR()``` function in  available in [```function_FDR.R```](https://github.com/Muthubioinfo/branch-site_FDR/blob/main/function_FDR.R), and calculate the likelihood ratio test or ```LRT```, and ```P-values```, and the ```q-values``` under BH-FDR (Benjamini and Hochberg, 1995) and ST-FDR (Storey, 2002) correction method. The P-values and q-values are significant at 5%. In the ```FDR``` function calculates the power of positive selection using both the BH-FDR (Benjamini-Hochberg, 1995) and ST-FDR (Storey, 2002) methods. 


## Realistic simulation using empirical data
For this section, I test the statistical properties of branch-site test and FDR under real data setting. The dataset is obtained from [Kosiol et al. (2008)](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000144), consisting of 9,566 genes with ancestral primate branch focussed as the foreground. Download the sequence alignments from [http://compgen.cshl.edu/projects/mammal-psg/](http://compgen.cshl.edu/projects/mammal-psg/).
The tree topology specified in the Kosiol dataset is fixed for all the genes. Evaluate the log-likelihood under the null ($\ell_{0}$) and alternative ($\ell_{1}$) model assuming F3X4 codon frequency model. The respective data of the parameter estimates such as branch lengths, kappa, are to be stored in a separate file (e.g. either in ```.xls```, ```.Rdata```, ```.csv``` or ```.txt```). 

The aim is to construct a realistic simulation using the existing primate dataset. For this analysis, I recommend to remove all the genes that has been inferred under positive selection such that the P-values evaluated are not significant, and ensure that all the alignments do not have any missing data (or missing species) within the specified tree topology. This ensures to avoid any unwanted noise in the experiment. 
For the primate dataset, the filtering resulted in 6903 genes. These genes are considered as neutral genes for the primate branch.

To construct a realistic simulation, randomly select $n = 500$ genes. You can use the following command to move files to another directory.

```
find . -mindepth 1 -maxdepth 1 -type f | shuf | head -n 10 | xargs -I{} mv {} dest_dir/modified/
```

Note: 'dest_dir' is the directory with all the 6903 genes or gene alignments. Create a directory called 'modified' before running the above command. It is best to keep a copy of the filenames of the randomly generated $n$ genes to avoid any mishaps. 

The respective parameters estimated from $n$ genes (such as branch lengths and transition-transversion ratio) are used in ```EVOLVER``` for simulating positive selected codons. For simulating positive selection, one uses alternative hypothesis with either $\omega_2 = 4$ (moderate selective pressure) or $\omega_2 = 10$ (strong selective pressure). Then, these newly simulated positive codons are concatenated at the end of $n$ gene alignment. Thus, the modified alignments now has the positively selected codons and are assumed to be under the alternative hypothesis. The remaining 6,403 unmodified genes are assumed to be under the null hypothesis.

Now, analyse all the alignments under the branch-site test. Use the ```function_FDR.R```, to calculate the LRT, and the P-values and q-value significant at 5% confidence level. 

## Simulation under misspecification of the null model






