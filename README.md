# Statistical properties of branch-site test and FDR

This repository contains the files to recreate the analysis in Chapter 2 - Statistical properties of the branch-site test and false discovery rate to detect episodic positive selection.

To run the analysis, you are required to download and compile ```EVOLVER``` and ```CODEML``` programs available in ```PAML``` package, see [http://abacus.gene.ucl.ac.uk/software/paml.html](http://abacus.gene.ucl.ac.uk/software/paml.html). Use bash/LINUX commands to run PAML programs.

## Simulating codon alignments under different levels of positive selection
Consider two phylogenetic trees, TREE I with 8 species and TREE II with 16 species. Both trees assume uniform codon frequencies (1/61) and a uniform branch length of 0.3 nucleotide substitutions per site. The transition-transversion ratio $\kappa$ is set to 2. We simulate TREE I and TREE II in two settings, each labelled at (1) alpha branch, (2) beta branch. Thus, we have four different tree scenarios, as shown below

![trees]([https://github.com/Muthubioinfo/branch-site_FDR/tree/main/branch-site-test/trees.pdf))



Simulate 10,000 protein-coding genes or codon alignments, with 90% genes representing null hypotheses (9000 genes) and the remaining 10% represent the alternative hypotheses, making up a phylogenome of 10,000 genes. 

We can improve the level of phylogenomic simulation by simulating under different levels of positive selection, as determined by the parameters of the branch-site model. 


|Parameters of branch-site test | Description | Parameter value at each phylogenomic simulation  |
|----------------------------------|-------------------------------------------------|----------------|
|        $\omega{_2}$           | Strength of positive selection | 1.5, 2, 3, 5, 8, 10 |
| $p_{2}$                        | Proportion of positively selected sites | 0.01, 0.05, 0.15, 0.20, 0.50 |
| $L_{c}$                     | Number of codon sites | 100, 200, 1000, 10000 |
| %TP                          | Percentage true positives among all simulations | 0.1, 0.2, 0.5, 1, 2, 5, 20, 50, 100 |


Now, a total of 24 different phylogenomic datasets can be simulated based on the table above.  The positive selection for the two trees (TREE I has 8 species and TREE II has 16 species) are tested by labelling either internal ($\alpha$) and external branch ($\beta$).

To simulate codon sequence alignments under different settings shown above, use the EVOLVER program and select option 6. 

```
evolver 6 MCcodonNSbranchsites.dat
```

Ensure that the appropriate control file called ```MCcodonNSbranchsites.dat``` is specified in the working directory. The control file for the respective tree (tree 1 and tree 2), the foreground branch (alpha and beta) and the parametric setting ($L_{c}$, $\omega$, $p_{2}$, %TP) used in phylogenomic simulations is in [MCbranchsite_files_simulations](https://github.com/Muthubioinfo/branch-site_FDR/tree/main/MCbranchsite_files_simulations). 

## Branch-site test for positive selection

### Log-likelihood estimation
For each simulated alignment, evaluated the log-likelihoods under the null model ($\ell_{0}$) and alternative model ($\ell_{1}$) using the ```CODEML``` program. Under the null model, the log-likelihood ($\ell_{0}$) is estimated by fixing the $\omega$ to [1](https://github.com/Muthubioinfo/branch-site_FDR/tree/main/codeml_ctl_positive_selection/null_model), whereas the omega and other free parameters are estimated under the alternative model. 
Use the control files from [codeml_ctl_positive_selection](https://github.com/Muthubioinfo/branch-site_FDR/tree/main/codeml_ctl_positive_selection) to evaluate the $\ell_{0}$ and $\ell_{1}$ under the null and alternative models respectively. The log-likelihood under alternative model is estimated for three iterations each with initial $\omega$ = [1](https://github.com/Muthubioinfo/branch-site_FDR/tree/main/codeml_ctl_positive_selection/alternate_iteration_1), [2](https://github.com/Muthubioinfo/branch-site_FDR/tree/main/codeml_ctl_positive_selection/alternate_iteration_2), and [4](https://github.com/Muthubioinfo/branch-site_FDR/tree/main/codeml_ctl_positive_selection/alternate_iteration_3). 

The ```CODEML``` returns the output file called ```rst1``` file, which contains the estimated parameters and maximum likelihood for each gene. We generate four rst1 files for the four codeml iterations (null model, alternative model 1, 2 and 3). The rst1 files after analysing under null hypothesis is named as ```rst1H0null```, ```rst1H0alt1```,```rst1H0alt2``` and ```rst1H0alt3```. Similarly, the rst1 files for simulations under the alternative hypothesis are named as ```rst1H1null```, ```rst1H1alt1```,```rst1H1alt2``` and ```rst1H1alt3```. 

### Likelihood ratio test, P-values and FDR
Import all the ```rst1``` files in R/Rstudio and apply the ```simFDR()``` function available in [```simFDR.R```](https://github.com/Muthubioinfo/branch-site_FDR/blob/main/function_FDR.R), and calculate the likelihood ratio test or ```LRT```, and ```P-values```, and the ```q-values``` under BH-FDR (Benjamini and Hochberg, 1995) and ST-FDR (Storey, 2002) correction method. Here, the P-values and q-values are significant at 5%. The classical power of branch-site test is compared 

In the ```simFDR``` function calculates the power of positive selection using both the BH-FDR (Benjamini-Hochberg, 1995) and ST-FDR (Storey, 2002) methods. 

## Real data analysis
For this section, I test the statistical properties of branch-site test and FDR under real data setting. The dataset is obtained from [Kosiol et al. (2008)](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000144), consisting of 9,566 genes with ancestral primate branch focussed as the foreground. 

The primate dataset from Kosiol et al. (2008) is available in [kosiol_primate_dataset.txt.zip](https://github.com/Muthubioinfo/branch-site_FDR/tree/main/real_data_files) to work out the experiment in this section. Use ```split``` the file into 9566 files each containing the 9566 codon alignment.

```
awk -v RS= '{print > ("s" NR ".txt")}' kosiol_primate_dataset.txt
```
Note: The above code splits the ```kosiol_primate_dataset.txt``` file into 9,566 files each containing a codon sequence alignment. The "s" in the above code generates files as ```s1.txt, s2.txt, s3.txt,......,s9566.txt```. 

The tree topology specified in the Kosiol dataset is fixed for all the genes. The control files for the four iterations of CODEML is available in [real_data_files](https://github.com/Muthubioinfo/branch-site_FDR/tree/main/real_data_files). Evaluate the log-likelihood under the null ($\ell_{0}$) and alternative ($\ell_{1}$) model assuming F3X4 codon frequency model. 

```
#Loop to link all files to the target directory for branch-site test
for i in {1..9566}
cd /target-directory/
mkdir s$i

ln -s /location-of-sequence/s$i.txt /target-directory/s$i/seq.txt
ln -s /location-phylip-tree/realdata.trees /target-directory/tree.trees
ln -s /location-CODEML-control-file/codeml.ctl /target-directory/codeml.ctl
cp /compiled-CODEML-program/codeml /target-directory/

#Execute the program
codeml codeml.ctl

cd ..
done
```

The above script is applied four times, for estimating log-likelihoods under null model, alternative model 1, alternative model 2, and alternative model 3. The respective data of the parameter estimates such as branch lengths, kappa, and the codon frequencies, are stored in a separate file (e.g. either in ```.xls```, ```.Rdata```, ```.csv``` or ```.txt```). You can use R to output all these from the ```mlc``` and ```rst1``` output files.


## Realistic simulation using empirical data
The aim is to construct a realistic simulation using the above primate dataset. 

Step 1:
For this analysis, the first step is to remove all the positively selected genes inferred under positive selection. This also includes the positively selected genes observed in literature (Kosiol et al. 2008). We ensure that the q-values evaluated are not significant at 5%. A total of 26 genes are removed. Also, use only the alignments do not have any missing data such as gaps (```---```). This is to avoid any unwanted noise in this experiment. Use the following commands,

```
#The file called psg.txt contains the list of positively selected genes in real data analysis
cd /directory_with_all_the_gene_alignments/
find -type f -name 'psg.txt' -delete
```

For the primate dataset, the above filtering steps resulted in 6903 genes. These genes are considered as neutral genes for the primate branch. To construct a realistic simulation, randomly select $n = 500$ genes. You can use the following command to move files to another directory.

```
find . -mindepth 1 -maxdepth 1 -type f | shuf | head -n 10 | xargs -I{} mv {} dest_dir/
```

Note: 'dest_dir' is the directory with all the randomly selected genes to be used in further steps. Create a directory called 'modified' before running the above command. It is a best practice to keep a list of the filenames of the randomly generated $n$ genes to avoid any mishaps during this experiment. See the list of gene names and the corresponding number sequence (from 1 to 9566). 

The respective parameters estimated from $n$ genes (such as branch lengths and transition-transversion ratio) are used in ```EVOLVER``` for simulating positive selected codons. For simulating positive selection, one uses alternative hypothesis with either $\omega_2 = 4$ (moderate selective pressure) or $\omega_2 = 10$ (strong selective pressure). Then, these newly simulated positive codons are concatenated at the end of $n$ gene alignment. Thus, the modified alignments now has the positively selected codons and are assumed to be under the alternative hypothesis. The remaining 6,403 unmodified genes are assumed to be under the null hypothesis.

Now, analyse all the alignments with ```CODEML``` to evaluate the log-likelihoods under the null and alternative models of branch-site test. Then, use the ```simFDR.R```, to calculate the LRT, and the P-values and q-value significant at 5% confidence level. The steps are similar to the one in simulation analysis.
