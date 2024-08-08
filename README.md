# Statistical properties of branch-site test and FDR

This repository contains the files to recreate the analysis in Chapter 2 - Statistical properties of the branch-site test and false discovery rate to detect episodic positive selection.

To run the analysis, users are required to download and compile ```EVOLVER``` and ```CODEML``` programs available in ```PAML``` package, see [http://abacus.gene.ucl.ac.uk/software/paml.html](http://abacus.gene.ucl.ac.uk/software/paml.html). Use bash/LINUX commands to run PAML programs.

## Simulating codon alignments under different levels of positive selection
To simulate alignments codon alignments, one can use EVOLVER program to PAML package ```EVOLVER```. 


## Branch-site test for positive selection

### Log-likelihood estimation
For the given sequence alignment, the log-likelihoods are evaluated under the null model (l0) and alternative model (l1). Under the null model, the log-likelihood (l0) is estimated by fixing the omega to 1, whereas the omega and other free parameters are estimated under the alternative model. 

The log-likelihood under alternative model is estimated for three iterations, each initialised with omega = 1, 2, and 4 (see initial_omega in codeml.ctl). And then take the largest log-likelihood among the all log-likelihood estimates under both the models (null and alternative). This is to ensure that the log-likelihood estimated under the alternative model is less than equal to log-likelihood estimated under the alternative model, (i.e., l1 <= l0)



### Likelihood ratio test, P-values and FDR


## Real dataset

## Simulating model misspecification






