
########################################################
##### Function to determine the power of the test ###### 
##### under the branch-site test and FDR          ######
######################################################## 


simFDR <- function(ln0,ln1A,ln1B,ln1C,lna0,lna1A,lna1B,lna1C) {

lnL0 <- c(lna0,ln0)

lnL1 <- pmax(c(lna0,ln0),c(lna1A,ln1A),c(lna1B,ln1B),c(lna1C,ln1C))
  
  #Likelihood ratio test statistic
  LRT <- 2 * ((lnL1)-(lnL0)) 
  
  #LRT == 0
  LRT_0 <- length(which(LRT == 0))
  
  #P-values combined
  pval <- pchisq(LRT, 1, lower.tail = FALSE) 
  
  #Dividing P-values by 2 under the mixture distribution
  p_halved <- ifelse(pval == 1, pval == pval, pval/2)  #Mixture distribution 
  
  #P-values significant within the 2.5% confidence interval
  p_sig_2.5 <- length(which(p_halved < 0.025))
  
  #P-values significant within the 5% confidence interval
  p_sig <- length(which(p_halved < 0.05))
  
  
  #Number of True positives are the simulations under alternative hypothesis
  z <- length(lna0)
  
  #True positive P-values - n (TP P < 0.05)
  true_positive_pval_signif <- sum(which(p_halved < 0.05) <= z)
  
  ### n (TP P = 1) ###
  true_positive_pval_1 <- sum(which(p_halved == 1) <= z)
  
  #CLASSICAL POWER
  Power <- (sum(which(p_halved < 0.05) <= z)/z)*100
  
  
  #######################################################
  ####False discovery rate (FDR)
  #CASE (1): BENJAMINI-HOCHBERG PROCEDURE
  
  #Loading qvalue 
  library('qvalue')
  
  #By default, the qvalue function uses BH-FDR method for the input P-values 
  BH.q.val <- qvalue::qvalue(p = p_halved)$qvalues
  
  BH.signif_q <- length(BH.q.val[BH.q.val < 0.05])  #number of q-value significance at 5%
  
  # Total positives are the total number of significant q-values
  BH.P <- length(which(BH.q.val < 0.05)) 
  
  #The true positives are the significant q-values that among the total positive hypothesis 
  # True positives are identified as
  BH.TP <- sum(which(BH.q.val < 0.05) <= z) 
  
  # Therefore, the false positives are
  BH.FP <- BH.P - BH.TP 
  
  #Proportion of false positives
  BH.PFP <- BH.FP/BH.P  #which is very close to 0.05!
  
  ### Percentage of FP ###
  BH.percent_FP <- BH.PFP*100
  
  #FDR POWER
  BH.FDRpower <- (BH.TP/z)*100 
  
  BH.output <- rbind(c(LRT_0, p_sig_2.5, p_sig, BH.signif_q, true_positive_pval_signif,
                 true_positive_pval_1, BH.TP, round(BH.percent_FP,digits=2), Power,  BH.FDRpower))
  
  colnames(BH.output) <- c("n(LRT=0)","n(P < 0.025)","n(P < 0.05)", "n(q < 0.05)",
                            "n(TP P < 0.05)","n(TP P = 1)","TP","percent_FP",
                            "classical_Power","FDR_power")
  rownames(BH.output) <- "BH-FDR"
  
  #(2) Storey's FDR
  #Checking if it is right or not
  row_pval <- paste("g",c(1:10000))
  dat_frame_pval <- data.frame(pval, row.names = row_pval)
  length(dat_frame_pval$pval[dat_frame_pval$pval < 0.05]) #Pvalue significant
  
  dat_frame_pval_lessthan1 <- dat_frame_pval[ which(dat_frame_pval$pval < 1) , ,drop = FALSE] #Length of P-values less than 1
  
  pval_new <- dat_frame_pval_lessthan1$pval
  
  ###By default, pi0est used lambda = seq(0.05,0.95,0.05), and 
  ####c("smoother","bootstrap") chooses the best method among the two methods
  pi0.estimate <- qvalue::pi0est(pval_new, pi0.method = "bootstrap")
  
  ###Proportion of true negatives over the P-values < 1 (under chi distribution)
  pi0_TN <- pi0.estimate$pi0     
  
  ####Total number of P-values: 10000
  ####The total proportion is 1
  pi1_P <- 1 - pi0_TN
  
  
  pi0est.total <- ((10000 - length(pval[pval < 1]))/10000) + (pi0_TN * length(pval_new)/10000)
  
  pi1est.total <- 1 - pi0est.total
  
  ###ST_q.val_case2 <- qvalue::qvalue(pval, fdr.level = 0.05, pi0 = pi0_est_total_case2)  
  ST.q.val <- qvalue::qvalue(pval, fdr.level = 0.05, pfdr = TRUE, pi0 = pi0est.total)
  
  ST.signif_q <- length(ST.q.val$qvalues[ST.q.val$qvalues < 0.05])
  
  # Total positives are the total number of significant q-values
  ST.P <- length(which(ST.q.val$qvalues < 0.05)) 
  
  #The true positives are the significant q-values that among the total positive hypothesis 
  # True positives are identified as
  ST.TP <- sum(which(ST.q.val$qvalues < 0.05) <= z) 
  
  # Therefore, there are 
  ST.FP <- ST.P - ST.TP #False positives
  
  #Proportion of false positives
  ST.percent.FP <- (ST.FP/ST.P)*100  #which is very close to 0.05!
  ST.FDRpower <- (ST.TP/z)*100 
  
  ##ST-FDR output
  ST.output <- rbind(c(length(pval_new), pi0est.total, pi1_P, pi0_TN, p_sig, 
                 ST.signif_q, true_positive_pval_signif, true_positive_pval_1, 
                 ST.TP, ST.percent.FP, Power, ST.FDRpower))
  
  colnames(ST.output) <- c("n(P < 1)", "pi0_dash", "p1", "pi0" , "n(P < 0.05)", 
                        "n(q < 0.05)", "n(TP P < 0.05)", "n(TP P = 1)", 
                        "TP","percent_FP","classical_Power","ST.FDR_power")
  rownames(ST.output) <- "ST-FDR"
  
  return(list(BH.output, ST.output))
  }

###Load the working directory containing the 
setwd("/Users/muthukumaranpanchaksaram/Documents/Proj20/paper/Results/Simresults/8sp/8a/8aw/8aw10/test")

############## Null hypothesis ################################
#Loading rst1 file from iteration 1
#Iteration 1 - Null model (H0) 
#Initial_omega=1 fix_omega = 1
ln0 <- read.table("rst1H0null")[,ncol(read.table("rst1H0null"))] #rst1 - lnL0 

#Loading rst1 file from iteration 2
#Alternate model: H1 
#estimate_omega=0 fix_omega = 1
ln1A <- read.table("rst1H0alt1")[,ncol(read.table("rst1H0alt1"))] #rst1 - lnL1 for 1st iteration

#Loading rst1 file from iteration 3
#Alternate model: H1
#estimate_omega=0 fix_omega = 2
ln1B <- read.table("rst1H0alt2")[,ncol(read.table("rst1H0alt2"))]

#Loading rst1 file from iteration 4
#Alternate model: H1
#estimate_omega=0 fix_omega = 4
ln1C <- read.table("rst1H0alt3")[,ncol(read.table("rst1H0alt3"))]

########Alternative hypothesis
#Loading rst1 file from iteration 1
lna0 <- read.table("rst1H1null")[,ncol(read.table("rst1H1null"))] #rst1 - lnL0 values

#Loading rst1 file from iteration 2
lna1A <- read.table("rst1H1alt1")[,ncol(read.table("rst1H1alt1"))] #rst1 - lnL1 for 1st iteration

#Loading rst1 file from iteration 3
lna1B <- read.table("rst1H1alt2")[,ncol(read.table("rst1H1alt2"))]

#Loading rst1 file from iteration 4
lna1C <- read.table("rst1H1alt3")[,ncol(read.table("rst1H1alt3"))]


#Applying simFDR() function
simFDR(ln0,ln1A,ln1B,ln1C,lna0,lna1A,lna1B,lna1C)





