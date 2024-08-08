
################################
####### FDR function ###########
################################ 


FDR <- function(ln0,ln1A,ln1B,ln1C,lna0,lna1A,lna1B,lna1C) {
  
  lnL0 <- ln0
  
  lnL1 <- pmax(ln0,ln1A,ln1B,ln1C)
  
  #LRT statistic for simulation under null hypothesis
  LRTnull <- 2 * ((lnL1)-(lnL0)) 
  
  ##Simulations under alternative hypothesis
  lnaL0 <- lna0
  
  lnaL1 <- pmax(lna0,lna1A,lna1B,lna1C)
  
  #LRT statistic for simulation under alternative hypothesis
  LRTalt <- 2 * ((lnaL1)-(lnaL0)) 
  
  
  #################################################################################################### 
  #Combining the first 1000 replicates simulated under alternative hypothesis with the subsequent 9000
  #simulated under null hypothesis####################################################################
  
  
  ###Set the number of alternative hypothesis
  #z represents the number of simulations under alternative hypothesis
  #z changes depending on the number of alternative hypothesis considered in the simulations
  z <- length(LRTalt)
  
  #Removing the first 1000 replicates from the vector LRTnull to get LRT_null
  LRT_null <- LRTnull[-(1:z)] #LRTnull is the LRT vector simulated for 10000 simulations under null
  
  #Adding the 1000 replicates from LRTalt combining to new LRT
  LRT <- c(LRTalt[1:z], LRT_null) #LRTalt is the LRT calculated for 1000 simulations under alternate hyp
  
  
  ###All the LRTs equal to 0
  LRT_0 <- length(LRT[LRT == 0])
  
  #P-values combined
  pval <- pchisq(LRT, 1, lower.tail = FALSE) 
  
  #Dividing P-values by 2 under the mixture distribution
  p_halved <- ifelse(pval == 1, pval == pval, pval/2)  #Mixture distribution 
  
  #P-values significant within the 2.5% confidence interval
  p_sig_2.5 <- length(which(p_halved < 0.025))
  
  #P-values significant within the 5% confidence interval
  p_sig <- length(which(p_halved < 0.05))
  
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
  
  BH.output <- c(LRT_0, p_sig_2.5, p_sig, BH.signif_q, true_positive_pval_signif,
              true_positive_pval_1, BH.TP, round(BH.percent_FP,digits=2), Power,  BH.FDRpower)
  
  
  #Checking if it is right or not
  row_pval <- paste("g",c(1:10000))
  dat_frame_pval <- data.frame(pval, row.names = row_pval)
  dat_frame_pval_not_1 <- dat_frame_pval[ which(dat_frame_pval$pval < 1) , ,drop = FALSE] #Length of P-values less than 1
  
  pval_new <- dat_frame_pval_not_1$pval
  
  ###By default, pi0est used lambda = seq(0.05,0.95,0.05), and 
  ####c("smoother","bootstrap") chooses the best method among the two methods
  pi0_estimate_case2 <- qvalue::pi0est(pval_new, pi0.method = "bootstrap")
  
  ###Proportion of true negatives over the P-values < 1 (under chi distribution)
  pi0_TN <- pi0_estimate_case2$pi0     
  
  ####Total number of P-values: 10000
  ####The total proportion is 1
  pi1_P <- 1 - pi0_TN
  
  pi0_est_total_case2 <- ((10000 - length(pval[pval < 1]))/10000) + (pi0_TN * length(pval_new)/10000)
  
  
  ###ST_q.val_case2 <- qvalue::qvalue(pval, fdr.level = 0.05, pi0 = pi0_est_total_case2)  
  ST_q.val_case2 <- qvalue::qvalue(pval, fdr.level = 0.05, pfdr = TRUE, pi0 = pi0_est_total_case2)
  
  dat_frame_qval <- data.frame(ST_q.val_case2$qvalues,row.names = row_pval)
  
  dat_frame_qval_rej <- dat_frame_qval[ which(dat_frame_qval$ST_q.val_case2.qvalues < 0.05), , drop = FALSE ] #Significance after applying FDR

  #FDR
  ST_P_case2 <- length(dat_frame_qval_rej$ST_q.val_case2.qvalues) # [1] Total Positives
  
  # However, only simulations 1 to 1,000 are true positives/considered under positive selection, thus:
  row_TP_case2 <- paste("g",c(1:z))
  ST_TP_frame_case2 <- dat_frame_qval_rej[rownames(dat_frame_qval_rej) %in% row_TP_case2, , drop = FALSE]
  ST_TP_case2 <- length(ST_TP_frame_case2$ST_q.val_case2.qvalues)
  
  
  # Therefore, there are 
  ST_FP_case2 <- ST_P_case2 - ST_TP_case2 #False positives
  
  #Proportion of false positives
  ST_PFP_case2 <- ST_FP_case2/ST_P_case2  #which is very close to 0.05!
  ST_FDRpower_case2<- (ST_TP_case2/z)*100 
  
  
  ####All 10,000 hypothesis under null
  length(which(pnull < 1))
  pi0_TN_null
  pi1_P_null
  pi0_est_total_case2_null
  length(ST_q.val_case2_null$pvalues[ST_q.val_case2_null$pvalues < 0.05])
  length(ST_q.val_case2_null$qvalues[ST_q.val_case2_null$qvalues < 0.05])
  ST_TP_case2_null
  ST_PFP_case2_null
  ST_FDRpower_case2_null
  
  #CASE (2): STOREY METHOD
  #Consider the P-values < 1
  pval_new <- pval[pval < 1] #Length of P-values less than 1
  
  ###By default, pi0est used lambda = seq(0.05,0.95,0.05), and 
  ####c("smoother","bootstrap") chooses the best method among the two methods
  #For example, when all the hypothesis is null, the pi0.method gives a pi0 
  # estimate of 90%.
  pi0_estimate <- qvalue::pi0est(pval_new, pi0.method="bootstrap")
  
  ###Proportion of true negatives over the P-values < 1 (under chi-square distribution)
  pi0_TN <- pi0_estimate$pi0     
  
  ####Total number of P-values: 10000
  ####The total proportion is 1
  pi1_P <- 1 - pi0_TN
  
  #The total pi0 estimate is calculated by the sum of pi0 estimated from two halves of the 
  #mixture distribution 
  pi0_est_total <- ((length(pval[pval == 1]))/10000) + (pi0_TN * length(pval_new)/10000)
  
  #Applying ST-FDR for P-values under mixture distribution at point mass 0 and chi-square
  #pi0 argument is specified to apply ST-FDR
  ST.q.val <- qvalue::qvalue(p_halved, fdr.level = 0.05, pfdr = TRUE, 
                             pi0 = pi0_est_total)$qvalues
  
  ###Significant q-values after ST-FDR ###
  ST.signif_q <- length(ST.q.val[ST.q.val < 0.05]) ### n (q < 0.05) ###
  
  # Total positives
  ST.P <- length(which(ST.q.val < 0.05)) 
  
  # True positives identified after FDR
  ST.TP <- sum(which(ST.q.val < 0.05) <= z) 
  
  ##False positives  
  ST.FP <- ST.P - ST.TP 
  
  #Proportion of false positives
  ST.PFP <- ST.FP/ST.P  #which is very close to 0.05!
  
  ### Percentage of FP ###
  ST.percent_FP <- ST.PFP*100
  
  #FDR POWER
  ST.FDRpower <- (ST.TP/z)*100 
  
  ##ST-FDR output
  ST.output <- c(LRT_0, p_sig_2.5, p_sig, ST.signif_q, true_positive_pval_signif,
                         true_positive_pval_1, ST.TP, round(ST.percent_FP,digits = 2), Power,  
                          ST.FDRpower)
  
  output <- rbind(BH.output,ST.output)
  colnames(output) <- c("n(LRT=0)","n(P < 0.025)","n(P < 0.05)", "n(q < 0.05)",
                        "n(TP P < 0.05)","n(TP P = 1)","TP","percent_FP",
                        "classical_Power","FDR_power")
  
  return(output)
  }

#Loading files######
setwd("/Users/muthukumaranpanchaksaram/Documents/Proj20/paper/Results/Simresults/8sp/8a/8aw/8aw8")

#Loading rst1 file from iteration 1
#Iteration 1 - Null model (H0) 
#Initial_omega=1 fix_omega = 1
ln0 <- read.table("rst1H0null")[,ncol(read.table("rst1H0null"))] #rst1 - lnL0 values

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

##Loading rst1 files output from CODEML program
lna0 <- read.table("rst1H1null")[,ncol(read.table("rst1H1null"))]  
lna1A <- read.table("rst1H1alt1")[,ncol(read.table("rst1H1alt1"))]  #lnL1 for 1st iteration
lna1B <- read.table("rst1H1alt2")[,ncol(read.table("rst1H1alt2"))] #for 2nd iteration
lna1C <- read.table("rst1H1alt3")[,ncol(read.table("rst1H1alt3"))] #lnL1 for 3rd iteration


#Applying function
FDR(ln0,ln1A,ln1B,ln1C,lna0,lna1A,lna1B,lna1C)




 
