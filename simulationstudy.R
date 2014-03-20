###########################################################################################################################
# SIMULATIONS: Why publishing everything is more effective than selective publishing of statistically significant results #
# 13-11-2013
###########################################################################################################################

##### DESCRIPTION SIMULATION STUDY #####
### GOAL: EXAMINE WHETHER PUBLISH EVERYTHING IS MORE EFFECTIVE THAN SELECTIVE PUBLISHING (REACTION ON DE WINTER & 
# HAPPEE (2013))
#
### DESCRIPTION SCENARIO 1:
# THE SAME SIMULATIONS AS IN THE PAPER BY DE WINTER AND HAPPEE WERE CONDUCTED, BUT WE USED THE SAMPLING 
# ERROR OF THE META-ANALYTIC ESTIMATE AFTER 40 PUBLICATIONS FOR BOTH CONDITIONS FOR EACH OF THE  5,000 REPLICATIONS INSTEAD 
# OF THE STANDARD DEVIATION OF THE META-ANALYTIC EFFECT AFTER 5,000 REPLICATIONS. IN BOTH CONDITIONS, PUBLISH EVERYTHING 
# (PE) AND SELECTIVE PUBLISHING (SP), A SINGLE REPLICATION WAS STOPPED WHEN 40 STUDIES GOT PUBLISHED. 
### SIMULATION CHARACTERISTICS:
# REPLICATIONS: 5,000
# EFFECT SIZE FOR GENERATING DATA (Etrue): 0.3
# SIGMA: 1
# ALPHA: .05
# SAMPLE SIZE PRIMARY STUDY (n): 50
### END DESCRIPTION SCENARIO 1
#
### DESCRIPTION SCENARIO 2:
# SCENARIO 2 IS SIMILAR TO SCENARIO 1 EXCEPT FOR THE EFFECT SIZE FOR GENERATING DATA (0 INSTEAD OF 0.3) AND THE STOPPING
# RULE. ACCORDING TO THE NEW STOPPING RULE EACH REPLICATION WAS STOPPED WHEN THE NULL HYPOTHESIS THAT THE POPULATION EFFECT
# SIZE IS AT LEAST SMALL (LARGER THAN OR EQUAL TO 0.2) IS REJECTED.
### SIMULATION CHARACTERISTICS:
# REPLICATIONS: 5,000
# EFFECT SIZE FOR GENERATING DATA (Etrue): 0
# SIGMA: 1
# ALPHA: .05
# SAMPLE SIZE PRIMARY STUDY (n): 50
### END DESCRIPTION SCENARIO 2
##### END DESCRIPTION SIMULATION STUDY #####

### Required packages
library(metafor)

###################################
##### SCENARIO 1: SIMULATIONS #####
###################################

### Settings for simulations
rep         <- 200 #5000
Etrue       <- 0.3
sigma       <- 1                
n           <- 50
alpha       <- 0.05
n.stu       <- 500 #1500
nrpubs      <- 20 #40
se          <- sigma/sqrt(n)
dsign       <- qnorm(.025,mean=0,sd=1,lower.tail=FALSE)
prec        <- 0.05

### Empty vectors/matrices to store results
d <- matrix(NA, nrow = rep, ncol = n.stu)

Emeta_SP <- matrix(NA, nrow = rep, ncol = n.stu)
Emeta_PE <- matrix(NA, nrow = rep, ncol = n.stu)

Epub_SP <- matrix(NA, nrow = rep, ncol = n.stu)
Epub_PE <- matrix(NA, nrow = rep, ncol = n.stu)

se_PE_RE <- matrix(NA, nrow = rep, ncol = n.stu)
se_SP_RE <- matrix(NA, nrow = rep, ncol = n.stu)

tau2_PE_RE <- matrix(NA, nrow = rep, ncol = n.stu)
tau2_SP_RE <- matrix(NA, nrow = rep, ncol = n.stu)

n.stu_PE <- vector(length = rep)
n.stu_SP <- vector(length = rep)

### Start simulations
for (r in 1:rep){
  print(r)
  ### Create new objects for storing scores and its ses
  scores <- matrix(NA, nrow = n, ncol = n.stu)
  se.scores <- numeric()
  
  for (i in 1:n.stu){
    
  # STEP 1: Generate data
    
    scores[,i] <- rnorm(n = n, mean = Etrue, sd = sigma)
    
    sd.scores <- sd(scores[ ,i]) 
    se.scores[i] <- sd.scores/sqrt(n)
    d[r,i] <- mean(scores[ ,i])
    
  # STEP 2: Establish published effect
    
    ### Publish everything (PE)
    Epub_PE[r,i] <- d[r,i]
    
    ### Publish only significant results, selective publishing (SP)
    if(i==1){
      Emeta_test <- 0
      ddiff <- Emeta_test-d[r,i]
    } else {
      ddiff <- Emeta_SP[r,i-1] - d[r,i]
    }
    
    ddiff_Z <- abs(ddiff/(1/sqrt(n)))
    
    if(ddiff_Z>dsign){
      Epub_SP[r,i] <- d[r,i]
    }
    
    ### Fill Emeta_SP with previous estimate if no study got published
    if(i==1){
        Emeta_SP[r,i] <- 0
      } else {
        Emeta_SP[r,i] <- Emeta_SP[r,i-1]
      }
    
  # STEP 3: Conduct random-effects meta-analysis (DerSimonian and Laird estimator for tau)
    
    ### Publish everything
    res_meta_PE <- rma(Epub_PE[r,1:i], sei = se.scores[1:i], method = "DL")
    se_PE_RE[r,i] <- res_meta_PE$se
    tau2_PE_RE[r,i] <- res_meta_PE$tau2
    Emeta_PE[r,i] <- res_meta_PE$b[1]
    
    ### Selective publishing
    if(!is.na(Epub_SP[r,i])) {
      res_meta_SP <- rma(Epub_SP[r, ][!is.na(Epub_SP[r, ])], sei = se.scores[!is.na(Epub_SP[r, ])], method = "DL")
      se_SP_RE[r,i] <- res_meta_SP$se 
      tau2_SP_RE[r,i] <- res_meta_SP$tau2
      Emeta_SP[r,i] <- res_meta_SP$b[1]
    }
    
  # STEP 4: Store the number of studies conducted
    
    ### Publish everything
    n.stu_PE[r] <- sum(!is.na(Epub_PE[r, ]))
    n.stu_SP[r] <- sum(!is.na(Epub_SP[r, ]))

    ### Break out loop if the number of published studies in selective publishing equals nrpubs
    if (n.stu_SP[r] == nrpubs)
      break
    
  }
}

################################
##### SCENARIO 1: ANALYSES #####
################################

### Remove NAs
tmp <- matrix(NA, nrow = rep, ncol = n.stu)

for(q in 1:n.stu) {
  for(m in 1:rep) {
    if(is.na(Epub_SP[m,q]) == FALSE) {
      tmp[m,q] <- Emeta_SP[m,q]
    } else {
      tmp[m,q] <- NA
    }
  }
}

test <- matrix(NA, nrow = rep, ncol = nrpubs)
for(z in 1:rep) {
	print(length(test[z,]))
	print(length(tmp[z, ][!is.na(tmp[z, ])]))
	if(length(test[z,])==length(tmp[z, ][!is.na(tmp[z, ])])){
		test[z, ] <- tmp[z, ][!is.na(tmp[z, ])]
	}
	else {
		test[z, ] <- c(tmp[z,!is.na(tmp[z, ])], rep(0,(length(test[z,])-length(tmp[z, ][!is.na(tmp[z, ])]))))
	}
}

### Mean of cumalative meta-analytic effect for PE
plot_Emeta_PE_mean <- colMeans(Emeta_PE[ ,1:nrpubs])

### Mean of cumalative meta-analytic effect for SP
plot_Emeta_SP_mean <- colMeans(test)

### Precision of the cumulative meta-effect (PE): Sd of cumulative meta-analytic effect for PE De Winter & Happee's method
plot_sd_PE_DWH <- apply(Emeta_PE[ ,1:nrpubs], 2, sd)

### Precision of the cumulative meta-effect (SP): Sd of cumulative meta-analytic effect for PE De Winter & Happee's method
plot_sd_SP_DWH <- apply(test, 2, sd)

### New way of calculating precision of the cumulative meta-effect (PE) (van Assen, van Aert, Nuijten, & Wicherts, 2013)
plot_sd_PE_new <- colMeans(se_PE_RE[ ,1:nrpubs])

### New way of calculating precision of the cumulative meta-effect (SP) (van Assen, van Aert, Nuijten, & Wicherts, 2013)
### Only select standard errors of published studies
se_SP_final <- matrix(numeric(),nrow=rep,ncol=nrpubs)

for(r in 1:rep){
	if(length(se_SP_final[r,])==length(se_SP_RE[r,!is.na(Epub_SP[r,])])){
		se_SP_final[r, ] <- se_SP_RE[r,!is.na(Epub_SP[r,])]
	}
	else {
		se_SP_final[r, ] <- c(se_SP_RE[r,!is.na(Epub_SP[r,])], 
			rep(0,length(se_SP_final[r, ])-length(se_SP_RE[r,!is.na(Epub_SP[r,])])))
	}
	  #se_SP_final[r, ] <- se_SP_RE[r,!is.na(Epub_SP[r,])]
}
plot_sd_SP_new <- colMeans(se_SP_final)

############################# PLOT 1 ########################################

### Plots similar to Figure 3 of De Winter and Happee
jpeg('rplot1.jpg',width=800,height=600)
plot(1:nrpubs,plot_Emeta_PE_mean,type="l",col="black",ylim=c(.15,.5),
     main="Cumulative Meta-Effect with Precision According to DWH")
#plot(1:nrpubs,plot_Emeta_PE_mean,type="l",col="black",ylim=c(.15,.5),
#     main="Cumulative Meta-Effect with Precision According to DWH")
print("plot1 done")

lines(1:nrpubs,plot_Emeta_SP_mean,col="red")

lines(1:nrpubs,plot_Emeta_PE_mean-plot_sd_PE_DWH,col="black",lty=2)
lines(1:nrpubs,plot_Emeta_PE_mean+plot_sd_PE_DWH,col="black",lty=2)

lines(1:nrpubs,plot_Emeta_SP_mean-plot_sd_SP_DWH,col="red",lty=2)
lines(1:nrpubs,plot_Emeta_SP_mean+plot_sd_SP_DWH,col="red",lty=2)
print("plot1 done")
dev.off()

############################# PLOT 2 ########################################

### Plots with new SE
jpeg('rplot2.jpg',width=800,height=600)

plot(1:nrpubs,plot_Emeta_PE_mean,type="l",col="black",ylim=c(.1,.5), lwd = 2, bty = "n", cex.lab = .9, cex.axis = .9,
     ylab = expression(paste(Mean ~ "\u00B1" ~ SE ~ meta-analysis)), xlab = "Publication number", yaxt = "n")
lines(1:nrpubs,plot_Emeta_SP_mean,col="red", lwd = 1.8)

axis(2, at = seq(.1, .5, .05), labels = seq(.1, .5, .05), las = 2, cex.axis = .9)
legend(x = 30, y = .5, legend = c("Publish everything", "Selective publishing"), lty = c(1,1), cex = .8, col = c("black", "red"))

lines(1:nrpubs,plot_Emeta_PE_mean-plot_sd_PE_new[1:nrpubs],col="black",lty=2)
lines(1:nrpubs,plot_Emeta_PE_mean+plot_sd_PE_new[1:nrpubs],col="black",lty=2)

lines(1:nrpubs,plot_Emeta_SP_mean-plot_sd_SP_new,col="red",lty=2)
lines(1:nrpubs,plot_Emeta_SP_mean+plot_sd_SP_new,col="red",lty=2)

### Standard error of the final cumulative meta-effect (PE)
se_PE_RE_final <- numeric()

for(r in 1:rep){
  se_PE_RE_final[r] <- se_PE_RE[r,nrpubs]
}

mean(se_PE_RE_final)

### tau^2 of the final cumulative meta-effect (PE)
tau2_PE_RE_final <- numeric()

for(r in 1:rep){
  tau2_PE_RE_final[r] <-tau2_PE_RE[r,nrpubs]
}

mean(tau2_PE_RE_final)

### How often tau^2 is 0 (in percentage) (PE)
(sum(tau2_PE_RE_final == 0))/rep

### Standard error of the final cumulative meta-effect (SP)
se_SP_RE_final <- numeric()

for(r in 1:rep){
  se_SP_RE_final[r] <- se_SP_RE[r,n.stu_PE[r]]
}

mean(se_SP_RE_final)

### tau^2 of the final cumulative meta-effect (SP)
tau2_SP_RE_final <- numeric()

for(r in 1:rep){
  tau2_SP_RE_final[r] <-tau2_SP_RE[r,n.stu_PE[r]]
}

mean(tau2_SP_RE_final)

### How often tau^2 is 0 (in percentage) (SP)
(sum(tau2_SP_RE_final == 0))/rep

###################################
##### SCENARIO 2: SIMULATIONS #####
###################################

### Clean workspace
rm(list=ls())   

### Settings for simulations
rep         <- 500 #5000
Etrue       <- 0
sigma       <- 1                
n           <- 50	
alpha       <- 0.1	#0.05
n.stu       <- 500	#1000
nrpubs      <- 40	#40
dsign       <- qnorm(.025,mean=0,sd=1,lower.tail=FALSE)

### Empty vectors/matrices to store results
d <- matrix(NA, nrow = rep, ncol = n.stu)

Emeta_SP <- matrix(NA, nrow = rep, ncol = n.stu)
Emeta_PE <- matrix(NA, nrow = rep, ncol = n.stu)

Epub_SP <- matrix(NA, nrow = rep, ncol = n.stu)
Epub_PE <- matrix(NA, nrow = rep, ncol = n.stu)

Epub_PE_pos <- matrix(NA, nrow = rep, ncol = n.stu)
Epub_PE_neg <- matrix(NA, nrow = rep, ncol = n.stu)

Epub_SP_pos <- matrix(NA, nrow = rep, ncol = n.stu)
Epub_SP_neg <- matrix(NA, nrow = rep, ncol = n.stu)

se_PE_RE <- matrix(NA, nrow = rep, ncol = n.stu)
se_SP_RE <- matrix(NA, nrow = rep, ncol = n.stu)

tau2_PE_RE <- matrix(NA, nrow = rep, ncol = n.stu)
tau2_SP_RE <- matrix(NA, nrow = rep, ncol = n.stu)

n.stu_PE <- vector(length = rep)
n.stu_SP <- vector(length = rep)

se_PE_RE_pos <- matrix(NA, nrow = rep, ncol = n.stu)
se_PE_RE_neg <- matrix(NA, nrow = rep, ncol = n.stu)

est_PE_RE_pos <- matrix(NA, nrow = rep, ncol = n.stu)
est_PE_RE_neg <- matrix(NA, nrow = rep, ncol = n.stu)

se_SP_RE_pos <- matrix(NA, nrow = rep, ncol = n.stu)
se_SP_RE_neg <- matrix(NA, nrow = rep, ncol = n.stu)

est_SP_RE_pos <- matrix(NA, nrow = rep, ncol = n.stu)
est_SP_RE_neg <- matrix(NA, nrow = rep, ncol = n.stu)

### Start simulations
for (r in 1:rep){
  
  ### Create new objects for storing scores and its ses
  scores <- matrix(NA, nrow = n, ncol = n.stu)
  se.scores <- numeric()
  
  for (i in 1:n.stu){
    
  # STEP 1: Generate data
    
    scores[ ,i] <- rnorm(n = n, mean = Etrue, sd = sigma)
    
    sd.scores <- sd(scores[ ,i]) 
    se.scores[i] <- sd.scores/sqrt(n)
    d[r,i] <- mean(scores[ ,i])
    
  # STEP 2: Establish published effect
    
    ### Publish everything (PE)
    Epub_PE[r,i] <- d[r,i]
    
    ### Publish only significant results, selective publishing (SP)
    if(i==1){
      Emeta_test <- 0
      ddiff <- Emeta_test-d[r,i]
    } else {
      ddiff <- Emeta_SP[r,i-1] - d[r,i]
    }
    
    ddiff_Z <- abs(ddiff/(1/sqrt(n)))
    
    if(ddiff_Z>dsign){
      Epub_SP[r,i] <- d[r,i]
    }
    
    ### Fill Emeta_SP with previous estimate if no study got published
    if(i==1){
      Emeta_SP[r,i] <- 0
    } else {
      Emeta_SP[r,i] <- Emeta_SP[r,i-1]
    }
    
  # STEP 3: Conduct random-effects meta-analysis (DerSimonian and Laird estim
    
    ### Publish everything
    res_meta_PE <- rma(Epub_PE[r,1:i], sei = se.scores[1:i], method = "DL")
    se_PE_RE[r,i] <- res_meta_PE$se
    tau2_PE_RE[r,i] <- res_meta_PE$tau2
    Emeta_PE[r,i] <- res_meta_PE$b[1]
    
    ### Split effects in first effect positive or negative (PE)
    if(Epub_PE[r,1] > 0) {
      Epub_PE_pos[r,i] <- Epub_PE[r,i]
    } else if(Epub_PE[r,1] < 0) {
      Epub_PE_neg[r,i] <- Epub_PE[r,i]
    }
    
    ### Random-effects meta-analysis for first effect is positive (PE)
    if(!is.na(Epub_PE_pos[r,i])) {
      res_meta_PE_pos <- rma(Epub_PE_pos[r, ][!is.na(Epub_PE_pos[r, ])], sei = se.scores, method = "DL")
      se_PE_RE_pos[r,i] <- res_meta_PE_pos$se 
      est_PE_RE_pos[r,i] <- res_meta_PE_pos$b[1]
    }
    
    ### Random-effects meta-analysis for first effect is negative (PE)
    if(!is.na(Epub_PE_neg[r,i])) {
      res_meta_PE_neg <- rma(Epub_PE_neg[r, ][!is.na(Epub_PE_neg[r, ])], sei = se.scores, method = "DL")
      se_PE_RE_neg[r,i] <- res_meta_PE_neg$se 
      est_PE_RE_neg[r,i] <- res_meta_PE_neg$b[1]
    }
    
    ### Selective publishing
    if(!is.na(Epub_SP[r,i])) {
      res_meta_SP <- rma(Epub_SP[r, ][!is.na(Epub_SP[r, ])], sei = se.scores[!is.na(Epub_SP[r, ])], method = "DL")
      se_SP_RE[r,i] <- res_meta_SP$se 
      tau2_SP_RE[r,i] <- res_meta_SP$tau2
      Emeta_SP[r,i] <- res_meta_SP$b[1]
    }
    
    ### Split effects in first effect positive or negative (SP)    
    if(!is.na(Epub_SP[r,i]) & Epub_SP[r,which(!is.na(Epub_SP[r, ]))][1] > 0) {
      Epub_SP_pos[r,i] <- Epub_SP[r,i]
    } else if(!is.na(Epub_SP[r,i]) & Epub_SP[r,which(!is.na(Epub_SP[r, ]))][1] < 0) {
      Epub_SP_neg[r,i] <- Epub_SP[r,i]
    }
    
    ### Random-effects meta-analysis for first effect is positive (SP)
    if(!is.na(Epub_SP_pos[r,i])) {
      res_meta_SP_pos <- rma(Epub_SP_pos[r, ][!is.na(Epub_SP_pos[r, ])], sei = se.scores[!is.na(Epub_SP_pos[r, ])], method = "DL")
      se_SP_RE_pos[r,i] <- res_meta_SP_pos$se 
      est_SP_RE_pos[r,i] <- res_meta_SP_pos$b[1]
    }
    
    ### Random-effects meta-analysis for first effect is negative (SP)
    if(!is.na(Epub_SP_neg[r,i])) {
      res_meta_SP_neg <- rma(Epub_SP_neg[r, ][!is.na(Epub_SP_neg[r, ])], sei = se.scores[!is.na(Epub_SP_neg[r, ])], method = "DL")
      se_SP_RE_neg[r,i] <- res_meta_SP_neg$se 
      est_SP_RE_neg[r,i] <- res_meta_SP_neg$b[1]
    }
    
  # STEP 4: Store the number of studies conducted
    
    ### Publish everything
    n.stu_PE[r] <- sum(!is.na(Epub_PE[r, ]))
    
    ### Selective publishing
    n.stu_SP[r] <- sum(!is.na(Epub_SP[r, ]))
     
    ### Break out loop if the null hypothesis that the population effect size is at least small (larger or equal to 0.2) is rejected
    if(is.na(Epub_SP[r,i]) == FALSE & pnorm(abs(Emeta_SP[r,i]), mean = 0.2, sd = se_SP_RE[r,i]) < alpha)
      break
    
  }
}

################################
##### SCENARIO 2: ANALYSES #####
################################

### Number of studies you need to CONDUCT & PUBLISH (PE)
n_stu_pub_PE <- numeric()

for(r in 1:rep){
  n_stu_pub_PE[r] <- which(pnorm(abs(Emeta_PE[r,]),mean=.2,sd=se_PE_RE[r,])<alpha)[1]
}

mean(n_stu_pub_PE)
sd(n_stu_pub_PE)

### Final cumulative meta-effect (PE)
Emeta_PE_final <- numeric()

for(r in 1:rep){
  Emeta_PE_final[r] <- Emeta_PE[r,n_stu_pub_PE[r]]
}

mean(Emeta_PE_final)

### Standard error of the final cumulative meta-effect (PE)
se_PE_RE_final <- numeric()

for(r in 1:rep){
  se_PE_RE_final[r] <- se_PE_RE[r,n_stu_pub_PE[r]]
}

mean(se_PE_RE_final)

### tau^2 of the final cumulative meta-effect (PE)
tau2_PE_RE_final <- numeric()

for(r in 1:rep){
  tau2_PE_RE_final[r] <-tau2_PE_RE[r,n_stu_pub_PE[r]]
}

mean(tau2_PE_RE_final)
sd(tau2_PE_RE_final)

### How often tau^2 is 0 (in percentage) (PE)
(sum(tau2_PE_RE_final == 0))/rep

### Final cumulative meta-effect for positive first study (PE)
est_PE_RE_pos_final <- numeric()

for(r in 1:rep){
  est_PE_RE_pos_final[r] <- est_PE_RE_pos[r,n_stu_pub_PE[r]]
}

mean(est_PE_RE_pos_final, na.rm = TRUE)

### Final cumulative meta-effect for negative first study (PE)
est_PE_RE_neg_final <- numeric()

for(r in 1:rep){
  est_PE_RE_neg_final[r] <- est_PE_RE_neg[r,n_stu_pub_PE[r]]
}

mean(est_PE_RE_neg_final, na.rm = TRUE)

### Standard error of the final cumulative meta-effect for positive first study (PE)
se_PE_RE_pos_final <- numeric()

for(r in 1:rep){
  se_PE_RE_pos_final[r] <- se_PE_RE_pos[r,n_stu_pub_PE[r]]
}

mean(se_PE_RE_pos_final, na.rm = TRUE)

### Standard error of the final cumulative meta-effect for negative first study (PE)
se_PE_RE_neg_final <- numeric()

for(r in 1:rep){
  se_PE_RE_neg_final[r] <- se_PE_RE_neg[r,n_stu_pub_PE[r]]
}

mean(se_PE_RE_neg_final, na.rm = TRUE)

### Selective publishing: number of studies you need to PUBLISH (SP)
n_stu_pub_SP <- apply(se_SP_RE,1,function(x) length(x[!is.na(x)]))
n_stu_cond_SP <- apply(se_PE_RE,1,function(x) length(x[!is.na(x)]))

mean(n_stu_cond_SP)
sd(n_stu_cond_SP)

mean(n_stu_pub_SP)
sd(n_stu_pub_SP)

### Final cumulative meta-effect (SP)
Emeta_SP_final <- numeric()

for(r in 1:rep){
  Emeta_SP_final[r] <- Emeta_SP[r,n_stu_cond_SP[r]]
}

mean(Emeta_SP_final)

### Standard error of the final cumulative meta-effect (SP)
se_SP_RE_final <- numeric()

for(r in 1:rep){
  se_SP_RE_final[r] <- se_SP_RE[r,n_stu_cond_SP[r]]
}

mean(se_SP_RE_final)

### tau^2 of the final cumulative meta-effect (SP)
tau2_SP_RE_final <- numeric()

for(r in 1:rep){
  tau2_SP_RE_final[r] <-tau2_SP_RE[r,n_stu_cond_SP[r]]
}

mean(tau2_SP_RE_final)
sd(tau2_SP_RE_final)

### How often tau^2 is 0 (in percentage) (SP)
(sum(tau2_SP_RE_final == 0))/rep

### Final cumulative meta-effect for positive first study (SP)
est_SP_RE_pos_final <- numeric()

for(r in 1:rep){
  est_SP_RE_pos_final[r] <- est_SP_RE_pos[r,n_stu_cond_SP[r]]
}

mean(est_SP_RE_pos_final, na.rm = TRUE)

# Final cumulative meta-effect for negative first study (SP)
est_SP_RE_neg_final <- numeric()

for(r in 1:rep){
  est_SP_RE_neg_final[r] <- est_SP_RE_neg[r,n_stu_cond_SP[r]]
}

mean(est_SP_RE_neg_final, na.rm = TRUE)

### Standard error of the final cumulative meta-effect for positive first study  (SP)
se_SP_RE_pos_final <- numeric()

for(r in 1:rep){
  se_SP_RE_pos_final[r] <- se_SP_RE_pos[r,n_stu_cond_SP[r]]
}

mean(se_SP_RE_pos_final, na.rm = TRUE)

### Standard error of the final cumulative meta-effect for negative first study (SP)
se_SP_RE_neg_final <- numeric()

for(r in 1:rep){
  se_SP_RE_neg_final[r] <- se_SP_RE_neg[r,n_stu_cond_SP[r]]
}

mean(se_SP_RE_neg_final, na.rm = TRUE)
dev.off()
quartz()