### smoking data
data<-read.table("nmesdata.txt", header=TRUE,sep=',')
data<-subset(data, packyears>0)  ## smoker
## complete case analysis
data<-subset(data, select=c(packyears, AGESMOKE, LASTAGE, MALE, RACE3, beltuse, 
				educate, marital,SREGION, POVSTALB,	HSQACCWT, TOTALEXP))
data<-na.omit(data)
#additional variables (quadratic transformations)
data$LASTAGE2<-data$LASTAGE^2
data$AGESMOKE2<-data$AGESMOKE^2
## factor variables
data$MALE<-factor(data$MALE, ordered=F)
data$RACE3<-factor(data$RACE3, ordered=F)
data$marital<-factor(data$marital, ordered=F)
data$SREGION<-factor(data$SREGION, ordered=F)
data$educate<-factor(data$educate, ordered=F)
data$beltuse<-factor(data$beltuse, ordered=F)
data$POVSTALB<-factor(data$POVSTALB, ordered=F)
data$logpack = log( data$packyears )

#consider the subset with just data$TOTALEXP > 0
data = subset(data, TOTALEXP > 0)
#Also remove packyears, which was transformed to logpack
data = subset(data, select = -c(packyears))

#relabel the treatment and outcome as a and y
names(data)[which(names(data) == "TOTALEXP")] = "y"
names(data)[which(names(data) == "logpack")] = "a"
#do a log transformation on y
data$y = log(data$y)

gaussianKernel = function(A, a, h){
	#the difference is
	u = A - a
	#then, the kernel is
	kernel = (1/sqrt(2*pi*h^2))*exp(-0.5*(u/h)^2)
	return(kernel)
}

#This function is useful for running a GAM
#for the PS model and outcome model.
#This simply increases the number
#of continuous variables from 4 (the default) to 5.
SL.gam.new <- function(..., cts.num = 5) {
	SL.gam(..., cts.num = 5)
}

#This function estimates the causal effect
#\psi_h(a), defined as the "smoothed causal effect"
#(without trimming).
#This is based on the EIF given in the paper.
#Requires a dataset with X, A, Y,
#a bandwidth h, and a set of treatments a.eval.
#Thus, this function estimates \psi_h(a.eval).
estPsih.split = function(data, h, a.eval,
               a.seq = seq(-2, 3, by = 0.05), nsplits = 2,
               sl.lib.ps = c("SL.earth", "SL.gam", 
                             "SL.glm", "SL.glm.interaction",
                             "SL.mean", "SL.ranger", "SL.rpart"),
               sl.lib.outcome = c("SL.earth", "SL.gam", 
                                  "SL.glm", "SL.glm.interaction",
                                  "SL.mean", "SL.ranger", "SL.rpart"),
               weights = NULL){
  
  #number rows
  n = nrow(data)
  
  #first, divide the data into
  #a train-test split. Define the split indexes:
  splitIndex = sample(rep(1:nsplits, ceiling(n/nsplits))[1:n])
  
  #we will ultimately compute the estimated \psi_h
  #for each a in a.eval.
  #This will be a matrix that collects
  #the estimated EIF values for each 
  #observation and treatment combination.
  #rows are observations; columns are treatment values.
  psih.ifvalues = matrix(nrow = n, ncol = length(a.eval))
  
  #now we will repeat this process for
  #each sample split...
  for(split in 1:nsplits){
    #divide into training and testing
    #which splitIndex is equal to this split?
    train <- splitIndex != split
    test <- splitIndex == split
    if (nsplits == 1) {
      train <- test
    }
    #then, the datasets are
    data.train = data[train,]
    data.test = data[test,]
    
    #it'll also be useful to have a covariate matrix
    #for the training data and for the test data.
    X.train = subset(data.train, select = -c(a,y))
    X.test = subset(data.test, select = -c(a,y))
    
    #first, estimate pi
    psMod.sl = SuperLearner(Y = data.train$a,
                            X = X.train, newX = X.test,
                            SL.library = sl.lib.ps,
                            obsWeights = weights[train])
    #thus, the mean estimate is (on the test data):
    meanA.est = psMod.sl$SL.predict
    #To flexibly estimate the variance, it'll be helpful
    #to first compute squared residuals from this model:
    meanA.est.train = as.numeric(predict(psMod.sl, newdata = X.train)$pred)
    piRes2 = (data.train$a - meanA.est.train)^2
    #Then, a flexible model for the variance is:
    pi2mod = SuperLearner(Y = piRes2,
                          X = X.train, newX = X.test,
                          SL.library = sl.lib.ps,
                          obsWeights = weights[train])
    varA.est = pi2mod$SL.predict
    #Then, \hat{\pi}(X, A) is:
    piHat.XA = dnorm(data.test$a, mean = meanA.est, sd = sqrt(varA.est))
    
    #meanwhile, the estimated mu is...
    muModel = SuperLearner(Y = data.train$y,
                           X = subset(data.train, select = -c(y)),
                           newX = subset(data.test, select = -c(y)),
                           SL.library = sl.lib.outcome,
                           obsWeights = weights[train])
    
    #First, compute \hat{\mu}(X, A): 
    muHat.XA = as.numeric(predict(muModel)$pred)
    
    #computing the IPW residuals
    ipwRes = (data.test$y - muHat.XA)/piHat.XA
    
    #It'll be useful to put this into an n * length(a.eval) matrix,
    #where each column is just ipwRes
    # ipwRes.mat = matrix(rep(ipwRes, length(a.eval)), 
    #                     ncol = length(a.eval))
    
    #now produce a matrix of estimates from the mu model,
    #where the rows correspond to subjects and the columns
    #correspond to the a.seq.
    #Thus, we must create a new data.frame, where
    #each row corresponds to a subject/a.seq combination.
    #first, create the covariate data
    mu.preddata = X.test[rep(seq_len(nrow(X.test)), each = length(a.seq)),]
    #now add the treatment data, corresponding to each a.seq.
    mu.preddata$a = rep(a.seq, times = nrow(data.test))
    
    #now predict across these a.seq
    #this will create a matrix, where
    #each row corresponds to a subject
    #and each column corresponds to a a.seq
    muHat.a0 = matrix(as.numeric(predict(muModel, newdata = mu.preddata)$pred),
                      ncol = length(a.seq), byrow = TRUE)
    #Note that this is \hat{mu}(X, a_0) in the paper,
    #which is used within the integral of the estimator.
    
    #when computing the integral manually,
    #we need to know the "gap" between the a0s
    #(in order to compute a Riemann sum)
    a.seq.len = a.seq[2] - a.seq[1]
    
    #Now we just need to compute the kernel terms,
    #which we'll store in matrices.
    #kernel.A will be n * length(a.eval), where each
    #column contains K_h(A) for a particular a.
    kernel.A = sapply(a.eval, gaussianKernel, A = data.test$a, h = h)
    #kernel.a0 will be length(a.seq) * length(a.eval), where each
    #column contains K_h(a0) for a particular a.
    kernel.a0 = sapply(a.eval, gaussianKernel, A = a.seq, h = h)
    #normalize the kernels
    kernel.A = t(apply(kernel.A, MARGIN = 1, FUN = function(x) x/colSums(kernel.a0*a.seq.len)))
    kernel.a0 = t(apply(kernel.a0, MARGIN = 1, FUN = function(x) x/colSums(kernel.a0*a.seq.len)))

      #the estimator consists of two terms.
      #The first term is:
      psih.1 = kernel.A*ipwRes
      
      #The second term then is:
      psih.2 = a.seq.len*muHat.a0%*%kernel.a0
      
      #then, the IF values are
      psih.ifvalues[test,] = psih.1 + psih.2
  }
  
  #now we'll construct point estimates and confidence intervals
  #We will do this pointwise, for each treatment value.
  
  #the point estimates are
  psih.est = colMeans(psih.ifvalues)
  
  #then, the estimated standard error is
  psih.se = apply(psih.ifvalues, MARGIN = 2, FUN = sd)/sqrt(n)
  #then, the confidence intervals are
  psih.ci.l = psih.est - qnorm(0.975)*psih.se
  psih.ci.u = psih.est + qnorm(0.975)*psih.se
  
  #return a list containing point estimates and CIs
  results = list(
    psih.est = psih.est,
    psih.ci.l = psih.ci.l,
    psih.ci.u = psih.ci.u)
  
  return(results)
}

#load required libraries
library(SuperLearner)

# ESTIMATE THE DOSE-RESPONSE CURVE WITHOUT TRIMMING
a.eval = seq(
	quantile(data$a, prob = 0.05),
	quantile(data$a, prob = 0.95),
	length = 10)

set.seed(123)
#fit for small bandwidth (h = 0.21)
psih.sl.h21 = estPsih.split(data = subset(data, select = -c(HSQACCWT)),
	h = 0.21,
	a.eval = a.eval,
	a.seq = seq(-4, 8, by = 0.1),
	sl.lib.ps = c("SL.earth", "SL.gam.new", "SL.glm", "SL.glm.interaction", 
        "SL.mean", "SL.ranger", "SL.rpart"),
	sl.lib.outcome = c("SL.earth", "SL.gam.new", "SL.glm", "SL.glm.interaction", 
        "SL.mean", "SL.ranger", "SL.rpart"),
	weights = data$HSQACCWT/sum(data$HSQACCWT),
	nsplits = 2)
set.seed(123)
#fit for large bandwidth (h = 0.92)
psih.sl.h92 = estPsih.split(data = subset(data, select = -c(HSQACCWT)),
	h = 0.92,
	a.eval = a.eval,
	a.seq = seq(-4, 8, by = 0.1),
	sl.lib.ps = c("SL.earth", "SL.gam.new", "SL.glm", "SL.glm.interaction", 
        "SL.mean", "SL.ranger", "SL.rpart"),
	sl.lib.outcome = c("SL.earth", "SL.gam.new", "SL.glm", "SL.glm.interaction", 
        "SL.mean", "SL.ranger", "SL.rpart"),
	weights = data$HSQACCWT/sum(data$HSQACCWT),
	nsplits = 2)

#To implement the trimming estimator,
#it'll be helpful to have functions for S_t and S'_t
getSt = function(piHat, epsilon = 10^(-4), t = 0.1){
	St = pnorm(piHat - t, mean = 0, sd = epsilon)
	return(St)
}
getSt.deriv = function(piHat, epsilon = 10^(-4), t = 0.1){
	St.deriv = dnorm(piHat - t, mean = 0, sd = epsilon)
	return(St.deriv)
}
# THIS IS THE PRIMARY FUNCTION TO COMPUTE
# OUR TWO "TRIMMED DOSE RESPONSE" ESTIMATORS
# WHEN THE TRIMMING THRESHOLD IS FIXED.
estPsih.trimmed.split = function(data,
	a.eval, t, h, epsilon = 10^(-2),
	a.seq = seq(-2, 3, by = 0.05),
	nsplits = 2,
	sl.lib.ps = c("SL.earth", "SL.gam", 
	    "SL.glm", "SL.glm.interaction",
	    "SL.mean", "SL.ranger", "SL.rpart"),
	sl.lib.outcome = c("SL.earth", "SL.gam", 
	    "SL.glm", "SL.glm.interaction",
	    "SL.mean", "SL.ranger", "SL.rpart"),
	weights = NULL){

	#sample size
	n = nrow(data)

	#first, divide the data into
	#a train-test split. Define the split indexes:
	splitIndex = sample(rep(1:nsplits, ceiling(n/nsplits))[1:n])
	
	#we will ultimately compute the estimated \psi_h
	#for each a in a.eval.
	#We'll consider two estimators:
	# 1) ratio(EIFs)
	# 2) EIF(ratio)
	#This will be a matrix that collects
	#the estimated EIF values for each 
	#observation and treatment combination.
	#rows are observations; columns are treatment values.
	#Note that, for the first estimator, we must store
	#the IF values for the numerator and denominator terms.
	psih.ifvalues1.num = matrix(nrow = n, ncol = length(a.eval))
	psih.ifvalues1.denom = matrix(nrow = n, ncol = length(a.eval))
	psih.ifvalues2 = matrix(nrow = n, ncol = length(a.eval))

	#Finally, we'll also consider a more simple EIF-based estimator,
	#where we compute the EIF for the dose-response curve,
	#but only within the trimmed sample.
	psih.naive.ifvalues = matrix(nrow = n, ncol = length(a.eval))
	#It'll be helpful to also have a matrix of trimming indicators.
	ind.plugin = matrix(nrow = n, ncol = length(a.eval))

	#now we will repeat this process for
	#each sample split...
	for(split in 1:nsplits){
		#divide into training and testing
		#which splitIndex is equal to this split?
		train <- splitIndex != split
		test <- splitIndex == split
		 if (nsplits == 1) {
                train <- test
            }
		#then, the datasets are
		data.train = data[train,]
		data.test = data[test,]

		#it'll also be useful to have a covariate matrix
		#for the training data and for the test data.
		X.train = subset(data.train, select = -c(a,y))
		X.test = subset(data.test, select = -c(a,y))

		#first, estimate pi
		psMod.sl = SuperLearner(Y = data.train$a,
		X = X.train, newX = X.test,
		SL.library = sl.lib.ps,
		obsWeights = weights[train])
		#thus, the mean estimate is (on the test data):
		meanA.est = psMod.sl$SL.predict
		#To flexibly estimate the variance, it'll be helpful
		#to first compute squared residuals from this model:
		meanA.est.train = as.numeric(predict(psMod.sl, newdata = X.train)$pred)
		piRes2 = (data.train$a - meanA.est.train)^2
		#Then, a flexible model for the variance is:
		pi2mod = SuperLearner(Y = piRes2,
			X = X.train, newX = X.test,
        	SL.library = sl.lib.ps,
        	obsWeights = weights[train])
    	varA.est = pi2mod$SL.predict
		#Then, \hat{\pi}(X, A) is:
		piHat.XA = dnorm(data.test$a, mean = meanA.est, sd = sqrt(varA.est))

		#meanwhile, the estimated mu(x,a) is...
		muModel = SuperLearner(Y = data.train$y,
			X = subset(data.train, select = -c(y)),
			newX = subset(data.test, select = -c(y)),
			SL.library = sl.lib.outcome,
			obsWeights = weights[train])
		
		#First, compute \hat{\mu}(X, A): 
		muHat.XA = as.numeric(muModel$SL.predict)

		#computing the IPW residuals
		#(again, using the test data)
		ipwRes = (data.test$y - muHat.XA)/piHat.XA

		#Now compute the smoothed indicator and its derivative
		St.A = getSt(piHat = piHat.XA, epsilon = epsilon, t = t)
		St.deriv.A = getSt.deriv(piHat = piHat.XA, epsilon = epsilon, t = t)

		#now produce a matrix of estimates from the mu model,
		#where the rows correspond to subjects and the columns
		#correspond to the a.seq.
		#Thus, we must create a new data.frame, where
		#each row corresponds to a subject/a.seq combination.
		#Note that we only use the test data here.
		#first, create the covariate data
		mu.preddata.a0 = X.test[rep(seq_len(nrow(X.test)), each = length(a.seq)),]
		#now add the treatment data, corresponding to each a.seq.
		mu.preddata.a0$a = rep(a.seq, times = nrow(data.test))

		#now predict across these a.seq
		#this will create a matrix, where
		#each row corresponds to a subject
		#and each column corresponds to a a.seq
		muHat.a0 = matrix(as.numeric(predict(muModel, newdata = mu.preddata.a0)$pred),
			ncol = length(a.seq), byrow = TRUE)
		#Note that this is \hat{mu}(X, a_0) in the paper,
		#which is used within the integral of the estimator.

		#now compute \hat{\pi}(a0|X) for each a0
		#this will again correspond to a matrix,
		#where each column is a particular a0
		piHat.a0 = sapply(a.seq, FUN = dnorm, mean = meanA.est, sd = sqrt(varA.est))
		#compute the corresponding St and St.deriv, which are matrices
		St.a0 = getSt(piHat.a0, epsilon = epsilon, t = t)
		St.deriv.a0 = getSt.deriv(piHat.a0, epsilon = epsilon, t = t)

		#Finally, we'll compute the simple trimming indicator
		piHat.a = sapply(a.eval, FUN = dnorm, mean = meanA.est, sd = sqrt(varA.est))
		ind.plugin[test,] = ifelse(piHat.a > t, 1, 0)

		#when computing the integral manually,
		#we need to know the "gap" between the a0s
		#(in order to compute a Riemann sum)
		a.seq.len = a.seq[2] - a.seq[1]

		for(a in a.eval){
			#Now compute the kernel term
			#(using only the test data)
			kernel.A = gaussianKernel(A = data.test$a, a = a, h = h)
			kernel.a0 = gaussianKernel(A = a.seq, a = a, h = h)
			#normalize the kernels
			kernel.A = kernel.A/sum(kernel.a0*a.seq.len)
			kernel.a0 = kernel.a0/sum(kernel.a0*a.seq.len)

			#First, the estimator will depend
			#on plug-in estimators:
			psihat.pi.num = rowSums(
				t(apply(muHat.a0*St.a0*a.seq.len,
				MARGIN = 1, FUN = function(x) x*kernel.a0)))
			psihat.pi.denom = rowSums(
				t(apply(St.a0*a.seq.len,
				MARGIN = 1, FUN = function(x) x*kernel.a0)))

			#Now note that, beyond the plug-in estimator parts,
			#\hat{\psi}_num as three terms,
			#and \hat{\psi}_denom as two terms.

			#numerator terms one and two
			psihat.num.term1 = kernel.A*ipwRes*St.A
			psihat.num.term2 = kernel.A*muHat.XA*St.deriv.A
			#denominator term one
			psihat.denom.term1 = kernel.A*St.deriv.A

			#then, the remaining terms for the numerator and denominator are:
			psihat.num.term3 = -rowSums(t(apply(muHat.a0*piHat.a0*St.deriv.a0*a.seq.len,
				MARGIN = 1, FUN = function(x) x*kernel.a0)))
			psihat.denom.term2 = -rowSums(t(apply(piHat.a0*St.deriv.a0*a.seq.len,
				MARGIN = 1, FUN = function(x) x*kernel.a0)))

			#Now compute the psihat.num and psihat.denom
			psihat.num = psihat.pi.num + 
				psihat.num.term1 + psihat.num.term2 + psihat.num.term3
			psihat.denom = psihat.pi.denom +
				psihat.denom.term1 + psihat.denom.term2

			#then, the particular estimate is
			#(specifically for the test set)
			psih.ifvalues1.num[test,which(a.eval == a)] = psihat.num
			psih.ifvalues1.denom[test,which(a.eval == a)] = psihat.denom

			#We'll also consider the EIF(ratio) estimator,
			#which requires us to compute the plug-in estimator.
			psihat.pi = mean(psihat.pi.num)/mean(psihat.pi.denom)

			#The alternative estimator would be
			#(specifically for the test set)
			psih.ifvalues2[test,which(a.eval == a)] = psihat.pi*(1 - psihat.denom/mean(psihat.pi.denom)) +
				psihat.num/mean(psihat.pi.denom)

			#Finally, we'll compute the "naive" estimator that just computes
			#the EIF for the dose-response curve within the trimmed set.
			#The estimator consists of two terms.
			#The first term is:
			psih.naive.1 = kernel.A*ipwRes

			#The second term then is:
			psih.naive.2 = rowSums(
				t(apply(muHat.a0*a.seq.len,
				MARGIN = 1, FUN = function(x) x*kernel.a0)))

			#then, the particular estimate is
			psih.naive.ifvalues[test,which(a.eval == a)] = psih.naive.1 + psih.naive.2
		}
	}
	#now we'll construct point estimates and confidence intervals
	#We will do this pointwise, for each treatment value.

	#the point estimates are
	psih.est1.num = colMeans(psih.ifvalues1.num)
	psih.est1.denom = colMeans(psih.ifvalues1.denom)
	psih.est1 = psih.est1.num/psih.est1.denom
	psih.est2 = colMeans(psih.ifvalues2)

	#Finally, the naive estimator will take the influence function
	#values, but only average across those with \hat{\pi} > t:
	psih.naive.est = colMeans(psih.naive.ifvalues*ind.plugin)/colMeans(ind.plugin)

	#the standard error for the first estimator is
	psih.se1 = vector(length = length(a.eval))
	for(j in 1:length(a.eval)){
		psih.se1[j] = sqrt( var((psih.ifvalues1.num[,j] - 
								psih.est1[j]*psih.ifvalues1.denom[,j])/psih.est1.denom[j])/n )
	}
	#and the standard error for the second estimator it is
    psih.se2 = sqrt(apply(psih.ifvalues2, MARGIN = 2, FUN = var)/n)

    #The standard error for the naive estimator
    #takes the empirical variance among the estimated IF values
    #where \hat{\pi} > t. Note that we also need to
    #compute the sample size in the non-trimmed set.
    psih.naive.se = vector(length = length(a.eval))
    for(j in 1:length(a.eval)){
    	ind.curr = ind.plugin[,j]
    	psih.naive.se[j] = sqrt( var(psih.naive.ifvalues[ind.curr == 1,j])/sum(ind.curr == 1) )
    }

    #then, the confidence intervals are
    psih.ci1.l = psih.est1 - qnorm(0.975)*psih.se1
    psih.ci1.u = psih.est1 + qnorm(0.975)*psih.se1
    psih.ci2.l = psih.est2 - qnorm(0.975)*psih.se2
    psih.ci2.u = psih.est2 + qnorm(0.975)*psih.se2
    psih.naive.ci.l = psih.naive.est - qnorm(0.975)*psih.naive.se
    psih.naive.ci.u = psih.naive.est + qnorm(0.975)*psih.naive.se

    #return a list containing point estimates and CIs
    results = list(
    	psih.est1.num = psih.est1.num,
    	psih.est1.denom = psih.est1.denom,
    	psih.est1 = psih.est1,
    	psih.est2 = psih.est2,
    	psih.naive.est = psih.naive.est,
    	psih.ci1.l = psih.ci1.l,
    	psih.ci1.u = psih.ci1.u,
    	psih.ci2.l = psih.ci2.l,
    	psih.ci2.u = psih.ci2.u,
    	psih.naive.ci.l = psih.naive.ci.l,
    	psih.naive.ci.u = psih.naive.ci.u)

    return(results)
}

#It can be useful to first understand what
#the estimated propensity scores, \pi(A|X), look like.

#estimating with SuperLearner
psMod.sl = SuperLearner(Y = data$a,
		X = subset(data, select = -c(a,y,HSQACCWT)),
		SL.library = c("SL.earth", "SL.gam.new", "SL.glm", "SL.glm.interaction", 
        "SL.mean", "SL.ranger", "SL.rpart"),
		obsWeights = data$HSQACCWT/sum(data$HSQACCWT))

#thus, the mean and SD estimates are...
#(note that we need to use the test data here)
meanA.est.sl = as.numeric(psMod.sl$SL.predict)
piRes2 = (data$a - meanA.est.sl)^2
pi2mod <- SuperLearner(Y = piRes2,
			X = subset(data, select = -c(a,y,HSQACCWT)),
        	SL.library = c("SL.earth", "SL.gam.new", "SL.glm", "SL.glm.interaction", 
        "SL.mean", "SL.ranger", "SL.rpart"),
        	obsWeights = data$HSQACCWT/sum(data$HSQACCWT))
pi2mod.vals <- pi2mod$SL.predict
#Then, \hat{\pi}(X, A) is:
piHat.XA.sl = dnorm(data$a, mean = meanA.est.sl, sd = sqrt(pi2mod.vals))
		
# VISUALIZATION OF ESTIMATED PROPENSITY SCORES
plot(data$a, piHat.XA.sl, pch = 16, cex = 0.25,
	xlab = "", ylab = "")
title(xlab = "Treatment A", line = 2.25, cex.lab = 1.25)
title(ylab = expression(paste("Estimated Propensity Score ", hat(pi), "(A|X)")), line = 2.25, cex.lab = 1.25)
#where Zhao et al. estimated treatment effects
abline(v = quantile(data$a, prob = 0.05))
abline(v = quantile(data$a, prob = 0.95))
abline(h = 0.05, col = "gray", lty = 2)


set.seed(123)
#trimming estimator with small bandwidth (h = 0.21)
psih.trimmed.sl.h21 = estPsih.trimmed.split(data = subset(data, select = -c(HSQACCWT)),
	h = 0.21,
	a.eval = a.eval,
	a.seq = seq(-4, 8, by = 0.1),
	t = 0.05,
	epsilon = 10^(-2),
	sl.lib.ps = c("SL.earth", "SL.gam.new", "SL.glm", "SL.glm.interaction", 
        "SL.mean", "SL.ranger", "SL.rpart"),
	sl.lib.outcome = c("SL.earth", "SL.gam.new", "SL.glm", "SL.glm.interaction", 
        "SL.mean", "SL.ranger", "SL.rpart"),
	weights = data$HSQACCWT/sum(data$HSQACCWT))
set.seed(123)
#trimming estimator with large bandwidth (h = 0.92)
psih.trimmed.sl.h92 = estPsih.trimmed.split(data = subset(data, select = -c(HSQACCWT)),
	h = 0.92,
	a.eval = a.eval,
	a.seq = seq(-4, 8, by = 0.1),
	t = 0.05,
	epsilon = 10^(-2),
	sl.lib.ps = c("SL.earth", "SL.gam.new", "SL.glm", "SL.glm.interaction", 
        "SL.mean", "SL.ranger", "SL.rpart"),
	sl.lib.outcome = c("SL.earth", "SL.gam.new", "SL.glm", "SL.glm.interaction", 
        "SL.mean", "SL.ranger", "SL.rpart"),
	weights = data$HSQACCWT/sum(data$HSQACCWT))

# PLOT OF TRIMMED AND NON-TRIMMED ESTIMATES

## PLOT USING h = 0.21
plot(data$a, data$y,
	pch = 16, col = "lightgray", cex = 0.25,
	xlab = "", ylab = "", ylim = c(0.5, 12),
	xlim = c(min(a.eval), max(a.eval)))
title(xlab = "Treatment A", line = 2.25, cex.lab = 1.25)
title(ylab = "Estimated Causal Effect", line = 2.25, cex.lab = 1.25)
#point estimate
lines(a.eval, psih.sl.h21$psih.est, col = "blue")
lines(a.eval, psih.trimmed.sl.h21$psih.est1, col = "purple")
#CI
lines(a.eval, psih.sl.h21$psih.ci.l, col = "blue", lty = 2)
lines(a.eval, psih.sl.h21$psih.ci.u, col = "blue", lty = 2)
lines(a.eval, psih.trimmed.sl.h21$psih.ci1.l, col = "purple", lty = 2)
lines(a.eval, psih.trimmed.sl.h21$psih.ci1.u, col = "purple", lty = 2)
legend("bottomright",
	legend = c("No Trimming", "Trimming"),
	lty = 1, col = c("blue", "purple"))

## PLOT USING h = 0.92
plot(data$a, data$y,
	pch = 16, col = "lightgray", cex = 0.25,
	xlab = "", ylab = "", ylim = c(0.5, 12),
	xlim = c(min(a.eval), max(a.eval)))
title(xlab = "Treatment A", line = 2.25, cex.lab = 1.25)
title(ylab = "Estimated Causal Effect", line = 2.25, cex.lab = 1.25)
#point estimate
lines(a.eval, psih.sl.h92$psih.est, col = "blue")
lines(a.eval, psih.trimmed.sl.h92$psih.est1, col = "purple")
#CI
lines(a.eval, psih.sl.h92$psih.ci.l, col = "blue", lty = 2)
lines(a.eval, psih.sl.h92$psih.ci.u, col = "blue", lty = 2)
lines(a.eval, psih.trimmed.sl.h92$psih.ci1.l, col = "purple", lty = 2)
lines(a.eval, psih.trimmed.sl.h92$psih.ci1.u, col = "purple", lty = 2)
legend("bottomleft",
	legend = c("No Trimming", "Trimming"),
	lty = 1, col = c("blue", "purple"))


#It's important for us to understand how the sample
#changes across the dose-response. 

#First, let's understand what the mean covariate is
#across different treatment values.
#This will involve kernel-smoothing, but not trimming.
#For simplicity, let's just focus on age:
age = data$LASTAGE

set.seed(123)
a0.seq = seq(-4, 8, by = 0.1)
a.seq.len = a0.seq[2] - a0.seq[1]
#We'll compute the mean age across a.eval
#for the two choices of bandwidth
age.smoothedMean.h21 = vector(length = length(a.eval))
age.smoothedMean.h92 = vector(length = length(a.eval))

#estimating the PS with SuperLearner
psMod.sl = SuperLearner(Y = data$a,
		X = subset(data, select = -c(a,y,HSQACCWT)),
		SL.library = c("SL.earth", "SL.gam.new", "SL.glm", "SL.glm.interaction", 
        "SL.mean", "SL.ranger", "SL.rpart"),
		obsWeights = data$HSQACCWT/sum(data$HSQACCWT))

#thus, the mean and SD estimates are...
#(note that we need to use the test data here)
meanA.est.sl = as.numeric(psMod.sl$SL.predict)
piRes2 = (data$a - meanA.est.sl)^2
pi2mod <- SuperLearner(Y = piRes2,
			X = subset(data, select = -c(a,y,HSQACCWT)),
        	SL.library = c("SL.earth", "SL.gam.new", "SL.glm", "SL.glm.interaction", 
        "SL.mean", "SL.ranger", "SL.rpart"),
        	obsWeights = data$HSQACCWT/sum(data$HSQACCWT))
pi2mod.vals <- pi2mod$SL.predict

#now compute \hat{\pi}(a0|X) for each a0
#this will again correspond to a matrix,
#where each column is a particular a0
piHat.a0 = sapply(a0.seq, FUN = dnorm, mean = meanA.est.sl, sd = sqrt(pi2mod.vals))
piHat.a = sapply(a.eval, FUN = dnorm, mean = meanA.est.sl, sd = sqrt(pi2mod.vals))

#compute the corresponding St and St.deriv, which are matrices
St.a0 = getSt(piHat.a0, epsilon = 10^(-2), t = 0.05)
St.a = getSt(piHat.a, epsilon = 10^(-2), t = 0.05)
#now we'll compute the smoothed (trimmed) mean age
for(a in a.eval){
	a.index = which(a.eval == a)
	#Compute the inner integral kernel term K_h(a_0 - a)
	kernel.curr.h21 = gaussianKernel(A = data$a, a = a, h = 0.21)
	kernel.curr.h92 = gaussianKernel(A = data$a, a = a, h = 0.92)
	St.a.curr = St.a[,a.index]
	age.smoothedMean.h21[a.index] = sum(age*kernel.curr.h21*St.a.curr)/sum(kernel.curr.h21*St.a.curr)
	age.smoothedMean.h92[a.index] = sum(age*kernel.curr.h92*St.a.curr)/sum(kernel.curr.h92*St.a.curr)
}

#Now let's conduct simple kernel smoothing
#for the ages.
age.kernelSmoothed.h21 = vector(length = length(a.eval))
age.kernelSmoothed.h92 = vector(length = length(a.eval))
for(a in a.eval){
	a.index = which(a.eval == a)
	kernel.curr.h21 = gaussianKernel(A = data$a, a = a, h = 0.21)
	kernel.curr.h92 = gaussianKernel(A = data$a, a = a, h = 0.92)
	age.kernelSmoothed.h21[a.index] = sum(age*kernel.curr.h21)/sum(kernel.curr.h21)
	age.kernelSmoothed.h92[a.index] = sum(age*kernel.curr.h92)/sum(kernel.curr.h92)
}

#smoothed age plot for h = 0.21
plot(data$a, age, xlab = "", ylab = "",
	col = "lightgray", xlim = c(min(a.eval), max(a.eval)),
	ylim = c(min(age), max(age)+10),
	cex = 0.25, pch = 16)
title(xlab = "Treatment A", ylab = "Age",
	line = 2.25, cex.lab = 1.25)
lines(a.eval, age.kernelSmoothed.h21, col = "blue")
lines(a.eval, age.smoothedMean.h21, col = "purple")
legend("topleft",
	legend = c("No Trimming", "Trimming"),
	lty = 1, col = c("blue", "purple"))
#smoothed age plot for h = 0.92
plot(data$a, age, xlab = "", ylab = "",
	col = "lightgray", xlim = c(min(a.eval), max(a.eval)),
	ylim = c(min(age), max(age)+10),
	cex = 0.25, pch = 16)
title(xlab = "Treatment A", ylab = "Age",
	line = 2.25, cex.lab = 1.25)
lines(a.eval, age.kernelSmoothed.h92, col = "blue")
lines(a.eval, age.smoothedMean.h92, col = "purple")
legend("topleft",
	legend = c("No Trimming", "Trimming"),
	lty = 1, col = c("blue", "purple"))



#Below we show how to estimate h for t = 0.05.

#The following function will be useful:
#this function is used to compute the estimators
#if already we have \hat{\pi}(x) and \hat{\mu}(x).
estPsih.trimmed.piMuGiven = function(data,
                                     t, h, epsilon = 10^(-4),
                                     a.eval, a0.seq,
                                     muHat.XA, piHat.XA,
                                     muHat.a0, piHat.a0){
  
  #sample size
  n = nrow(data)
  #This will be a matrix that collects
  #the estimated EIF values for each 
  #observation and treatment combination.
  psih.ifvalues1.num = matrix(nrow = n, ncol = length(a.eval))
  psih.ifvalues1.denom = matrix(nrow = n, ncol = length(a.eval))
  
  #compute smoothed indicator and its derivative.
  St.A = getSt(piHat = piHat.XA, epsilon = epsilon, t = t)
  St.deriv.A = getSt.deriv(piHat = piHat.XA, epsilon = epsilon, t = t)
  St.a0 = getSt(piHat.a0, epsilon = epsilon, t = t)
  St.deriv.a0 = getSt.deriv(piHat.a0, epsilon = epsilon, t = t)
  #Compute IPW term
  ipwRes = (data$y - muHat.XA)/piHat.XA
  
  #when computing the integral manually,
  #we need to know the "gap" between the a0s
  #(in order to compute a Riemann sum)
  a.seq.len = a0.seq[2] - a0.seq[1]
  #Now we just need to compute the kernel terms,
  #which we'll store in matrices.
  #kernel.A will be n * length(a.eval), where each
  #column contains K_h(A) for a particular a.
  kernel.A = sapply(a.eval, gaussianKernel, A = data$a, h = h)
  #kernel.a0 will be length(a.seq) * length(a.eval), where each
  #column contains K_h(a0) for a particular a.
  kernel.a0 = sapply(a.eval, gaussianKernel, A = a0.seq, h = h)
  #normalize the kernels
  kernel.A = t(apply(kernel.A, MARGIN = 1, FUN = function(x) x/colSums(kernel.a0*a.seq.len)))
  kernel.a0 = t(apply(kernel.a0, MARGIN = 1, FUN = function(x) x/colSums(kernel.a0*a.seq.len)))
    
    #First, the estimator will depend
    #on plug-in estimators:
  psihat.pi.num = (a.seq.len*muHat.a0*St.a0)%*%kernel.a0
  psihat.pi.denom = (a.seq.len*St.a0)%*%kernel.a0

    #Now note that, beyond the plug-in estimator parts,
    #\hat{\psi}_num as three terms,
    #and \hat{\psi}_denom as two terms.
    
    #numerator terms one and two
    psihat.num.term1 = kernel.A*ipwRes*St.A
    psihat.num.term2 = kernel.A*muHat.XA*St.deriv.A
    #denominator term one
    psihat.denom.term1 = kernel.A*St.deriv.A
    
    #then, the remaining terms for the numerator and denominator are:
    psihat.num.term3 = -(a.seq.len*muHat.a0*piHat.a0*St.deriv.a0)%*%kernel.a0
    psihat.denom.term2 = -(a.seq.len*piHat.a0*St.deriv.a0)%*%kernel.a0

    #Now compute the psihat.num and psihat.denom
    psihat.num = psihat.pi.num + 
      psihat.num.term1 + psihat.num.term2 + psihat.num.term3
    psihat.denom = psihat.pi.denom +
      psihat.denom.term1 + psihat.denom.term2
    
    #then, the particular estimate is
    #(specifically for the test set)
    psih.ifvalues1.num = psihat.num
    psih.ifvalues1.denom = psihat.denom
  
  #now we'll construct point estimates and confidence intervals
  #We will do this pointwise, for each treatment value.
  
  #the point estimates are
  psih.est1.num = colMeans(psih.ifvalues1.num)
  psih.est1.denom = colMeans(psih.ifvalues1.denom)
  psih.est1 = psih.est1.num/psih.est1.denom
  
  #return the point estimate
  return(psih.est1)
}

estPsih.piMuGiven = function(data,
    h, epsilon = 10^(-4),
    a.eval, a0.seq,
    muHat.XA, piHat.XA, muHat.a0){

    #sample size
    n = nrow(data)
    #This will be a matrix that collects
    #the estimated EIF values for each 
    #observation and treatment combination.
    psih.ifvalues = matrix(nrow = n, ncol = length(a.eval))

    #Compute IPW term
    ipwRes = (data$y - muHat.XA)/piHat.XA

    #when computing the integral manually,
    #we need to know the "gap" between the a0s
    #(in order to compute a Riemann sum)
    a.seq.len = a0.seq[2] - a0.seq[1]

  #Now we just need to compute the kernel terms,
  #which we'll store in matrices.
  #kernel.A will be n * length(a.eval), where each
  #column contains K_h(A) for a particular a.
  kernel.A = sapply(a.eval, gaussianKernel, A = data$a, h = h)
  #kernel.a0 will be length(a.seq) * length(a.eval), where each
  #column contains K_h(a0) for a particular a.
  kernel.a0 = sapply(a.eval, gaussianKernel, A = a0.seq, h = h)
  #normalize the kernels
  kernel.A = t(apply(kernel.A, MARGIN = 1, FUN = function(x) x/colSums(kernel.a0*a.seq.len)))
  kernel.a0 = t(apply(kernel.a0, MARGIN = 1, FUN = function(x) x/colSums(kernel.a0*a.seq.len)))
    
	#The first term of the EIF is:
    psih.term1 = kernel.A*ipwRes

    #The second term then is:
    psih.term2 = a.seq.len*muHat.a0%*%kernel.a0

    #then, the particular estimate is
    psih.ifvalues = psih.term1 + psih.term2
    
    #the point estimates are
    psih.est = colMeans(psih.ifvalues)
    
    #return the point estimate
    return(psih.est)
}

#this function computes the risk for a set of bandwidths
#The data needs to be split into halves.
computeRisk = function(data1, data2,
                       t, h.seq, epsilon = 10^(-4),
                       a.eval, a0.seq, sl.lib = c("SL.glm")){
  
  #the covariate matrices are
  X1 = subset(data1, select = -c(y,a,HSQACCWT))
  X2 = subset(data2, select = -c(y,a,HSQACCWT))
  #the weights are
  weights1 = data1$HSQACCWT/sum(data1$HSQACCWT)
  
  #Estimate propensity scores and outcome regression
  #on the first third.
  #estimate pi
  psMod.sl = SuperLearner(Y = data1$a,
                          X = X1, newX = X2,
                          SL.library = sl.lib,
                          obsWeights = weights1)
  #thus, the mean estimate is (on data2):
  meanA.est = psMod.sl$SL.predict
  #To flexibly estimate the variance, it'll be helpful
  #to first compute squared residuals from this model
  #(this is on data1)
  meanA.est.1 = as.numeric(predict(psMod.sl, newdata = X1)$pred)
  piRes2 = (data1$a - meanA.est.1)^2
  #Then, a flexible model for the variance is:
  pi2mod = SuperLearner(Y = piRes2,
                        X = X1, newX = X2,
                        SL.library = sl.lib,
                        obsWeights = weights1)
  varA.est = pi2mod$SL.predict
  
  #Then, \hat{\pi}(X, A) is (on data2)
  piHat.XA = dnorm(data2$a, mean = meanA.est, sd = sqrt(varA.est))
  
  #meanwhile, the estimated mu(x,a) is...
  muModel = SuperLearner(Y = data1$y,
                         X = subset(data1, select = -c(y, HSQACCWT)),
                         newX = subset(data2, select = -c(y, HSQACCWT)),
                         SL.library = sl.lib,
                         obsWeights = weights1)
  
  #Now compute \hat{\mu}(X, A) (on data2)
  muHat.XA = as.numeric(muModel$SL.predict)
  
  #Now compute \hat{\pi}(a0|X) and \hat{\mu}(X,a0) on data2
  piHat.a0 = sapply(a0.seq, FUN = dnorm, mean = meanA.est, sd = sqrt(varA.est))
  mu.preddata.a0 = X2[rep(seq_len(nrow(X2)), each = length(a0.seq)),]
  #now add the treatment data, corresponding to each a0.seq.
  mu.preddata.a0$a = rep(a0.seq, times = nrow(data2))
  muHat.a0 = matrix(as.numeric(predict(muModel, newdata = mu.preddata.a0)$pred),
                    ncol = length(a0.seq), byrow = TRUE)
  
  #Also compute \hat{\pi}(a|X) and \hat{\mu}(a|X) on data2.
  piHat.a = sapply(a.eval, FUN = dnorm, mean = meanA.est, sd = sqrt(varA.est))
  mu.preddata.a = X2[rep(seq_len(nrow(X2)), each = length(a.eval)),]
  #now add the treatment data, corresponding to each a0.seq.
  mu.preddata.a$a = rep(a.eval, times = nrow(data2))
  muHat.a = matrix(as.numeric(predict(muModel, newdata = mu.preddata.a)$pred),
                     ncol = length(a.eval), byrow = TRUE)
  #We also need to compute smoothed indicator.
  St.a = getSt(piHat = piHat.a, t = t, epsilon = epsilon)

  #Now we'll compute the estimated risk.
  #To do this, we'll need to compute integrals.
  #Note that these integrals are over a.eval !
  a.seq.len = a.eval[2] - a.eval[1]
  
  # THIS PART DEPENDS ON THE PARTICULAR BANDWIDTH
  #vector of risk values across bandwidths
  risk.est = vector(length = length(h.seq))
  #Now, for each h in bandwidth, compute the risk.
  for(h in h.seq){
    h.index = which(h.seq == h)
    
    #Now compute \hat{\psi}_h(a) on data2.
    psiHat.a = estPsih.trimmed.piMuGiven(data = data2,
                                         t = t, h = h, epsilon = epsilon,
                                         a.eval = a.eval, a0.seq = a0.seq,
                                         muHat.XA = muHat.XA, piHat.XA = piHat.XA,
                                         muHat.a0 = muHat.a0, piHat.a0 = piHat.a0)

    #Compute the "denominator" risk
    risk.denom.est = mean(a.seq.len*St.a%*%psiHat.a^2)
    #Compute the "numerator" risk
    risk.num.est = mean((a.seq.len*muHat.a*St.a)%*%psiHat.a)

    #Then, the overall risk estimate is:
    risk.est[h.index] = risk.denom.est - 2*risk.num.est
    print(h)
  }
  
  return(risk.est)
}

library(SuperLearner)
#bandwidths we'll consider
h.seq = seq(0.05, 2, by = 0.01)
#points where treatment is estimated
a.eval = seq(
  quantile(data$a, prob = 0.05),
  quantile(data$a, prob = 0.95),
  length = 10)

#To compute the "plug-in" risk, we just need two datasets:
set.seed(123)
n = nrow(data)
#split indicator
split.ind = sample(rep(1:2, ceiling(n/2))[1:n])
#the datasets are:
data1 = data[split.ind == 1,]
data2 = data[split.ind == 2,]

#then, the risk across bandwidths is:
set.seed(123)
estimatedRisk12 = computeRisk(
	data1 = data1, data2 = data2,
	t = 0.05, h.seq = h.seq, epsilon = 10^(-2),
	a.eval = a.eval,
	a0.seq = seq(-4, 8, by = 0.1),
	sl.lib = c("SL.earth", "SL.gam.new", "SL.glm", "SL.glm.interaction", 
        "SL.mean", "SL.ranger", "SL.rpart"))
set.seed(123)
estimatedRisk21 = computeRisk(
	data1 = data2, data2 = data1,
	t = 0.05, h.seq = h.seq, epsilon = 10^(-2),
	a.eval = a.eval,
	a0.seq = seq(-4, 8, by = 0.1),
	sl.lib = c("SL.earth", "SL.gam.new", "SL.glm", "SL.glm.interaction", 
        "SL.mean", "SL.ranger", "SL.rpart"))

estimatedRisk = (estimatedRisk12 + estimatedRisk21)/2

#plot the estimated risk for the trimming estimator
#across different bandwidth choices.
#This shows that h = 0.92 minimizes the risk
#for the trimming estimator.
plot(h.seq, estimatedRisk,
	xlab = "", ylab = "", type = "l")
title(xlab = "Bandwidth h", ylab = "Estimated Risk",
	line = 2.25, cex.lab = 1.25)
abline(v = 0.92, lty = 2)


#this function computes the risk for the simpler
#non-trimmed estimator \hat{\psi}_h(a).
computeRisk.noTrimming = function(data1, data2,
	h.seq, epsilon = 10^(-2),
    a.eval, a0.seq, sl.lib = c("SL.glm")){

	#the covariate matrices are
	X1 = subset(data1, select = -c(y,a,HSQACCWT))
	X2 = subset(data2, select = -c(y,a,HSQACCWT))
	#the weights are
	weights1 = data1$HSQACCWT/sum(data1$HSQACCWT)

	#Estimate propensity scores and outcome regression
	#on the first third.
	#estimate pi
	psMod.sl = SuperLearner(Y = data1$a,
		X = X1, newX = X2,
		SL.library = sl.lib,
		obsWeights = weights1)
	#thus, the mean estimate is (on data2):
	meanA.est = psMod.sl$SL.predict
	#To flexibly estimate the variance, it'll be helpful
	#to first compute squared residuals from this model
	#(this is on data1)
	meanA.est.1 = as.numeric(predict(psMod.sl, newdata = X1)$pred)
	piRes2 = (data1$a - meanA.est.1)^2
	#Then, a flexible model for the variance is:
	pi2mod = SuperLearner(Y = piRes2,
		X = X1, newX = X2,
		SL.library = sl.lib,
		obsWeights = weights1)
	varA.est = pi2mod$SL.predict

	#Then, \hat{\pi}(X, A) is (on data2)
	piHat.XA = dnorm(data2$a, mean = meanA.est, sd = sqrt(varA.est))

	#meanwhile, the estimated mu(x,a) is...
	muModel = SuperLearner(Y = data1$y,
		X = subset(data1, select = -c(y, HSQACCWT)),
		newX = subset(data2, select = -c(y, HSQACCWT)),
		SL.library = sl.lib,
		obsWeights = weights1)

	#Now compute \hat{\mu}(X, A) (on data2)
	muHat.XA = as.numeric(muModel$SL.predict)

	#Now compute \hat{\pi}(a0|X) and \hat{\mu}(X,a0) on data2
	mu.preddata.a0 = X2[rep(seq_len(nrow(X2)), each = length(a0.seq)),]
	#now add the treatment data, corresponding to each a0.seq.
	mu.preddata.a0$a = rep(a0.seq, times = nrow(data2))
	muHat.a0 = matrix(as.numeric(predict(muModel, newdata = mu.preddata.a0)$pred),
	ncol = length(a0.seq), byrow = TRUE)

	#Now compute \hat{\pi}(a|X) and \hat{\mu}(a|X) on data2.
	mu.preddata.a = X2[rep(seq_len(nrow(X2)), each = length(a.eval)),]
	#now add the treatment data, corresponding to each a0.seq.
	mu.preddata.a$a = rep(a.eval, times = nrow(data2))
	muHat.a = matrix(as.numeric(predict(muModel, newdata = mu.preddata.a)$pred),
		ncol = length(a.eval), byrow = TRUE)

	#Now we'll compute the estimated risk.
	#To do this, we'll need to compute integrals.
	#Note that these integrals are over a.eval !
	a.seq.len = a.eval[2] - a.eval[1]

	# THIS PART DEPENDS ON THE PARTICULAR BANDWIDTH
	#vector of risk values across bandwidths
	risk.est = vector(length = length(h.seq))
	#Now, for each h in bandwidth, compute the risk.
	for(h in h.seq){
		h.index = which(h.seq == h)

		#compute \hat{\psi}_h(a) on data2.
		psiHat.a = estPsih.piMuGiven(data = data2,
			h = h, epsilon = epsilon,
			a.eval = a.eval, a0.seq = a0.seq,
			muHat.XA = muHat.XA, piHat.XA = piHat.XA,
			muHat.a0 = muHat.a0)

		#The first term for the risk is
		risk.term1 = sum(psiHat.a^2*a.seq.len)
		#The second term for the risk is
		risk.term2 = -2*a.seq.len*muHat.a%*%psiHat.a

		#Then, the overall risk estimate is:
		risk.est[h.index] = mean(risk.term1 + risk.term2)
		print(h)
	}

	return(risk.est)
}

set.seed(123)
n = nrow(data)
#split indicator
split.ind = sample(rep(1:2, ceiling(n/2))[1:n])
#the three data subsets are:
data1 = data[split.ind == 1,]
data2 = data[split.ind == 2,]

library(SuperLearner)
h.seq = seq(0.05, 2, by = 0.01)
set.seed(123)
estimatedRiskNoTrimming12 = computeRisk.noTrimming(data1 = data1, data2 = data2,
	h.seq = h.seq, epsilon = 10^(-2),
	a.eval = a.eval,
	a0.seq = seq(-4, 8, by = 0.1),
	sl.lib = c("SL.earth", "SL.gam.new", "SL.glm", "SL.glm.interaction", 
        "SL.mean", "SL.ranger", "SL.rpart"))
set.seed(123)
estimatedRiskNoTrimming21 = computeRisk.noTrimming(data1 = data2, data2 = data1,
	h.seq = h.seq, epsilon = 10^(-2),
	a.eval = a.eval,
	a0.seq = seq(-4, 8, by = 0.1),
	sl.lib = c("SL.earth", "SL.gam.new", "SL.glm", "SL.glm.interaction", 
        "SL.mean", "SL.ranger", "SL.rpart"))

estimatedRiskNoTrimming = (estimatedRiskNoTrimming12 + estimatedRiskNoTrimming21)/2

#plot the estimated risk for the estimator without trimming
#across different bandwidth choices.
#This shows that h = 0.21 minimizes the risk
#for the estimator without trimming.
plot(h.seq, estimatedRiskNoTrimming, type = "l",
	xlab = "Bandwidth h", ylab = "Estimated Risk",
	main = "Risk with No Trimming")
abline(v = 0.21, lty = 2)
h.seq[which.min(estimatedRiskNoTrimming)]

#We can also consider selecting epsilon via entropy:
#This function estimates entropy for a given epsilon.
estEntropy = function(data,
	a.eval, t, epsilon.seq,
	sl.lib = c("SL.earth", "SL.gam", 
	    "SL.glm", "SL.glm.interaction",
	    "SL.mean", "SL.ranger", "SL.rpart")){

	#the covariate matrices are
	X = subset(data, select = -c(y,a,HSQACCWT))
	#the weights are
	weights = data$HSQACCWT/sum(data$HSQACCWT)

	#Estimate propensity scores and outcome regression
	#on the first third.
	#estimate pi
	psMod.sl = SuperLearner(Y = data$a, X = X,
		SL.library = sl.lib,
		obsWeights = weights)
	#thus, the mean estimate is:
	meanA.est = psMod.sl$SL.predict
	#To flexibly estimate the variance, it'll be helpful
	#to first compute squared residuals from this model
	piRes2 = (data$a - meanA.est)^2
	#Then, a flexible model for the variance is:
	pi2mod = SuperLearner(Y = piRes2, X = X,
		SL.library = sl.lib,
		obsWeights = weights)
	varA.est = pi2mod$SL.predict

	#Now compute \hat{\pi}(a|X)
  	piHat.a = sapply(a.eval, FUN = dnorm, mean = meanA.est, sd = sqrt(varA.est))
  	
    #now compute entropy across epsilon values
    entropy.est = vector(length = length(epsilon.seq))
    for(epsilon in epsilon.seq){
    	ep.index = which(epsilon.seq == epsilon)
    	#compute smoothed indicator
		St.a = getSt(piHat = piHat.a, t = t, epsilon = epsilon)
		#then, the two terms for the entropy are:
		H1 = St.a*log(St.a, base = 2)
		H2 = (1-St.a)*log(1-St.a, base = 2)
		#then, the entropy is:
		entropy.curr = -(H1 + H2)
		#There may be extreme values when S = 0 or S = 1.
		#In this case, the estimated entropy should be 0.
		entropy.curr = ifelse(St.a == 0 | St.a == 1, 0, entropy.curr)
		#Then, the estimated entropy is:
		entropy.est[ep.index] = mean(colMeans(entropy.curr))
	}

    return(entropy.est)
}

set.seed(123)
n = nrow(data)
#split indicator
split.ind = sample(rep(1:2, ceiling(n/2))[1:n])
#the three data subsets are:
data1 = data[split.ind == 1,]
data2 = data[split.ind == 2,]

#we'll consider the follow epsilon values:
set.seed(123)
epsilon.seq = c(10^(-seq(1, 5, by = 0.2)))
entropy = estEntropy(data = data,
	a.eval = a.eval, t = 0.05,
	epsilon.seq = epsilon.seq,
	sl.lib = c("SL.earth", "SL.gam.new", "SL.glm", "SL.glm.interaction", 
        "SL.mean", "SL.ranger", "SL.rpart"))
#Plot the entropy across log_{10}(\epsilon) values.
plot(log(epsilon.seq, base = 10), entropy,
	xlab = "", ylab = "", pch = 16)
title(xlab = expression(log[10](epsilon)),
	ylab = "Estimated Entropy", line = 2.25, cex.lab = 1.25)
abline(h = 0.05, lty = 2)
