require(mice)
require(mosaic)
require(dplyr)
require(tidyr)
require(magrittr)

#' generateComplete
#'
#' This function generate complete datasets.
#' @param outcomeType Is the outcome variable binary or continuous?
#' @param numSamp Number of observations in the sample
#' @param beta0 Intercept
#' @param beta1 Treatment effect
#' @param beta2 Slope of auxiliary variable
#' @param xmean Mean of auxiliary variable
#' @param xsd Standard deviation of auxiliary variable
#' @param errorsd Parameter to add noise. Standard deviation of errors

generateComplete <- function(outcomeType, numSamp, beta0,
                             beta1, beta2, xmean, xsd,
                             errorsd) {

  x <- rnorm(numSamp, mean = xmean, sd = xsd) #Continuous variable
  t <- c(rep(x = -1, times = numSamp/2), rep(x = 1, times = numSamp/2)) #Treatment variable. Half assigned to treatment group (1)
  #Half assigned to control group (-1)
  error <- rnorm(numSamp, mean = 0, sd = errorsd) #Error term

  if(outcomeType == "binary") {
    m <- beta0 + beta1 * t + beta2 * x + error #Relationship between the explanatory variables (x & t) and response variable (y)
    sigmam = exp(m)/(1 + exp(m)) #Inverse logistic function to get back probability
    y <- (runif(numSamp) < sigmam) #Randomly determine value based on probability calculated in previous step
  }

  else if(outcomeType == "continuous") {
    y <- beta0 + beta1 * t + beta2 * x + error
  }

  else{
    print("Not a valid outcome type")
    return()
  }

  df <- data.frame(cbind(y = y, x = x, t = t))
  return(df)
}

#' makeMissing
#'
#' This function take a complete dataset and generates missingness
#' @param outcomeType Is the outcome variable binary or continuous?
#' @param mechanism What is the missing data mechanism?
#' @param numSamp Number of observations in the sample
#' @param propMiss Expected proportion of missingness
#' @param rr Risk ratio
#' @param df The complete dataframe
#' @param beta0 Intercept
#' @param beta1 Treatment effect
#' @param beta2 Slope of auxiliary variable
#' @param xmean Mean of auxiliary variable

makeMissing <- function(outcomeType, mechanism, numSamp, propMiss, rr, df, beta0, beta1, beta2, xmean) {
  #Regardless of missing data mechanism want to achieve an overall proportion of missing equal to p
  if(mechanism == "mcar") { #Idea is to use a logistic regression model to determine missingness
    g0 <- log(propMiss/(1-propMiss)) #Three different cases for the different missing data mechanisms
    g1 <- 0
    g2 <- 0
  }

  else if(mechanism == "mar") {
    z <- 2*propMiss/(1+rr)
    g0 <- log(rr*z^2/((1-z)*(1-rr*z)))/2
    g1 <- log(rr*(1-z)/(1-rr*z))/2
    g2 <- 0
  }

  else if(mechanism == "mnar"){

    if(outcomeType == "binary") {
      w <- prop(~y, df)
      z <- propMiss/(1-w + rr*w)
    }

    else if(outcomeType == "continuous") {
      z <- propMiss/(.5 + rr*.5)
    }

    else{
      print("Not a valid outcome type")
      return()
    }

    g0 <- log(z/(1-z))
    g1 <- 0
    g2 <- log(rr*(1-z)/(1-rr*z))
  }

  else {
    print("Not a valid missing data mechanism")
    return()
  }

  if(outcomeType == "binary") {
    m <- g0 + g1*df$t + g2*df$y
  }

  else if(outcomeType == "continuous") {
    m <- g0 + g1*df$t + g2*((df$y - beta0 + beta1 - beta2*xmean)/(2*beta1))
  }

  sigmam <- exp(m) / (1+ exp(m))

  temp <- runif(numSamp)
  r <- ifelse(sigmam > temp, 0, 1)
  df$r <- r
  df <- df %>% mutate(obsy = ifelse(r == 0, NA, y))
  return(df)
}

#' generateData
#'
#' This function creates a complete dataset and then generates missingness
#' @param mechanism What is the missing data mechanism?
#' @param outcomeType Is the outcome variable binary or continuous?
#' @param propMiss Expected proportion of missingness. Defaults to 0.1.
#' @param numSamp Number of observations in the sample. Defaults to 1000.
#' @param beta0 Intercept. Defaults to 1.
#' @param beta1 Treatment effect. Defaults to 2.
#' @param beta2 Slope of auxiliary variable. Defaults to 1.
#' @param xmean Mean of auxiliary variable. Defaults to 0.
#' @param xsd Standard deviation of auxiliary variable. Defaults to 1.
#' @param errorsd Standard deviation of error term. Defaults to 1.
#' @param rr Risk ratio. Defaults to 5. Defaults to 5.

generateData <- function(mechanism, outcomeType, propMiss = .1, numSamp = 1000,
                         beta0 = 1, beta1 = 2, beta2 = 1, xmean = 0, xsd = 1,
                         errorsd = .1, rr = 5) {

  df <- generateComplete(outcomeType, numSamp, beta0, #Creates the data
                         beta1, beta2, xmean, xsd, errorsd)

  df <- makeMissing(outcomeType = outcomeType, mechanism = mechanism,
                    numSamp = numSamp, df = df, propMiss = propMiss, rr = rr,
                    beta0 = beta0, beta1 = beta1, beta2 = beta2, xmean = xmean) #Creates missingness
}

#' runImputation
#'
#' This function creates a models via multiple imputation and returns model parameters
#' @param df The dataset with missing values generated by generateData()
#' @param m Number of imputations.
#' @param method The method of imputing univariate missingness
#' @param interaction Want to add interaction to the imputation model? Defaults to FALSE.
#' @param donors Number of donors if using predictive mean matching. Defaults to 3.
#' @param outcomeType Is the outcome variable binary or continuous?

runImputation <- function(df, m = 5, method, interaction = FALSE, donors = 3, outcomeType) {
  if(method == "logreg" || method == "logreg.boot") {
    df$obsy <- as.factor(df$obsy)
  }

  if(interaction == TRUE) {
    df <- df %>% mutate(tx = (t-mean(t)) * (x-mean(x)))
    dummy <- mice(df, method = method, maxit = 0)
    meth <- dummy$meth
    meth["tx"] <- paste("~I(", "(t - mean(t)) * (x - mean(x)))", sep = "")
    pred <- dummy$pred
    pred["obsy",] <- c(0,1,1,0,0,1)
    imp <- mice(df, method = method, m = m, printFlag = FALSE, predictorMatrix = pred)
  }

  else {
    pred <- matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0), nrow = 5, ncol = 5, byrow = T)
    if(method == "pmm") {
      imp <- mice(df, method = method, m = m, printFlag = FALSE, predictorMatrix = pred, donors = donors)
    }
    else {
      imp <- mice(df, method = method, m = m, printFlag = FALSE, predictorMatrix = pred)
    }
  }

  if(outcomeType == "binary") {
    mods <- with(imp, glm(obsy ~ t + x, family = binomial(logit)))
  }

  else if(outcomeType == "continuous") {
    mods <- with(imp, lm(obsy ~ t + x))
  }

  info <- summary(pool(mods))
  return(as.vector(t(info[, c("est", "se", "lo 95", "hi 95")])))
}

#' runTrue
#'
#' This function creates a models using the true values of y
#' @param df The dataset with missing values generated by generateData()
#' @param outcomeType Is the outcome variable binary or continuous?

runTrue <- function(df, outcomeType) {
  if(outcomeType == "binary") {
    mod <- glm(y ~ t + x, data = df, family = binomial(logit))
  }

  else if(outcomeType == "continuous") {
    mod <- lm(y ~ t + x, data = df)
  }

  else  {
    print("Outcome type is not valid")
    return()
  }

  info <- summary(mod)$coefficients[, c("Estimate", "Std. Error")]
  info <- cbind(info, confint.default(mod))
  return(as.vector(t(info)))
}

#' runCC
#'
#' This function creates a models using complete case analysis
#' @param df The dataset with missing values generated by generateData()
#' @param outcomeType Is the outcome variable binary or continuous?

runCC <- function(df, outcomeType) {
  if(outcomeType == "binary") {
    mod <- glm(obsy ~ t + x, df, na.action = na.exclude, family = binomial(logit))
  }

  else if(outcomeType == "continuous") {
    mod <- lm(obsy ~ t + x, df, na.action = na.exclude)
  }

  else  {
    print("Outcome type is not valid")
    return()
  }

  info <- summary(mod)$coefficients[, c("Estimate", "Std. Error")]
  info <- cbind(info, confint.default(mod))
  return(as.vector(t(info)))
}

#' #' runSimulation
#'
#' This function runs a simulations fitting models using various methods of handling missing data. Passes back an array of model parameter values.
#' @param numSim Number of simulations. Defaults to 1000.
#' @param m Number of imputations. Defaults to 5.
#' @param mechanism What is the missing data mechanism? Defaults to mcar.
#' @param propMiss Expected proportion of missingness. Defaults to 0.1.
#' @param numSamp Number of observations in the sample. Defaults to 1000.
#' @param beta0 Intercept. Defaults to 1.
#' @param beta1 Treatment effect. Defaults to 2.
#' @param beta2 Slope of auxiliary variable. Defaults to 1.
#' @param xmean Mean of auxiliary variable. Defaults to 0.
#' @param xsd Standard deviation of auxiliary variable. Defaults to 1.
#' @param errorsd Standard deviation of error term. Defaults to 0.5.
#' @param rr Risk ratio. Defaults to 5. Defaults to 5.
#' @param outcomeType Is the outcome variable binary or continuous?

runSimulation <- function(numSim = 1000, m = 5, mechanism = "mcar",
                          propMiss = .1, numSamp = 1000, beta0 = 1,
                          beta1 = 2, beta2 = 1, xmean = 0,
                          xsd = 1, errorsd = .5, rr = 5, outcomeType) {

  tab <- array(NA, dim = c(8, numSim, 12))
  dimnames(tab) <- list(c("Truth", "Complete Case", "Bayesian MI", "Bootstrap MI", "PMM (d = 1)",
                          "PMM (d = 3)", "PMM (d = 10)", "PMM (d = 20)"),
                        as.character(1:numSim),
                        c("beta0", "se0", "lowci0", "highci0",
                          "beta1", "se1", "lowci1", "highci1",
                          "beta2", "se2", "lowci2", "highci2"))

  for (i in 1:numSim){
    df <- generateData(outcomeType = outcomeType, mechanism = mechanism, propMiss = propMiss, numSamp = numSamp,
                       beta0 = beta0, beta1 = beta1, beta2 = beta2, xmean = xmean,
                       xsd = xsd, errorsd = errorsd, rr = rr)

    tab[1, i,] <- runTrue(df = df, outcomeType = outcomeType)
    tab[2, i,] <- runCC(df = df, outcomeType = outcomeType)

    if(outcomeType == "binary") {
      meth1 = "logreg"
      meth2 = "logreg.boot"
    }
    else{
      meth1 = "norm"
      meth2 = "norm.boot"
    }
    tab[3, i,] <- runImputation(df = df, m = m, method = meth1, outcomeType = outcomeType)
    tab[4, i,] <- runImputation(df = df, m = m, method = meth2, outcomeType = outcomeType)

    tab[5, i,] <- runImputation(df = df, m = m, method = "pmm", donors = 1, outcomeType = outcomeType)
    tab[6, i,] <- runImputation(df = df, m = m, method = "pmm", donors = 3, outcomeType = outcomeType)
    tab[7, i,] <- runImputation(df = df, m = m, method = "pmm", donors = 10, outcomeType = outcomeType)
    tab[8, i,] <- runImputation(df = df, m = m, method = "pmm", donors = 20, outcomeType = outcomeType)
  }
  return(tab)
}
