
# Packages ----------------------------------------------------------------

require(tidyverse)
require(magrittr)
require(R2jags)
require(bayesplot)
require(TeachingDemos)
require(factoextra)
require(highcharter)
require(dplyr)

# Data --------------------------------------------------------------------

dat <- read_csv("./heart.csv")
dat <- unique(dat) # remove any duplicate presents

# Describing the features  ------------------------------------------------

summary(dat)

# analyzing the qualitative and quantitative data

# quantitative data
print_column_chart <- function(name, title){
  keeps <- c(name, 'output')
  frame = dat[ , (names(dat) %in% keeps)]
  
  return(hchart(frame[, (names(frame) %in% name)], type = "column") %>%
           hc_title(text= title) %>%
           #hc_add_theme(hc_theme_google()) %>%
           hc_xAxis(title = name) %>%
           hc_chart(options3d=list(enabled=TRUE, alpha=2, beta=-10, 
                                   depth=100, viewDistance=25)) %>% 
           hc_plotOptions(column=list(depth= 100)))
  
}

print_column_chart('age', 'Persons Age')

# prefer removing less significance data and not informative (in my opinion variables)

# for visualization
res.pca <- prcomp(dat, scale = FALSE) # compute PCA
fviz_eig(res.pca) # Show the percentage of variances explained by each principal component

# patients  with a similar profile are grouped together
fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,     # Avoid text overlapping
             label=FALSE
)

fitted <- glm(y ~ x2+x3+x4+x8+x9+x10+x11+x12+x13, family = binomial(),data=dat[, -c(1, 5, 6, 7, 14)]); summary(fitted)

hchart(cor(dat))

# Load JAGS data \w Logistic Regression model ----------------------------------------------------------

# Writing model for jags
model <- function(){
  # Likelihood
  for (i in 1:N){
    y[i] ~ dbern(p[i])
    logit(p[i]) <-  beta2*x2[i] + beta3*x3[i] + beta4*x4[i] + beta8*x8[i] + beta9*x9[i] + beta10*x10[i] + beta11*x11[i] + beta12*x12[i] + beta13*x13[i] # omit intercept and not useful variables
  }
  
  # Defining the priors
  beta2 ~ dnorm(-1, 1.0E-6)
  beta3 ~ dnorm(0.5, 1.0E-6)
  beta4 ~ dnorm(0, 1.0E-6)
  beta8 ~ dnorm(0, 1.0E-6)
  beta9 ~ dnorm(-0.5, 1.0E-6)
  beta10 ~ dnorm(0.5, 1.0E-6)
  beta11 ~ dnorm(0.5, 1.0E-6)
  beta12 ~ dnorm(-0.5, 1.0E-6)
  beta13 ~ dnorm(-0.5, 1.0E-6)
}

# Preparing data for JAGS
N <- nrow(dat) # number of obs

# Dependent Variable
y <- as.vector(dat$output) # Target variable, risks heart attack?

# Independent Variables
# x1 <- as.vector(dat$age) # Age
x2 <- as.vector(dat$sex) # Sex
x3 <- as.vector(dat$cp) # Chest Pain Type: 1) Typical Angina 2) Atypical Angina 3) Non-Anginal Pain 4) Asymptomatic
x4 <- as.vector(dat$trtbps) # resting blood pressure (in mm Hg)
# x5 <- as.vector(dat$chol) # cholestoral in mg/dl fetched via BMI sensor
# x6 <- as.vector(dat$fbs) # (fasting blood sugar > 120 mg/dl) (1 = true; 0 = false)
# x7 <- as.vector(dat$restecg) # resting electrocardiographic results 0) Normal 1) having ST-T wave abnormality (T wave inversions and/or ST elevation or depression of > 0.05 mV) 2) showing probable or definite left ventricular hypertrophy by Estes' criteria 
x8 <- as.vector(dat$thalachh) # maximum heart rate achieved
x9 <- as.vector(dat$exng) # exercise induced angina (1 = yes; 0 = no)
x10 <- as.vector(dat$oldpeak) # Previous peak
x11 <- as.vector(dat$slp) # Slope
x12 <- as.vector(dat$caa) # number of major vessels (0-3)
x13 <- as.vector(dat$thall) # Thal rate

data.jags <- list("y" = y, "N" = N,
                  "x2" = x2, "x3" = x3, "x4" = x4, "x8" = x8, "x9" = x9, "x10" = x10, "x11" = x11, "x12" = x12, "x13" = x13)

# Defining parameters of interest
mod.params <- c("beta2", "beta3", "beta4", "beta8", "beta9", "beta10", "beta11", "beta12", "beta13")

# Run JAGS
# set.seed(123)
n.chains <- 3
mod.fit <- jags(data = data.jags,                                               # DATA
                model.file = model,                                             # MODEL
                parameters.to.save = mod.params,                                # TRACKING
                n.chains = n.chains, n.iter = 10000, n.burnin = 1000, n.thin=10)       # MCMC

# Results -------------------------------------------------

mod.fit 
mod.fit$BUGSoutput$summary # Rhat is a measure of the convergence through these 3 sims of MC, equal 1 is pretty good

# Diagnostic
traceplot(mod.fit)

# We can get better diagnostics with other packages

# Plots with BayesPlot
chainArray <- mod.fit$BUGSoutput$sims.array

# considering to split each couples of parameters
bayesplot::mcmc_combo(chainArray, pars = c("deviance", "beta2", "beta3"))
bayesplot::mcmc_combo(chainArray, pars = c("beta4", "beta8", "beta9", "beta10"))
bayesplot::mcmc_combo(chainArray, pars = c("beta11", "beta12", "beta13"))

bayesplot::mcmc_acf(chainArray, pars = c("deviance", "beta2", "beta3"))
bayesplot::mcmc_acf(chainArray, pars = c("beta4", "beta8", "beta9", "beta10"))
bayesplot::mcmc_acf(chainArray, pars = c("beta11", "beta12", "beta13"))

autocorr.diag(as.mcmc(mod.fit))

# Diagnostic with coda
coda.fit <- coda::as.mcmc(mod.fit)

coda::acfplot(coda.fit)
coda::densplot(coda.fit)

coda::raftery.diag(coda.fit)

coda::geweke.diag(coda.fit)
coda::geweke.plot(coda.fit)

coda::gelman.diag(coda.fit)
coda::gelman.plot(coda.fit)

coda::heidel.diag(coda.fit)


# Cumeans -----------------------------------------------------------------

df <- as.data.frame(mod.fit$BUGSoutput$sims.array)
for (i in seq(1, length(df), by=3)) {
  param <- unlist(str_split(colnames(df)[i], "\\."))[[2]]
  
  hc <- df %>%
    hchart('line', hcaes(x = 1:nrow(df), y = cummean(df[, i])), color = "red", name = "First Chain") %>%
      hc_add_series( cummean(df[, i+1]), type = "line", color = "blue", name = "Second Chain") %>%
        hc_add_series( cummean(df[, i+2]), type = "line", color = "green", name = "Thrid Chain") %>%
          hc_title(text = paste("The Empirical Means of ", param, sep="")) %>%
            hc_xAxis(title = list(text = "Iteration")) %>%
              hc_yAxis(title = list(text = 'Cumulative Mean'))
  print(hc)
}

# Approximation errors ----------------------------------------------------

mcse_dataframe <- data.frame(MCSE_Chain1 = rep(NA, length(colnames(mod.fit$BUGSoutput$sims.matrix))), MCSE_Chain2 = rep(NA, length(colnames(mod.fit$BUGSoutput$sims.matrix))), MCSE_Chain3 = rep(NA, length(colnames(mod.fit$BUGSoutput$sims.matrix))), MCSE_mean = rep(NA, length(colnames(mod.fit$BUGSoutput$sims.matrix))))
rownames(mcse_dataframe) <- colnames(mod.fit$BUGSoutput$sims.matrix)
for(colname in colnames(mod.fit$BUGSoutput$sims.matrix)){
  for(j in 1:n.chains)
    mcse_dataframe[colname, j] <- LaplacesDemon::MCSE(mod.fit$BUGSoutput$sims.array[, j, colname])
  mcse_dataframe[colname, "MCSE_mean"] <- mean(unlist(mcse_dataframe[colname, c(1:3)]))
}

mcse_dataframe

# Joining chains --------------------------------------------------------------

chainMat <- mod.fit$BUGSoutput$sims.matrix

# Point estimates
(p.hat.jags <- colMeans(chainMat))

# Intervals
cred <- 0.95
(p.ET.jags <- apply(chainMat, 2, quantile, 
                    prob=c((1-cred)/2, 1-(1-cred)/2)))

# What about the HPD?
(p.HPD.jags <- coda::HPDinterval(as.mcmc(chainMat)))


# Second Model ------------------------------------------------------------

# Writing model for jags
model2 <- function(){
  # Likelihood
  for (i in 1:N){
    y[i] ~ dbern(p[i])
    cloglog(p[i]) <-  beta2*x2[i] + beta3*x3[i] + beta4*x4[i] + beta8*x8[i] + beta9*x9[i] + beta10*x10[i] + beta11*x11[i] + beta12*x12[i] + beta13*x13[i] # omit intercept and not useful variables
  }
  
  # Defining the priors
  beta2 ~ dnorm(-1, 1.0E-6)
  beta3 ~ dnorm(0.5, 1.0E-6)
  beta4 ~ dnorm(0, 1.0E-6)
  beta8 ~ dnorm(0, 1.0E-6)
  beta9 ~ dnorm(-0.5, 1.0E-6)
  beta10 ~ dnorm(0.5, 1.0E-6)
  beta11 ~ dnorm(0.5, 1.0E-6)
  beta12 ~ dnorm(-0.5, 1.0E-6)
  beta13 ~ dnorm(-0.5, 1.0E-6)
}

# Preparing data for JAGS
N <- nrow(dat) # number of obs

# Dependent Variable
y <- as.vector(dat$output) # Target variable, risks heart attack?

# Independent Variables
# x1 <- as.vector(dat$age) # Age
x2 <- as.vector(dat$sex) # Sex
x3 <- as.vector(dat$cp) # Chest Pain Type: 1) Typical Angina 2) Atypical Angina 3) Non-Anginal Pain 4) Asymptomatic
x4 <- as.vector(dat$trtbps) # resting blood pressure (in mm Hg)
# x5 <- as.vector(dat$chol) # cholestoral in mg/dl fetched via BMI sensor
# x6 <- as.vector(dat$fbs) # (fasting blood sugar > 120 mg/dl) (1 = true; 0 = false)
# x7 <- as.vector(dat$restecg) # resting electrocardiographic results 0) Normal 1) having ST-T wave abnormality (T wave inversions and/or ST elevation or depression of > 0.05 mV) 2) showing probable or definite left ventricular hypertrophy by Estes' criteria 
x8 <- as.vector(dat$thalachh) # maximum heart rate achieved
x9 <- as.vector(dat$exng) # exercise induced angina (1 = yes; 0 = no)
x10 <- as.vector(dat$oldpeak) # Previous peak
x11 <- as.vector(dat$slp) # Slope
x12 <- as.vector(dat$caa) # number of major vessels (0-3)
x13 <- as.vector(dat$thall) # Thal rate

data.jags <- list("y" = y, "N" = N,
                  "x2" = x2, "x3" = x3, "x4" = x4, "x8" = x8, "x9" = x9, "x10" = x10, "x11" = x11, "x12" = x12, "x13" = x13)

# Defining parameters of interest
mod.params <- c("beta2", "beta3", "beta4", "beta8", "beta9", "beta10", "beta11", "beta12", "beta13")

# Run JAGS
# set.seed(123)
n.chains <- 3
mod.fit2 <- jags(data = data.jags,                                               # DATA
                model.file = model,                                             # MODEL
                parameters.to.save = mod.params,                                # TRACKING
                n.chains = n.chains, n.iter = 10000, n.burnin = 1000, n.thin=10)       # MCMC

mod.fit2

# First Model \w prediction ------------------------------------------------------------

new_record <- list(x2 = sample(x2, 1, replace = T), x3 = sample(x3, 1, replace = T), x4 = sample(x4, 1, replace = T), x8 = sample(x8, 1, replace = T), x9 = sample(x9, 1, replace = T),
                   x10 = sample(x10, 1, replace = T), x11 = sample(x11, 1, replace = T), x12 = sample(x12, 1, replace = T), x13 = sample(x13, 1, replace = T))

x <- mod.fit$BUGSoutput$sims.matrix[,"beta2"]*new_record$x2 + mod.fit$BUGSoutput$sims.matrix[,"beta3"]*new_record$x3 + mod.fit$BUGSoutput$sims.matrix[,"beta4"]*new_record$x4 + mod.fit$BUGSoutput$sims.matrix[,"beta8"]*new_record$x8 + mod.fit$BUGSoutput$sims.matrix[,"beta9"]*new_record$x9 + mod.fit$BUGSoutput$sims.matrix[,"beta10"]*new_record$x10 + mod.fit$BUGSoutput$sims.matrix[,"beta11"]*new_record$x11 + mod.fit$BUGSoutput$sims.matrix[,"beta12"]*new_record$x12 + mod.fit$BUGSoutput$sims.matrix[,"beta13"]*new_record$x13	
x_mu <- 1/(1+exp(-x))

y_pred <- unlist(lapply(x_mu, function(x) rbinom(n = 1, size = 1, prob = x)))
sd(y_pred)
prop.table(table(y_pred))

summary(y_pred)[4]
hchart(y_pred, type = "column", name = "model1", color = randomcoloR::randomColor()) %>%
  hc_title(text = "Predictions with model 1") %>%
  hc_xAxis(title = "model1") %>%
  hc_chart(options3d=list(enabled=TRUE, alpha=2, beta=-10, 
                          depth=100, viewDistance=25)) %>% 
  hc_plotOptions(column=list(depth= 100))

# second model fit
z <- mod.fit2$BUGSoutput$sims.matrix[,"beta2"]*new_record$x2 + mod.fit2$BUGSoutput$sims.matrix[,"beta3"]*new_record$x3 + mod.fit2$BUGSoutput$sims.matrix[,"beta4"]*new_record$x4 + mod.fit2$BUGSoutput$sims.matrix[,"beta8"]*new_record$x8 + mod.fit2$BUGSoutput$sims.matrix[,"beta9"]*new_record$x9 + mod.fit2$BUGSoutput$sims.matrix[,"beta10"]*new_record$x10 + mod.fit2$BUGSoutput$sims.matrix[,"beta11"]*new_record$x11 + mod.fit2$BUGSoutput$sims.matrix[,"beta12"]*new_record$x12 + mod.fit2$BUGSoutput$sims.matrix[,"beta13"]*new_record$x13	
z_mu <- 1/(1+exp(-z)) 

z_pred <- unlist(lapply(z_mu, function(z) rbinom(n = 1, size = 1, prob = z)))
sd(z_pred)
prop.table(table(z_pred))

summary(z_pred)[4]
hchart(z_pred, type = "column", name = "model2", color = randomcoloR::randomColor()) %>%
  hc_title(text = "Predictions with model 2") %>%
  hc_xAxis(title = "model2") %>%
  hc_chart(options3d=list(enabled=TRUE, alpha=2, beta=-10, 
                          depth=100, viewDistance=25)) %>% 
  hc_plotOptions(column=list(depth= 100))

## better the first model, is more concetrated

hchart( density(y_pred), type = "area", color = randomcoloR::randomColor(), name = "Model1 - logit") %>%
  hc_add_series( density(z_pred), type = "area", color = randomcoloR::randomColor(), name = "Model2 - cloglog") %>%
  hc_title(text = "Densities") %>%
  hc_xAxis(title = "Comparison")

data.frame(sd_Ypred = sd(y_pred), sd_Zpred = sd(z_pred))
data.frame(mean_Ypred = mean(y_pred), mean_Zpred = mean(z_pred))

prop.table(table(y))
prop.table(table(y_pred))
prop.table(table(z_pred))

# Bring the report --------------------------------------------------------

library(ggmcmc)
S <- ggs(as.mcmc(mod.fit))
ggmcmc(S)

