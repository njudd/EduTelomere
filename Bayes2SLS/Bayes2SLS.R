# Nicholas Judd
# 2024-10-19



# trying a bunch of different approaches for fuzzy RD

pacman::p_load(tidyverse, rstan, causaldata, fixest, modelsummary)


vet <- causaldata::mortgages


# Create an "above-cutoff" variable as the instrument
vet <- vet %>% mutate(above = qob_minus_kw > 0)

vet <- vet %>%  filter(abs(qob_minus_kw) < 12)


# DV ~ running*treated

# treated ~ running


stage1 <- lm(above ~ qob_minus_kw, data = vet)

stage2 <- lm(home_ownership ~ stage1$fitted.values*qob_minus_kw, data = vet)


m <- feols(home_ownership ~
             nonwhite  | # Control for race
             bpl + qob | # fixed effect controls
             qob_minus_kw*vet_wwko ~ # Instrument our standard RDD
             qob_minus_kw*above, # with being above the cutoff
           se = 'hetero', # heteroskedasticity-robust SEs
           data = vet) 

# And look at the results
summary(m)

m_simple <- feols(home_ownership ~
             qob | # fixed effect controls
             qob_minus_kw*vet_wwko ~ # Instrument our standard RDD
             qob_minus_kw*above, # with being above the cutoff
           se = 'hetero', # heteroskedasticity-robust SEs
           data = vet) 


summary(m_simple)




# stan links from MEA
# https://mc-stan.org/users/documentation/
# https://bruno.nicenboim.me/bayescogsci/ch-introstan.html#stan-syntax




# https://khakieconomics.github.io/2017/11/26/Bayesian_iv.html

# https://rpsychologist.com/adherence-analysis-IV-brms


# https://david-salazar.github.io/posts/bayesian-statistics/2020-06-03-bayesian-instrumental-variable-regression.html




# https://www.r-bloggers.com/2014/01/instrumental-variables-simulation/
library(MASS)
# we are really generating x* and c and using a common variance
xStarAndC <- mvrnorm(1000, c(20, 15), matrix(c(1, 0.5, 0.5, 1), 2, 2))
xStar <- xStarAndC[, 1]
c <- xStarAndC[, 2]

z <- rnorm(1000)
x <- xStar + z

# now lets simulate the response var
# using 1 makes it easy to estimate how 'wrong' an estimator is and toss
# some noise on y
y <- 1 + x + c + rnorm(1000, 0, 0.5)



cor(x, c); cor(z, c)


lm(y ~ x + c)

lm(y ~ x) # incorrectly estimated...


# now lets use the IV estimator
xHat <- lm(x ~ z)$fitted.values
lm(y ~ xHat)

















set.seed(42)
# Load the stan library
library(rstan)
options(mc.cores = parallel::detectCores())

# Compile the model
compiled_model <- stan_model("~/projects/EduTelomere/Bayes2SLS/classic_iv_hierarchical_prior.stan")

# Let's make some fake data

N <- 1000 # Number of observations
PX <- 1 # Number of exogenous variables
PZ <- 1 # Number of instruments
J <- 10 # Number of previous studies

# Previous study parameters, drawn from the prior
beta_hat <- rnorm(1, 0, 1)
sigma_beta <- truncnorm::rtruncnorm(1, a = 0)
beta_j <- rnorm(J, beta_hat, sigma_beta)
se_j <- truncnorm::rtruncnorm(J, a = 0)
b_j <- rnorm(J, beta_j, se_j)

# Exogenous variables (make them correlated with a random correlation matrix)
X_exog <- MASS::mvrnorm(N, rep(0, PX), cor(matrix(rnorm(PX*(PX+5)), PX+5, PX)))

# We have to feed it some dummy endogenous variables but these won't make a difference
X_endog <- rnorm(N)

# Some fake instruments
Z <- MASS::mvrnorm(N, rep(0, PZ), cor(matrix(rnorm(PZ*(PZ+5)), PZ+5, PZ)))

Y_outcome <- rnorm(N)

data_list <- list(N = N, PX = PX, PZ = PZ, J = J,
                  b_j = b_j, se_j = se_j, X_exog = X_exog, X_endog = X_endog, 
                  Z = Z, Y_outcome = Y_outcome, 
                  run_estimation = 0)




draws_from_model <- sampling(compiled_model, data = data_list, iter = 50, chains = 1)



# Let's use the next to last non-warm-up draw as our fake data (draw 24)
y_sim <- extract(draws_from_model, pars = "y_sim")[[1]][24,]
x_endog <- extract(draws_from_model, pars = "x_endog")[[1]][24,]
true_beta <- extract(draws_from_model, pars = "beta_ours")[[1]][24]

# Now let's make a new data list
data_list_2 <- data_list
data_list_2$X_endog <- x_endog
data_list_2$Y_outcome <- y_sim
data_list_2$run_estimation <- 1






























































# https://gist.github.com/rmcelreath/87ab316cfb1be10fb1057c47b20317d3


pacman::p_load(rethinking, dagitty, shape, brms)

# the structural model
# W: wages
# Q: quarter of year of birth
# E: education
# U: unobserved confound

the_dag <- dagitty("dag{
    Q -> E -> W
    E <- U -> W
}")
coordinates(the_dag) <- list(x=c(Q=0,E=1,W=2,U=1.5),y=c(U=0,Q=1,E=1,W=1))
drawdag(the_dag)

# find instruments
instrumentalVariables( the_dag , exposure="E" , outcome="W" )

# simulate
# bEW: effect of education on wages

bEW <- 0

set.seed(73)
N <- 500
U_sim <- rnorm( N )
Q_sim <- sample( 1:4 , size=N , replace=TRUE )
E_sim <- rnorm( N , U_sim + Q_sim )
W_sim <- rnorm( N , U_sim + bEW*E_sim )
dat_sim <- list( 
  W=standardize(W_sim) , 
  E=standardize(E_sim) , 
  Q=standardize(Q_sim) )


# standard regression of W on E (is biased by U)
m14.4 <- ulam(
  alist(
    W ~ dnorm( mu , sigma ),
    mu <- aW + bEW*E,
    aW ~ dnorm( 0 , 0.2 ),
    bEW ~ dnorm( 0 , 0.5 ),
    sigma ~ dexp( 1 )
  ) , data=dat_sim , chains=4 , cores=4 )
precis( m14.4 )

# instrumental variable regression (unbiased)
# note covariance between W and E measured by Rho[1,2]

m14.5 <- ulam(
  alist(
    c(W,E) ~ multi_normal( c(muW,muE) , Rho , Sigma ),
    muW <- aW + bEW*E,
    muE <- aE + bQE*Q,
    c(aW,aE) ~ normal( 0 , 0.2 ),
    c(bEW,bQE) ~ normal( 0 , 0.5 ),
    Rho ~ lkj_corr( 2 ),
    Sigma ~ exponential( 1 )
  ), data=dat_sim , chains=4 , cores=4 )
precis( m14.5 , depth=3 )



### BRMS implimetnation

#https://rpsychologist.com/adherence-analysis-IV-brms

dat_sim <- as.data.frame(dat_sim)

# W: wages
# Q: quarter of year of birth
# E: education


f1 <- bf(E ~ Q)
f2 <- bf(W ~ E)



# default_prior(f1 + f2, data = dat_sim)
# define some priors
# bprior <- c(prior_string("normal(0,10)", class = "b"))


nlprior <- c(prior(normal(0, 1), nlpar = "E_Intercept"),
             prior(normal(0, 1), nlpar = "W_Intercept"),
             prior(normal(0, 1), nlpar = "E_Q"))


nlprior <- c(set_prior("normal(0,1)", class = "b", coef = ""))


# fit the model
IV_brm <- brm(f1 + f2, data = dat_sim, cores = 4, prior = nlprior)


prior_summary(IV_brm)

summary(IV_brm)



bayestestR::bayesfactor_parameters(IV_brm)



