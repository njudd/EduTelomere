# Nicholas Judd
# 2024-10-19
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



