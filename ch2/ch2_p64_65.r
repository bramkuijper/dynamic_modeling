# dynamic modeling exercise
# Mangel & Clark chapter 'Variable Handling Times'

library("raster")
library("ggplot2")


#' Represents the algorithm on p53-54 of Mangel & Clark
#' with a modification of variable handling times, as described
#' on p64-65
#' @param params_patches list with patch-specific parameters
#' as described on page 45, 54:
#'    beta (mortality prob)
#'    alpha (cost of just being alive)
#'    lambda (food encounter prob)
#'    Y (increment in state variable)
#' @param params_other list with parameters regarding the
#' algorithm:
#'    xc (minimum resource level before death)
#'    C maximum possible level of resources
#'    R maximum increase in state between timesteps
#'      (see p64, experiment 6)
#'    Tmax maximum time steps per iteration
#'    max_tau maximum tau value
#' @param phi terminal fitness function
#' @export
iterate.sys <- function(params_patches
                        ,params_other
                        ,phi=NA)
{
  # minimum resource level
  # below which individual dies
  xc <- params_other$xc

  # max dimension C (see Mangel & Clark p53)
  C <- params_other$C

  # max increase in resource between timesteps
  max.increase <- params_other$R

  # maximum handling time
  max.tau <- params_other$max_tau

  # values of x
  x.vals <- seq(0, C, 1)

  # terminal survival function
  if (is.na(phi))
  {
      phi <- function(x) 
      {
        return(ifelse(x > xc,1,0))
      }
  }

  # maximum iteration time (note that
  # we perform backwards calculations)
  Tmax <- params_other$Tmax

  # vector representing F(x, t+1, T)
  # initialize with values x0 = xc to xmax = C
  F1.x <- sapply(X=x.vals, FUN=phi)

  # make sequence of tau timesteps
  tau.vals <- seq(from=1,to=max.tau,by=1)
  F1.tau <- tau.vals

  # make the F1 vector in x and tau, in which
  # the column F1 contains
  # the actual values of F(x',t+1,T)
  F1.df <- expand.grid(x=F1.x, tau=F1.tau)
  F1.df$F1 <- NA
  
  # initially, one sets F1(x,1) = phi(x)
  # and F1(x,tau) = 0 for tau > 1
  F1.df[F1.df$tau == 1,"F1"] <- phi(F1.df[F1.df$tau == 1,"x"])

  # auxiliary variable to store all the data
  data.all <- NULL

  # loop in a reversed fashion
  # from t=Tmax to t=1
  for (t in seq(Tmax,1,-1))
  {
    # data.frame giving reproductive value
    # and patch choice for each level of x and tau
    
    # first allocate x column to 0
    F0_D <- data.frame(x=rep(0,times=length(F1.x)))
    F0_D$Fx_t <- 0
    F0_D$patch <- NA
    F0_D$tau <- NA

    # now step 2 of the algorithm, cycle over the
    # resource values (x) and then for each value of x
    # cycle over the different patch combinations
    # asking the question:
    # - what patch best to choose at time t
    # when I have x resources
    # - what handling time best to choose in that patch
    # when I have x resources
    for (iter.x.i in 1:length(x.vals))
    {
      x.val.i <- x.vals[[iter.x.i]]

      # F(x,t,T) = 0 for x<=x_{c}
      if (x.val.i <= xc)
      {
        next
      }

      # reset auxiliary variables
      # which keep track of which patch and handling
      # time has the maximum payoff
      Fx_t_max <- 0
      patch_max <- 0
      handling_time_max <- 0

      # loop through all the values of tau (handling time)
      # tau.i == 0 corresponds to t'=t in F(x,t',T) on p65, 2nd para
      for (tau.i in seq(from=0,to=tau_max,by=1))
      {
        # now check which patch gives maximal gains
        for (patch.i in 1:nrow(params))
        {
          alpha.i  <- params[params$patch == patch.i,"alpha"]
          Y.i <- params[params$patch == patch.i,"Y"]
          lambda.i <- params[params$patch == patch.i,"lambda"]
          beta.i <- params[params$patch == patch.i,"beta"]

          # calculate x' according to eq. (2.13) on page 64
          xprime <- clamp(x=x.val.i - tau.i * alpha.i + Y.i
                         ,lower=xc
                         ,upper=C
                         )
          # see whether max increment has been exceeded
          if (xprime - x.val.i > max.increase)
          {
            xprime <- x.val.i + max.increase
          }

          # calculate x'' according to eq. (2.13) on page 64
          xprimeprime <- clamp(x=x.val.i - alpha.i
                               ,lower=xc
                               ,upper=C)

          # see whether max increment has been exceeded
          if (xprimeprime - x.val.i > max.increase)
          {
            xprimeprime <- x.val.i + max.increase
          }

          # check if the values of F1 exist for this level of x
          stopifnot(!is.na(F1[[match(xprime,x.vals)]]))
          stopifnot(!is.na(F1[[match(xprimeprime,x.vals)]]))

          # calculate Vi  
          Fx_t.unclamped <-
                (1.0 - beta.i)^(tau.i) * (
                lambda.i * F1.df[F1.df$tau == tau.i & F1.df$x == xprime,"F1"]
                + (1.0 - lambda.i) * F1.df[F1.df$tau == 1 & F1.df$x == xprime,"F1"]
                )

          Fx_t <- clamp(x=Fx_t.unclamped, lower=0, upper=1)

          if (Fx_t > Fx_t_max)
          {
            Fx_t_max <- Fx_t
            patch_max <- patch.i
            handling_time_max <- tau.i
          }
        } # end for patch.i
      } # end for tau.i

      F0_D[F0_D$x == x.val.i,c("Fx_t","patch","tau")] <-
        c(Fx_t_max, patch_max, handling_time_max)
    } # end for x.val.i

    # step 4 of the algorithm
    # F1(x,tau) = F1(x,tau-1)
    # F1(x,tau=1) = F0(x)
    for (tau.i in seq(from=max.tau, to=2, by=-1))
    {
      # oh shit we need to do something about the
      # specific values of x...
      F1.df[F1.df$tau == tau.i,"F1"] <-
        F1.df[F1.df$tau == tau.i - 1,"F1"]
    }

    F1.df[F1.df$tau == 1] <- F0_D$Fx_t_max



    F1 <- F0_D$Fx_t

      # make another column with the value of t
      # and add that to the F0 column
      F0_D$t <- t

      # append data
      if (is.null(data.all)) {
        data.all <- F0_D
      } else {
        data.all <- rbind(data.all,F0_D)
      }
  } # end for (t in seq(T,0,-1))

  return(data.all)
} # end function iterate sys.

# set all the parameters for this example
params <- data.frame(
  patch=1:3
  ,beta=c(0,0.004,0.02)
  ,alpha=1
  ,lambda=c(0,0.4,0.6)
  ,Y=c(0,3,5)
)

params.run <- list(
  xc=3
  ,C=100
  ,R=1
  ,Tmax=40)

data.iter <- iterate.sys(
  params_patches=params
  ,params_other=params.run)

data.iter$C <- 100

params.run["C"] <- 200

data.iter2 <- iterate.sys(
  params_patches=params
  ,params_other=params.run)

data.iter2$C <- params.run$C

data.all <- rbind(data.iter,data.iter2)
data.all$t <- as.factor(data.all$t)

data.all <- data.all[data.all$t %in% c(1,10,params.run$Tmax),]

ggplot(data=data.all[data.all$x>=params.run$xc,], aes(x=x,y=patch)) +
  geom_line(aes(color=t)) +
  facet_grid(. ~ C) +
  xlab("Individual's resource level, x") +
  ylab("Risk vs non-risk patch choice, x") +

ggsave("page64_experiment6_output_graph.pdf")
