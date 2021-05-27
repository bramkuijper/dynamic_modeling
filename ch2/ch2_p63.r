library("raster")
library("ggplot2")


#' Represents the algorithm on p53-54 of Mangel & Clark
#' with the modification that a state variable increases
#' at a maximum rate
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
#' @param phi terminal fitness function
#' @export

iterate.sys <- function(params_patches, params_other)
{
  # cutoff below which individual dies
  xc <- params_other$xc

  # max dimension C (see Mangel & Clark p53)
  C <- params_other$C

  # maximum increase in state between timesteps
  max.increase <- params_other$R

  # values of x
  x.vals <- seq(0, C, 1)

  # vector representing F(x, t+1, T)
  # initialize with values x0 = xc to xmax = dim
  F1 <- ifelse(x.vals > xc,1,0)

  Tmax <- params_other$Tmax

  data.all <- NULL

  # loop in a reversed fashion
  # from Tmax to larger timespans
  for (t in seq(Tmax,1,-1))
  {
    # data.frame giving reproductive value
    # and patch choice for each level of x
    F0_D <- data.frame(x=x.vals
                       ,Vmax=0
                       ,patch=NA)

      # now step 2 of the algorithm, cycle over the
      # resource values (x) and then for each value of x
      # cycle over the different patch combinations
      # asking the question: what patch to choose at time t
      # when I have x resources
      for (iter.x.i in 1:length(x.vals))
      {
        x.val.i <- x.vals[[iter.x.i]]

        if (x.val.i <= xc)
        {
          next
        }

        Vi_max <- 0
        patch_max <- 0

        # now check which patch gives maximal gains
        for (patch.i in 1:nrow(params))
        {
          alpha.i  <- params[params$patch == patch.i,"alpha"]
          Y.i <- params[params$patch == patch.i,"Y"]
          lambda.i <- params[params$patch == patch.i,"lambda"]
          beta.i <- params[params$patch == patch.i,"beta"]

          # calculate x'
          xprime <- clamp(x.val.i - alpha.i + Y.i
                         ,xc
                         ,C
                         )

          if (xprime - x.val.i > max.increase)
          {
            xprime <- x.val.i + max.increase
          }
          else
          {
            print("whatevz")
          }

          # calculate x''
          xprimeprime <- clamp(x.val.i - alpha.i
                               ,xc
                               ,C)

          print(paste0("x'': ",xprimeprime
                      ,"x': ",xprime
                      ,"max.increase: ",max.increase))

          # check if the values of F1 exist for this level of x
          stopifnot(!is.na(F1[[match(xprime,x.vals)]]))
          stopifnot(!is.na(F1[[match(xprimeprime,x.vals)]]))

          # calculate Vi
          Vi <- clamp(
                (1.0 - beta.i) * (
                lambda.i * F1[[match(xprime,x.vals)]]
                + (1.0 - lambda.i) * F1[[match(xprimeprime,x.vals)]]
              )
              ,0
              ,1
          )

          if (Vi > Vi_max)
          {
            Vi_max <- Vi
            patch_max <- patch.i
          }
        } # end for patch.i

        F0_D[F0_D$x == x.val.i,c("Vmax","patch")] <- c(Vi_max, patch_max)
      } # end for x.val.i

      # step 4 of the algorithm
      F1 <- F0_D$Vmax

      # add the current value of t to the data.frame
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
