# dynamic modeling exerise 3
# Mangel & Clark experiment 3, chapter 2 page 60
# we need to implement terminal fitness
# function 
# (x - xc) / (x + x0)
# hence even when x > xc, 
# survival is less than one
# and if x0 is large, survival is 
# unlikely to be larger than 0.5
# hence, for relatively low x you now
# see that there is way more risk taking
# in the last time step relative to earlier
# ones (this was previously the other way around).
# Having too little resources in the beginning is 
# unlikely to be important, this only becomes important 
# for terminal survival
# so hence one can take it easy on the risk at the beginning.
# However, close to the end (t close to T)
# the amount of resources becomes v important, as it determines
# terminal survival


# parameters
library("raster")
library("ggplot2")

# let's first find out how x, x0 and xc relate
x.var.graph <- seq(1,10,0.1)
x0.graph <- c(1,2,3,4)
xc.graph <- c(1,2,3,4)

term.fit <- function(row)
{
  x = row["x"]
  x0 = row["x0"]
  xc = row["xc"]
  return((x-xc)/(x+x0))
}

vars.graph <- expand.grid(x=x.var.graph
                          ,x0=x0.graph
                          ,xc=xc.graph)

vars.graph[,"y"] <- apply(X=vars.graph
                          ,MARGIN=1
                          ,FUN=term.fit)

ggplot(data=vars.graph,mapping=aes(x=x, y=y)) +
  geom_line() +
  facet_grid(x0 ~ xc, labeller=label_both) +
  xlab("Individual's resource level, x") +
  ylab("Terminal fitness function, y")

ggsave("page60_term_fit_func.pdf")

params <- data.frame(
  patch=1:3
  ,beta=c(0,0.004,0.02)
  ,alpha=1
  ,lambda=c(0,0.4,0.6)
  ,Y=c(0,3,5)
)

xc <- 3

# max dimension C (see Mangel & Clark p53)
C <- 10

# values of x
x.vals <- seq(0, C, 1)

x0 <- 4.0 

# define terminal fitness function
phi <- function(x) 
{
  if (x <= xc)  
  {
    return(0)
  }
  
  return((x - xc)/(x + x0))
}

# vector representing F(x, t+1, T)
F1 <- sapply(X=x.vals, FUN=phi)




Tmax <- 40

# the total data frame of the whole iteration
# this will gradually be filled
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
  # values of x and then for each value of x
  # cycle over the different patch combinations
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
      
      # calculate x''
      xprimeprime <- clamp(x.val.i - alpha.i 
                           ,xc
                           ,C)
      
      # get the survival values from the next
      # timestep as we are doing backwards iteration
      print(F1)
      
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
      
      print(paste("t:"
                  ,t
                  ,"Vi:"
                  ,Vi
                  ,"patch:"
                  ,patch.i
                  ))
      
      if (Vi > Vi_max)
      {
        Vi_max <- Vi
        patch_max <- patch.i
      }
    } # end for patch.i
    
    
    F0_D[F0_D$x == x.val.i,c("Vmax","patch")] <- c(Vi_max, patch_max)
  } # end for x.val.i

  print(F0_D)  
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

data.all$t <- as.factor(data.all$t)

data.all <- data.all[data.all$t %in% c(1,10,Tmax-3,Tmax),]

ggplot(data=data.all[data.all$x>=xc,], aes(x=x,y=patch)) +
  geom_line(aes(color=t)) +
  xlab("Individual's resource level, x") +
  ylab("Risk vs non-risk patch choice, x") +

ggsave("page60_output_graph.pdf")