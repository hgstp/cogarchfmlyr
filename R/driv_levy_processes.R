#' Simulation of compound Poisson process
#' 
#' comp_pois_proc simulates the trajectory of a compound Poisson process
#' 
#' The jumps of the compound Poisson process arrive accordingly to a poisson process.
#' The sizes of jumps are random, with a specified probability distribution. 
#' 
#' @param start start point of the time interval [start,end]
#' @param end  end point of the time interval [start, end]
#' @param lambda  intensity of compound Poisson process
#' @param distribution expression defining the distribution of the jump sizes
#'   in the form r'dist'
#' @param dparams further parameters of the jump size distribution  
#' @param obs_times observation times of the process in the intervall [start, end]
#' 
#' @return output list containing the jump times, the jump sizes, the observation times
#'   and the compound Poisson process
#'  
#' @examples
#' comp_pois_proc(distribution = rnorm, dparams = list(mean = 0, sd = 2))
#' comp_pois_proc(end = 100, obs_times = seq(1/10, to = 100, by = 1/10),
#'                  distribution = rt, dparams = list(df = 2))    

comp_pois_proc <- function(start = 0, end = 10, lambda = 1, 
                            distribution, dparams = list(), obs_times = start:end){

#########################  
# see Algorithm 6.2 in Cont & Tankov  

  length_interval = end - start

# number of jumps in the interval [start,end]
  N <- rpois(1, lambda = lambda * length_interval)

# jump times in the interval [start, end]
  U <- runif(N, min = start, max = end)

# jump sizes of the N jumps
  Y <- do.call(what = distribution, args = c(list(n = N), dparams))

# trajectory

  support_fun <- function(obt, jt, jumps){
    vapply(jt, FUN = function(w){w < obt}, FUN.VALUE = seq(along = obt))
  }
# support_fun checks which jump times are smaller than the respective
# oberservation times -> matrix with as many rows as observation times 
# and as many columns as there are jump times
#
# we multiple the transposed of support_fun with the corresponding jump sizes
# and calculate the columnwise sums

  X <- colSums(Y * t(support_fun(obs_times, U)))
  

# plot trajetory with ggplot2

# ggplot(data.frame(time = obs_times, X = X), aes(x = time, y = X)) + geom_step()

# output  with the jump times, jump sizes, the observation times and the value
# of the process at the observation times

list(jump_times = sort(U), jump_sizes = Y[order(U)], obs_times = obs_times, 
               process = X)  

}

########################################################################################


#' Simulation of Variance Gamma process
#' 
#' var_gam_proc simulates the trajectory of a variance gamma process 
#' with parameters sigma>0, nu>0 and theta
#' 
#' The process is defined through a time change TG, which is done wrt a gamma process.
#' W_TG is then a time changed standard Brownian Motion and the variance gamma process is
#' defined through
#'      V = theta*TG + sigma * W_TG 
#' 
#' @param sigma standard deviation of the Brownian motion
#' @param nu scale parameter of the Gamma distribution defining the increments
#'   of the Gamma process
#' @param theta  drift parameter
#' @param obs_times observation times of the process
#' 
#' @return output list containing the jump times, the jump sizes, the observation times
#'   and the Variance Gamma process
#'  
#'
#

  var_gam_proc <- function(sigma,nu,theta,obs_times){
    
    # time grid
    diff_obs_times <- diff(obs_times)
    
    # the increments of the gamma process TG
    dTG <- rgamma(length(diff_obs_times), shape = (1/nu) * diff_obs_times, scale = nu)
    
    # the increments of the time changed brownian motion with drift given by theta*dTG
    dV <- rnorm(length(diff_obs_times), mean = theta * dTG, sd = sigma * sqrt(dTG))
    
    list(obs_times = obs_times, V = c(0, cumsum(dV)))  
    
  }
