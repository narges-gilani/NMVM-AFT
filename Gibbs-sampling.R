#-------------------------Gibbs.sampleing---------------------------------------

Gibbs.sampleing = function(cc,
                           y,
                           x,
                           n.iter,
                           n.thin,
                           burnin,
                           cens = c("Left", "Right"),
                           hyper_set = hyper_set,
                           prior = c("Exp", "Hierar_1", "Unif"),
                           hyper = hyper,
                           Family = c("AL", "GHST", "VG", "H", "Normal", "NMVBS", "NIG", "NMVL"),
                           n.chains = n.chains,
                           Print = T,
                           Prog.Bar = F) {
  n = length(y)
  Repli = 0
  x = as.matrix(x)
  p = ncol(x)
  n = nrow(x)
  
  if (cens == "Left") {
    lower = rep(-Inf, n); upper = y
  }
  
  if (cens == "Right") {
    lower = y; upper = rep(Inf, n)
  }
  
  #-------------------------------------------------------------------------------
  #-------------------------- (Required Functions) -------------------------------
  fGH = function(x, mu, sigma2, tau, kappa, chi, psi, log.val = F){
    chii = chi + (x-mu)^2 / sigma2
    psii = psi + tau^2 / sigma2
    out1 = BesselK(sqrt(chii * psii), kappa - 0.5, log.val = T) -
      BesselK(sqrt(chi * psi), kappa, log.val = T)
    out2 = (kappa/2 - 1/4) * log(chii /  psii)
    out3 = (x - mu) * tau / sigma2 - 0.5 * log(sigma2 * 2 * pi)
    F1 = out1 + out2 + out3
    ifelse(log.val == T, return(F1), return(exp(F1)))
  }
  
  prob <- function(y, mu, sigma2, tau, nu){
    dout = fGH(y, mu, sigma2, tau, -0.5, nu^(-2), nu^(-2), log.val = T) - 
      fGH(y, mu, sigma2, tau, 0.5, nu^(-2), nu^(-2), log.val = T)
    1 / (1 + exp(dout))
  }
  
  prob2 <- function(y,mu, sigma2, tau, nu){
    psi = tau ^2 / sigma2 + 2  * nu 
    chi = (y - mu)^2 / sigma2
    dout = 0.5 * log(chi/psi) + BesselK(sqrt(psi*chi), 1.5, log.val = T) - BesselK(sqrt(psi*chi), 0.5, log.val = T)
    1 / (1 + exp(dout))
  }
  
#-------------------------------------------------------------------------------
  
  r.normal.multiv = function(n = 1,
                             mu = NA,
                             Sigma = NA,
                             method = c("eigen", "chol"))
  {
    massage = cli::combine_ansi_styles("red", "bold")
    
    if (!isSymmetric(Sigma,
                     tol = sqrt(.Machine$double.eps),
                     check.attributes = FALSE)) {
      stop(massage("Sigma is not a symmetric matrix"))
    }
    if (length(mu) != nrow(Sigma))
      stop(massage("different size of mu and Sigma"))
    method <- match.arg(method)
    
    R = if (method == "eigen") {
      ev <- eigen(Sigma, symmetric = TRUE)
      if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
        warning(massage("Sigma is numerically not positive semidefinite"))
      }
      t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 0))))
    }
    else if (method == "chol") {
      R = chol(Sigma)
    }
    out = matrix(rnorm(n * ncol(Sigma)), nrow = n, byrow = TRUE) %*% R
    out <- sweep(out, 2, mu, "+")
    out
  }
  
  truncgamma <- function(n, lower, upper, rate, scale) {
    u = runif(n, min = pgamma(lower, rate, scale), 
              max = pgamma(upper, rate, scale))
    qgamma(u, rate, scale)
  }
  
#-------------------------------------------------------------------------------
#-------------------------- (Hyper-parameter setting)---------------------------
#-------------------------------------------------------------------------------
  if (!is.na(hyper_set)) {
    mu0 = rep(hyper_set[1], p)
    Sigma0 = hyper_set[2] * diag(p)
    alfa1 = hyper_set[3]
    alfa2 = hyper_set[4]
    a0 = hyper_set[5]
    b0 = hyper_set[6]
  }
  
  if (is.na(hyper_set)) {
    mu0 = matrix(0, p, 1)
    Sigma0 = 1000 * diag(p)
    alfa1 = 0.01
    alfa2 = 0.01
    a0 = 0
    b0 = 1
  }
  
  c_2 = 0.01 
  d_2= ifelse (Family == "GHST", 0.5, 1000)
  
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
  
  sigma2 = tau = nu = lambda = Aux = c()
  mu = W = pr = matrix(0, n.iter * n.chains, n)  
  beta = matrix(0, n.iter * n.chains, p)
  
  if (Family == "Normal") {
    W = matrix(1, n.iter * n.chains, n)
    tau = nu = rep(0, n.iter * n.chains)
  }
  
  if (Family == "AL")   nu = rep(0, n.iter * n.chains)
  
  set1 = set = c()
  
  while (Repli < n.chains) {
    Repli = Repli + 1
    Cont = ceiling(n.iter / 10)
    n.iter = Cont * 10
    Lin = n.iter * (Repli - 1)
    G = seq(Cont + Lin, n.iter * Repli, Cont)
    
    cat(cli::rule(line = 2, line_col = "orange"), "\n")
    massage = cli::combine_ansi_styles("orange", "bold")
    nbld = "Initialization"
    cli::cli_alert_info(massage("{.emph {nbld}} is on process. Please wait..."))
    
    reg = lm(y ~ x[, 2:p])
    
    if (Family != "Normal")  tau[1 + Lin] = 0.5
    
    beta[1 + Lin,] = as.vector(coefficients(reg), mode = "numeric")
    sigma2[1 + Lin] = sum((y - x %*% (beta[1 + Lin,])) ^ 2) / (n - p)
    
    if (Family != "Normal" && Family != "AL") {
      nu[1 + Lin] = 5
      lambda[1 + Lin] = hyper  
    }
    
    mu[1 + Lin,] = x %*% (beta[1 + Lin,])
    
    if (Family == "GHST") 
      W[1 + Lin, ] = 1 / rgamma(n, shape = nu[1 + Lin]/2, rate = nu[1 + Lin]/2)
    
    if (Family == "AL") 
      W[1 + Lin, ] <- rgamma(n, 1, 1/2)
    
    if (Family == "H") 
      W[1 + Lin, ] = ghyp :: rgig(n, 1,  nu[1 + Lin], nu[1 + Lin])
    
    if (Family == "NIG") 
      W[1 + Lin, ] = ghyp :: rgig(n, -0.5, 1, nu[1 + Lin]^2)
    
    if (Family == "VG") 
      W[1 + Lin, ] = rgamma(n, shape = nu[1 + Lin], rate = nu[1 + Lin]/2)
    
    if (Family == "NMVBS") {
      Zgen = rnorm(n)
      W[1 + Lin, ] = 0.25 * (nu[1 + Lin] * Zgen + sqrt(( nu[1 + Lin] * Zgen) ^ 2 + 4 ) ) ^ 2 
    }
    
    if (Family == "NMVL") {
      W[1 + Lin, ] = rlindley(n, theta=nu[1 + Lin], mixture=TRUE)
    }
    
    massages = cli::combine_ansi_styles("green", "bold")
    cli::cli_alert_success(massages("{.emph {nbld}} is finished."))
    
    if (Print) {
      massage = cli::combine_ansi_styles("orange", "bold")
      cli::cli_alert_info(massage("Gibbs sampling"))
      cli::cli_progress_bar(
        total = n.iter,
        format = "{cli::pb_bar} {cli::pb_percent} [{cli::pb_elapsed}]"
      )
    }
    
    if (Prog.Bar) {
      Start <-
        tcltk::tkProgressBar(
          title = paste("Progress bar of Bayesian inference of the ", Family, cens, "censored linear regression"),
          min = 0,
          max = n.iter,
          width = 1000
        )
    }
    
    cat("\n")
    Fil = 2 + Lin
    Col = n.iter * Repli
    for (j in Fil:Col) {
      if (Prog.Bar) {
        tcltk::setTkProgressBar(Start,
                                j / Repli,
                                label = paste(
                                  round(j / (Repli * n.iter) * 100, 0),
                                  "% of the maximum iteration is done! Please wait..."
                                ))
      }
      
      if (Print) {
        Sys.sleep(2/n.iter)
        cli::cli_progress_update()
      }
      
      y[cc == 0] = truncnorm :: rtruncnorm(1,
                                           lower[cc == 0],
                                           upper[cc == 0],
                                           mean = mu[j - 1, cc == 0] + tau[j - 1] * W[j - 1, cc == 0],
                                           sd = sqrt((W[j - 1, cc == 0]) * sigma2[j - 1]))
      
      rate = sum(W[j - 1, ] ^ (-1) * (y - mu[j - 1,] - tau[j - 1] * W[j - 1, ]) ^ 2)
      sigma2[j] = 1 / rgamma(1, shape = alfa1 + n / 2, rate = rate / 2 + alfa2)
      
      xast = sqrt(W[j - 1, ] ^ (-1)) * x
      yast = sqrt(W[j - 1, ] ^ (-1)) * (y - tau[j - 1] * W[j - 1,])
      SigmaA = solve(solve(Sigma0) + t(xast) %*% (xast) / sigma2[j])
      muA = SigmaA %*% (solve(Sigma0) %*% mu0 + t(xast) %*% yast / sigma2[j])
      beta[j,] = r.normal.multiv(1, muA, SigmaA, method = "eigen")
      
      mu[j,] = x %*% (beta[j, ])
      
      if (Family != "Normal") {
        sigma2.tau = (1 / b0 + sum(W[j - 1,]/ sigma2[j])) ^ (-1)
        mu.tau = sigma2.tau * (sum((y - mu[j, ]) / sigma2[j]) + a0 / b0)
        tau[j] = rnorm(1, mean = mu.tau, sd = sqrt(sigma2.tau))
      }
      
      if (Family == "AL") {
        chi = c((y - mu[j, ]) ^ 2 / sigma2[j])
        psi = tau[j]^2 / sigma2[j] + 1
        W[j, ] =  apply(cbind(chi, psi), 1, function(x) ghyp :: rgig(1, 0.5, x[1], x[2]))
      }
      
      if (Family != "Normal" && Family != "AL") {
        if (prior == "Hierar_1") {
          lambda[j] <- truncgamma(1, lower = c_2,  upper = d_2, rate = 2, scale = nu[j - 1])
        }else{
          lambda[j] = hyper
        }
        
        
        if (Family != "Normal" && Family != "AL") 
          nu[j] = MHnu(nu[j - 1], W[j - 1, ], lambda = lambda[j],Family=Family, prior = prior, hyper=hyper)
      }
      
      if (Family == "GHST") {
        chi = c((y - mu[j, ]) ^ 2 / sigma2[j] + nu[j])
        psi = tau[j]^2 / sigma2[j]
        W[j, ] = apply(cbind(chi, psi), 1, function(x) ghyp :: rgig(1, -(nu[j] + 1) / 2, x[1], x[2])) 
      }
      
      if (Family == "H") {
        chi = c((y - mu[j, ]) ^ 2 / sigma2[j] + nu[j])
        psi = tau[j]^2 / sigma2[j] + nu[j]
        W[j, ] =  apply(cbind(chi, psi), 1, function(x) ghyp :: rgig(1, 0.5, x[1], x[2])) 
      } 
      
      if (Family == "NIG") {
        chi = c((y - mu[j, ]) ^ 2 / sigma2[j] + nu[j])
        psi = tau[j]^2 / sigma2[j] + nu[j]
        W[j, ] =  apply(cbind(chi, psi), 1, function(x) ghyp :: rgig(1, -1, x[1], x[2])) 
      }
      
      if (Family == "VG") {
        chi = c((y - mu[j, ]) ^ 2 / sigma2[j])
        psi = tau[j]^2 / sigma2[j] + nu[j]
        W[j, ] =  apply(cbind(chi, psi), 1, function(x) ghyp :: rgig(1, (2 * nu[j] - 1) / 2, x[1], x[2]))       
      }
      
      if (Family == "NMVBS") {
        chi <- c((y - mu[j,]) ^ 2 / sigma2[j] + nu[j] ^ (-2))
        psi <- tau[j] ^ 2 / sigma2[j] + nu[j] ^ (-2)
        pr<- c(prob(y, mu[j,], sigma2[j], tau[j], nu[j]))
        KAPPa = c(0, -1)
        for(i in 1:n){
          kk = rbinom(1, 1, 1 - pr[i]) + 1
          W[j, i] <- ghyp :: rgig(1, KAPPa[kk], chi[i], psi)
        }
      }
      
      if (Family == "NMVL") {
        chi <- c((y - mu[j,]) ^ 2 / sigma2[j])
        psi <- tau[j] ^ 2 / sigma2[j] + 2 * nu[j]
        pr <- c(prob2(y, mu[j,], sigma2[j], tau[j], nu[j]))
        KAPPa = c(0.5, 1.5)
        for(i in 1:n){
          kk = rbinom(1, 1, 1 - pr[i]) + 1
          W[j, i] <- ghyp :: rgig(1, KAPPa[kk], chi[i], psi)
        }
      }
      
    }
    
    if (Print) {
      cli::cli_alert_success(massages("Sampling done!"))
    }
    
    if (Prog.Bar) close(Start)
    cat("\r")
    initial = burnin
    set1 = seq(initial + n.thin + n.iter * (Repli - 1), n.iter *
                 Repli, n.thin)
    set = c(set, set1)
  }
  
  Out = list( beta = beta[set, ], sigma2 = sigma2[set], tau = tau[set],
              W = W[set, ], nu = nu[set], mu = mu[set, ], lambda = lambda[set])    
  
  return( Out )
}

rlindley <- function(n, theta, mixture = TRUE)
{
  stopifnot(theta > 0)
  if(mixture)
  {
    p <- rbinom(n, size = 1, prob = theta / (1 + theta))
    p * rgamma(n, shape = 1, rate = theta) + (1 - p) * rgamma(n, shape = 2, rate = theta)
  }
  else
  {
    qlindley(p = runif(n), theta, lower.tail = TRUE, log.p = FALSE)
  }
}

