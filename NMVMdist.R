
#---------------------------- likelihood function ------------------------------

BesselK = function(x, kappa, log.val = T){
  F1 = log( besselK(x, kappa, expon.scaled = TRUE) ) - x
  F2 = besselK(x, kappa)
  ifelse(log.val == T, return(F1), return(F2))
}

Log.like = function(cc,
                    y,
                    mu,
                    sigma2,
                    tau,
                    nu,
                    Family = c("AL", "GHST", "VG", "H", "Normal", "NMVBS", "NMVL", "NIG"),
                    cens = cens,
                    result = c("PDF", "likelihood", "Log.like")) {
  
  IND = which(cc == 0)
  
  if (Family == "Normal") {
    PDF = dnorm(y, mu, sqrt(sigma2), log = T)
    for (k in IND) {
      if (cens == "Left") PDF[k] = pnorm(y[k], mu[k], sqrt(sigma2), lower.tail = T, log.p = T)
      if (cens == "Right") PDF[k] = pnorm(y[k], mu[k], sqrt(sigma2), lower.tail = F, log.p = T)
    }
  }
  
  if (Family == "GHST") {
    dGHST = function(y, mu, sigma2, tau, nu, log = F) {
      z =  (y - mu) / sqrt(sigma2)
      out1 = (nu / 2) * log(nu) + (1 - nu / 2) * log(2) - lgamma(nu / 2) - 0.5 * log(2 * pi * sigma2) +
        log(abs(tau) / sqrt(sigma2)) * (nu + 1) / 2
      out2 = z * tau / sqrt(sigma2) - ((nu + 1) / 4) * log(nu + z^2)
      out3 = BesselK(sqrt((nu + z^2) * tau^2 / sigma2), (nu + 1) / 2, log.val = T)
      
      out = out1 + out2 + out3
      ifelse(log == T, return(out), return(exp(out)))
    }
    
    pGHST = function(y, mu, sigma2, tau, nu, lower.tail = T) {
      z = y - mu
      dinvgamma = function(x, shape, rate) {
        out = shape * log(rate) - lgamma(shape) - (shape + 1) * log(x) - (rate/x)
        return(out)
      }
      f1 = function(w) {
        exp( dinvgamma(w, nu / 2, nu / 2) + pnorm(z, mean = tau * w, sd = sqrt(w * sigma2),
                                                  log.p = T, lower.tail = lower.tail))
      }
      resp = integrate(f1, 0, Inf, rel.tol = .Machine$double.eps^0.5,
                       stop.on.error = FALSE)$value
      return(resp)
    }
    
    PDF = dGHST(y, mu, sigma2, tau, nu, log = T)
    for (k in IND) {
      if (cens == "Left") PDF[k] = log(pGHST(y[k], mu[k], sigma2, tau, nu, lower.tail = T))
      if (cens == "Right") PDF[k] = log(pGHST(y[k], mu[k], sigma2, tau, nu, lower.tail = F))
    }
  }
  
  if (Family == "VG") {
    dVG = function(y, mu, sigma2, tau, nu, log = F) {
      z =  (y - mu) / sqrt(sigma2)
      out1 = nu * log(nu) + (1/2 - nu) * log(2) - lgamma(nu ) - 
        0.5 * log( pi * sigma2) + (nu - 1 / 2) * log(abs(z))
      out2 = z * tau / sqrt(sigma2) - (nu/2 - 1/4) * log(nu + tau^2/sigma2)
      out3 = BesselK(sqrt((nu + tau^2 / sigma2)*z^2),nu - 1 / 2, log.val = T)
      
      out = out1 + out2 + out3
      ifelse(log == T, return(out), return(exp(out)))
    }
    
    pVG = function(y, mu, sigma2, tau, nu, lower.tail = T) {
      z = y - mu
      dgamma = function(x, shape, rate) {
        out = shape * log(rate) - lgamma(shape) + (shape - 1) * log(x) - (rate*x)
        return(out)
      }
      f1 = function(w) {
        exp( dgamma(w, nu, nu / 2) + pnorm(z, mean = tau * w, sd = sqrt(w * sigma2),
                                           log.p = T, lower.tail = lower.tail))
      }
      resp = integrate(f1, 0, Inf, rel.tol = .Machine$double.eps^0.5,
                       stop.on.error = FALSE)$value
      return(resp)
    }
    
    PDF = dVG(y, mu, sigma2, tau, nu, log = T)
    for (k in IND) {
      if (cens == "Left") PDF[k] = log(pVG(y[k], mu[k], sigma2, tau, nu, lower.tail = T))
      if (cens == "Right") PDF[k] = log(pVG(y[k], mu[k], sigma2, tau, nu, lower.tail = F))
    }
  }
  
  if (Family == "H") {
    dH = function(y, mu, sigma2, tau, nu, log = F) {
      z =  (y - mu) / sqrt(sigma2)
      out1 = -0.5 * log(nu*sigma2 + tau ^ 2) - log(2) - BesselK(nu, 1, log.val = T)
      out2 = z * tau / sqrt(sigma2) - sqrt((tau^2/sigma2 + nu)*(z ^ 2 + nu))
      out = out1 + out2 
      ifelse(log == T, return(out), return(exp(out)))
    }
    
    pH = function(y, mu, sigma2, tau, nu, lower.tail = T) {
      z = y - mu
      dgig = function(x,kappa,psi,chi) {
        out = (1/2) * kappa * log(psi/chi) + (kappa - 1) * log(x) - 
          1 / 2 * (chi / x + psi * x) - log(2) - BesselK(sqrt(psi * chi), kappa, log.val = T)
        return(out)
      }
      f1 = function(w) {
        exp( dgig(w, kappa = 1, psi = nu, chi = nu) + pnorm(z, mean = tau * w, sd = sqrt(w * sigma2),
                                                            log.p = T, lower.tail = lower.tail))
      }
      resp = integrate(f1, 0, Inf, rel.tol = .Machine$double.eps^0.5,
                       stop.on.error = FALSE)$value
      return(resp)
    }
    
    PDF = dH(y, mu, sigma2, tau, nu, log = T)
    for (k in IND) {
      if (cens == "Left") PDF[k] = log(pH(y[k], mu[k], sigma2, tau, nu, lower.tail = T))
      if (cens == "Right") PDF[k] = log(pH(y[k], mu[k], sigma2, tau, nu, lower.tail = F))
    }
  }
  
  if (Family == "AL") {
    dAL = function(y, mu, sigma2, tau, log = F){
      kappa = sqrt(tau^2 / sigma2 + 1)
      z = y - mu
      out = (tau * z / sigma2 - kappa * abs(z) / sqrt(sigma2)) - log(sqrt(sigma2) * kappa) - log(2)
      ifelse(log == T, return(out), return(exp(out)))
    }
    
    pAL = function(y, mu, sigma2, tau, lower.tail = T) {
      z = y - mu
      dexp = function(x) {
        out = -log(2) - 0.5 * x
        return(out)
      }
      
      f1 = function(w) {
        exp( dexp(w) + pnorm(z, mean = tau * w, sd = sqrt(w * sigma2),
                             log.p = T, lower.tail = lower.tail))
      }
      resp = integrate(f1, 0, Inf, rel.tol = .Machine$double.eps^0.5,
                       stop.on.error = FALSE)$value
      return(resp)
    }
    
    PDF = dAL(y, mu, sigma2, tau, log = T)
    for (k in IND) {
      if (cens == "Left")  PDF[k] = log(pAL(y[k], mu[k], sigma2, tau, lower.tail = T))
      if (cens == "Right") PDF[k] = log(pAL(y[k], mu[k], sigma2, tau, lower.tail = F))
    }
  }
  
  if (Family == "NIG") {
    
    fGH = function(x, mu, sigma2, tau, kappa, chi, psi, log.val = F){
      chii = chi + (x-mu)^2 / sigma2
      psii = psi + tau^2 / sigma2
      out1 = BesselK(sqrt(chii * psii), kappa - 0.5, log.val = T) -
        BesselK(sqrt(chi * psi), kappa, log.val = T)
      out2 = (kappa/2 - 1/4) * log(chii /  psii) + (kappa / 2) * log(psi/chi)
      out3 = (x - mu) * tau / sigma2 - 0.5 * log(sigma2 * 2 * pi)
      F1 = out1 + out2 + out3
      ifelse(log.val == T, return(F1), return(exp(F1)))
    }
    
    dNIG = function(x, mu, sigma2, lambda, alpha, log = F) {
      out = fGH(x, mu, sigma2, lambda, -0.5, 1, alpha^(2), log.val = F)
      ifelse(log == F, return(out), return(log(out)))
    }
    
    pNIG= function(y, mu, sigma2, tau, nu, lower.tail = T) {
      z = y - mu
      dgig = function(x,kappa,psi,chi) {
        out = (1/2) * kappa * log(psi/chi) + (kappa - 1) * log(x) - 
          1 / 2 * (chi / x + psi * x) - log(2) - BesselK(sqrt(psi * chi), kappa, log.val = T)
        return(out)
      }
      f1 = function(w) {
        exp( dgig(w, kappa = -0.5, psi = nu^2, chi = 1) + pnorm(z, mean = tau * w, sd = sqrt(w * sigma2),
                                                                log.p = T, lower.tail = lower.tail))
      }
      resp = integrate(f1, 0, Inf, rel.tol = .Machine$double.eps^0.5,
                       stop.on.error = FALSE)$value
      return(resp)
    }
    
    PDF = dNIG(y, mu, sigma2, tau, nu, log = T)
    for (k in IND) {
      if (cens == "Left") PDF[k] = log(pNIG(y[k], mu[k], sigma2, tau, nu, lower.tail = T))
      if (cens == "Right") PDF[k] = log(pNIG(y[k], mu[k], sigma2, tau, nu, lower.tail = F))
    }
  }
  
  if (Family == "NMVBS") {
    
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
    
    dNMVBS = function(x, mu, sigma2, lambda, alpha, log = F) {
      out = 0.5 * (fGH(x, mu, sigma2, lambda, 0.5, alpha^(-2), alpha^(-2), log.val = F)+
                     fGH(x, mu, sigma2, lambda, -0.5, alpha^(-2), alpha^(-2), log.val = F))
      ifelse(log == F, return(out), return(log(out)))
    }
    
    pNMVBS = function(y, mu, sigma2, tau, nu, lower.tail = T) {
      z = y - mu
      
      a.Fun = function(x, A, B){
        (1/A) * (sqrt(x / B) - sqrt(B / x))
      }
      A.Fun = function(x, A, B){
        (x + B)/(2 * A * sqrt(B) * sqrt(x ^ 3))  
      }
      dBS = function(x, alpha, beta) {
        out = log(A.Fun(x, alpha, beta)) + dnorm(a.Fun(x, alpha, beta), log = T)
        return(out)
      }
      
      f1 = function(w) {
        exp(dBS(w, nu, 1) + pnorm(z, mean = tau * w, sd = sqrt(w * sigma2),
                                  log.p = T, lower.tail = lower.tail))
      }
      resp = integrate(f1, 0, Inf, rel.tol = .Machine$double.eps^0.5,
                       stop.on.error = FALSE)$value
      return(resp)
    }
    
    PDF = dNMVBS(y, mu, sigma2, tau, nu, log = T)
    for (k in IND) {
      if (cens == "Left") PDF[k] = log(pNMVBS(y[k], mu[k], sigma2, tau, nu, lower.tail = T))
      if (cens == "Right") PDF[k] = log(pNMVBS(y[k], mu[k], sigma2, tau, nu, lower.tail = F))
    }
  }
  if (Family == "NMVL") {
    dNMVL = function(x, mu, sigma2, lambda, alpha, log = F) {
      z = (y - mu) / sqrt(sigma2)
      chi = z ^ 2 
      psi  = tau ^ 2 / sigma2 + 2 * nu
      out1 = 2 * nu ^ 2 * exp(tau * z /sqrt(sigma2)) *(chi /psi)^(1/4)/((1+nu)*sqrt(2*pi*sigma2))
      out2 = BesselK(sqrt(chi* psi),  0.5, log.val = F) + sqrt(chi / psi)*
        BesselK(sqrt(chi * psi), 1.5, log.val = F)
      
      out = out1 * out2    
      ifelse(log == F, return(out), return(log(out)))
    }
    
    
    pNMVL= function(y, mu, sigma2, tau, nu, lower.tail = T) {
      z = y - mu
      
      dL = function(x, alpha) {
        out = 2 * log(alpha) + log(1 + x) - alpha * x -log (1 + alpha)
        return(out)
      }
      
      f1 = function(w) {
        exp(dL(w, nu) + pnorm(z, mean = tau * w, sd = sqrt(w * sigma2),
                              log.p = T, lower.tail = lower.tail))
      }
      resp = integrate(f1, 0, Inf, rel.tol = .Machine$double.eps^0.5,
                       stop.on.error = FALSE)$value
      return(resp)
    }
    
    PDF = dNMVL(y, mu, sigma2, tau, nu, log = T)
    for (k in IND) {
      if (cens == "Left") PDF[k] = log(pNMVL(y[k], mu[k], sigma2, tau, nu, lower.tail = T))
      if (cens == "Right") PDF[k] = log(pNMVL(y[k], mu[k], sigma2, tau, nu, lower.tail = F))
    }
  }
  
  if (result == "PDF") Out = exp(PDF)
  if (result == "Log.like") Out = sum(PDF)
  if (result == "likelihood") Out = prod(exp(PDF))
  return(Out)
}


