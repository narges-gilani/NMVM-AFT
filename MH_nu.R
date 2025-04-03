
#--------------------------------MHnu-------------------------------------------

MHnu = function(last, W, lambda, Family=Family, prior =prior, hyper=hyper) {
  n = length(W)
  f1=function(nu, Family=Family){
    if(Family== "GHST") 
      resp= 0.5 * n * nu * log(nu/2)  - n * lgamma(nu/2)
    if(Family == "VG")   
      resp= n * nu * log(nu/2) - n * lgamma(nu)
    if(Family== "NMVBS")
      resp= - n * (log(nu)-1/nu^2) 
    if(Family == "NMVL") 
      resp= n * (2*log(nu) - log(1 + nu))
    if(Family == "NIG") 
      resp= n * nu
    if(Family == "H")
      resp= - n * BesselK(nu, 1, log.val = T)
    return(resp)
  }
  
  f2=function(nu, W, Family=Family){
    if(Family== "GHST") 
      resp= - (0.5 * nu) * (log(W) + 1/W)
    
    if(Family == "VG")   
      resp= nu * (log(W) - W/2)
    
    if(Family== "NMVBS")
      resp= - 0.5 * (W + 1 / W) / nu^2
    
    if(Family == "NMVL") 
      resp= - nu * (W)
    if(Family == "NIG") 
      resp= -0.5 * W * nu^2
    if(Family == "H") 
      resp= -0.5 * nu * (W + 1 / W)
    return(resp)
  }
  
  d=function(nu, prior=prior){
    if (prior=="Hierar_1") {resp= -lambda * nu}
    if (prior=="Exp")      {resp= -hyper * nu}
    if (prior=="Unif")     {resp= 0}
    return(resp)
  }
  
  g= function(nu, W, lambda, Family, prior) {
    n = length(W)
    ff = d(nu,prior) + f1(nu, Family) + sum(f2(nu,W, Family))
    return(ff)
  }
  
  if(Family == "GHST"){
    if(prior=="Hierar_1"){
      Fonseca1 = deriv(~ (-lambda * nu) + 0.5 * n * nu * log(nu/2) - n * lgamma(nu/2),
                       c("nu"), function(nu) {
                       }, hessian = TRUE)
    }
    
    if(prior=="Exp"){
      Fonseca1 = deriv(~ (-hyper * nu) + 0.5 * n * nu * log(nu/2) - n * lgamma(nu/2),
                       c("nu"), function(nu) {
                       }, hessian = TRUE)
    }
    
    if(prior=="Unif"){
      Fonseca1 = deriv(~ 0.5 * n * nu * log(nu/2) - n * lgamma(nu/2),
                       c("nu"), function(nu) {
                       }, hessian = TRUE)
    }
    Fonseca2 = deriv(~ -(0.5 * nu) * (log(W) + 1/W),
                     c("nu"), function(W,nu) {
                     }, hessian = TRUE)
  }
  
  if(Family == "VG"){
    if(prior=="Hierar_1"){
      Fonseca1 = deriv(~ (-lambda * nu) +  n * nu * log(nu/2) - n * lgamma(nu),
                       c("nu"), function(nu) {
                       }, hessian = TRUE)
    }
    
    if(prior=="Exp"){
      Fonseca1 = deriv(~ (-hyper * nu) +  n * nu * log(nu/2) - n * lgamma(nu),
                       c("nu"), function(nu) {
                       }, hessian = TRUE)
    }
    
    if(prior=="Unif"){
      Fonseca1 = deriv(~ n * nu * log(nu/2) - n * lgamma(nu),
                       c("nu"), function(nu) {
                       }, hessian = TRUE)
    }
    Fonseca2 = deriv(~ nu * (log(W) - W/2),
                     c("nu"), function(W,nu) {
                     }, hessian = TRUE)
  }
  
  if(Family == "NMVBS"){
    if(prior=="Hierar_1"){
      Fonseca1 = deriv(~ (-lambda * nu) - n * (log(nu)-1/nu^2),
                       c("nu"), function(nu) {
                       }, hessian = TRUE)
    }
    
    if(prior=="Exp"){
      Fonseca1 = deriv(~(-hyper * nu) - n * (log(nu)-1/nu^2),
                       c("nu"), function(nu) {
                       }, hessian = TRUE)
    }
    
    if(prior=="Unif"){
      Fonseca1 = deriv(~ - n * (log(nu)-1/nu^2),
                       c("nu"), function(nu) {
                       }, hessian = TRUE)
    }
    Fonseca2 = deriv(~ -0.5 * (W + 1 / W) / nu^2,
                     c("nu"), function(W,nu) {
                     }, hessian = TRUE)
  }
  
  if(Family == "NMVL"){
    if(prior=="Hierar_1"){
      Fonseca1 = deriv(~ (-lambda * nu) + n * (2*log(nu) - log(1 + nu)),
                       c("nu"), function(nu) {
                       }, hessian = TRUE)
    }
    
    if(prior=="Exp"){
      Fonseca1 = deriv(~ (-hyper * nu) + n * (2*log(nu) - log(1 + nu)),
                       c("nu"), function(nu) {
                       }, hessian = TRUE)
    }
    
    if(prior=="Unif"){
      Fonseca1 = deriv(~ n * nu * log(nu/2) - n * lgamma(nu),
                       c("nu"), function(nu) {
                       }, hessian = TRUE)
    }
    Fonseca2 = deriv(~ - nu * (W),
                     c("nu"), function(W,nu) {
                     }, hessian = TRUE)
  }
  
  if(Family == "NIG"){
    if(prior=="Hierar_1"){
      Fonseca1 = deriv(~ (-lambda * nu) +  n * nu,
                       c("nu"), function(nu) {
                       }, hessian = TRUE)
    }
    
    if(prior=="Exp"){
      Fonseca1 = deriv(~ (-hyper * nu) +  n * nu,
                       c("nu"), function(nu) {
                       }, hessian = TRUE)
    }
    
    if(prior=="Unif"){
      Fonseca1 = deriv(~ n * nu,
                       c("nu"), function(nu) {
                       }, hessian = TRUE)
    }
    Fonseca2 = deriv(~ -0.5 * nu^2 * (W),
                     c("nu"), function(W,nu) {
                     }, hessian = TRUE)
  }
  
#-------------------------------------------------------------------------------  
#-------------------------------------------------------------------------------
  
  if(Family == "H"){
    R=function(nu) exp(BesselK(nu, 2, log.val = T) - BesselK(nu, 1, log.val = T))
    
    S=function(nu) exp(BesselK(nu, 3, log.val = T) - BesselK(nu, 1, log.val = T)) -
      exp(2*(BesselK(nu, 2, log.val = T) - BesselK(nu, 1, log.val = T)))
    
    f_prime=function(nu) 1 / nu - R(nu)
    
    f_zegond=function(nu) S(nu)- R(nu)/nu - 1/nu^2
    
    if(prior=="Hierar_1"){
      Fonseca1 =function(nu) c(-lambda * nu - n * BesselK(nu, 1, log.val = T), -lambda - n * f_prime(nu),-n * f_zegond(nu))
    }
    
    if(prior=="Exp"){
      Fonseca1 =function(nu) c(-hyper * nu - n * BesselK(nu, 1, log.val = T), - hyper - n * f_prime(nu), -n * f_zegond(nu))
      
    }
    
    if(prior=="Unif"){
      Fonseca1 =function(nu) c(- n * BesselK(nu, 1, log.val = T), - n * f_prime(nu), -n * f_zegond(nu))
      
    }
    
    Fonseca2 = deriv(~ -0.5 * nu * (W + 1 / W),
                     c("nu"), function(W,nu) {
                     }, hessian = TRUE)
  }
  
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
  
  aux1 = Fonseca1(last)  
  aux2 = Fonseca2(W,last)
  if(Family== "H"){
    q1= aux1[2]+ sum(attr(aux2, "gradient"))
    q2= aux1[3]+ sum(attr(aux2, "hessian"))
  }else{
    q1 = attr(aux1, "gradient")[1] + sum(attr(aux2, "gradient"))
    q2 = attr(aux1, "hessian")[1] + sum(attr(aux2, "hessian"))
  }
  
  
  bw = max(0.001, -1/q2)
  aw = last + q1 * bw
  
  if(Family=="GHST") {DB=2.001}else{DB=0.001}
  
  cand = truncnorm :: rtruncnorm(1, a = DB, b = 50, mean = aw, sd = sqrt(bw)) 
  
  aux1 = Fonseca1(cand)
  aux2 = Fonseca2(W, cand)
  
  if(Family== "H"){
    q11=  aux1[2]+ sum(attr(aux2, "gradient"))
    q21=  aux1[3] + sum(attr(aux2, "hessian"))
  }else{
    q11 = attr(aux1, "gradient")[1] + sum(attr(aux2, "gradient"))
    q21 = attr(aux1, "hessian")[1] + sum(attr(aux2, "hessian"))
  }
  
  bw1 = max(0.0001, -1/q21)
  aw1 = cand + q11 * bw1
  
  alfa <- (exp(g(cand, W, lambda,Family, prior) - g(last, W, lambda,Family,prior))) * 
    (truncnorm :: dtruncnorm(last, a = DB, b = 50, mean = aw1, sd = sqrt(bw1))/
       truncnorm :: dtruncnorm(cand, a = DB, b = 50, mean = aw, sd = sqrt(bw)))
  last = ifelse(runif(1) < min(alfa, 1), cand, last)
  return(last)
}
