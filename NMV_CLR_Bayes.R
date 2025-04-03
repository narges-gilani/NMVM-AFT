#------------------------------------NMV.CLR.Bayes------------------------------

NMV.CLR.Bayes = function(y,
                         x,
                         Family = c("AL", 
                                    "GHST", 
                                    "VG", 
                                    "H", 
                                    "Normal", 
                                    "NMVBS",
                                    "NMVL",
                                    "NIG"),
                         cens = c("Left", "Right"),
                         cc,
                         influence = FALSE,
                         spacing = NULL,
                         hyper_set= NA,
                         prior = c("Exp",
                                   "Hierar_1",
                                   "Unif"),
                         
                         hyper = hyper,
                         n.thin = n.thin,
                         burnin = burnin,
                         n.iter = n.iter,
                         n.chains = n.chains,
                         sample.out = TRUE,
                         table = NULL,
                         level = 95,
                         Prog.Bar = TRUE,
                         Print = FALSE,
                         plot = FALSE)
{
  packages <- c("cli", "knitr")
  
#----------------Install packages not yet installed-----------------------------
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages])
  }
  
#---------------------------Packages loading------------------------------------
  invisible(lapply(packages[-1], library, character.only = TRUE))
  massage = cli::combine_ansi_styles("red", "bold")
  
#-------------------------------------------------------------------------------
#------------------------------- Required Functions ----------------------------
#-------------------------------------------------------------------------------
  
  hpd = function(x, alpha) {
    n = length(x)
    x = sort(x)
    m = floor(alpha * n)
    i = 1:(n - m)
    L = x[i + m] - x[i]
    best = which.min(L)
    structure(c(x[best], x[best + m]), names = c("Lower Bound", "Upper Bound"))
    #return(c(lower = x[best], upper = x[best + m]))
  }
  
  criteria = function(y,
                      x,
                      espac = 20,
                      out,
                      Family,
                      influence,
                      cens = c("Left", "Right")) {
    n = nrow(x)
    p = ncol(x)
    
    n.iter = length(out$sigma2)/espac
    likapprox = matrix(0, n, n.iter)
    inner.CPO = matrix(0, n, n.iter)
    
    beta = colMeans(out$beta)
    sigma2 = mean(out$sigma2)
    tau = mean(out$tau)
    mu = x %*% (beta)
    nu = mean(out$nu)
    
    loglike = Log.like(cc, y, mu, sigma2, tau, nu, Family = Family, cens = cens, result = "Log.like")
    
    for (k in 1:n.iter) {
      beta = out$beta[k, ]
      sigma2 = out$sigma2[k]
      tau = out$tau[k]
      mu = out$mu[k, ]
      nu = out$nu[k]
      likapprox[, k] = Log.like(cc, y, mu, sigma2, tau, nu, Family = Family, cens = cens, result = "PDF")
      
      inner.CPO[, k] = 1 / likapprox[, k]
    }
    
    n.para = p + 3  
    if (Family == "Normal") n.para = p + 1
    if (Family == "AL") n.para = p + 2
    
    mean.log.p = rowMeans(log(likapprox))
    mean.p = rowMeans(likapprox)
    
    var.log <- matrix(0, n, n.iter)
    for (k in 1:n.iter) {
      var.log[, k] <- (log(likapprox[, k]) - mean.log.p)^2
    }
    var.log.p <- rowSums(var.log)/(n.iter - 1)
    
    CPO = 1/(rowMeans(inner.CPO))
    LPML = sum(log(CPO))
    
    pdic = 2 * (loglike - sum(mean.log.p))
    DIC = 2 * pdic - 2 * loglike
    EAIC = -2 * loglike + 2 * n.para
    EBIC = -2 * loglike + log(n) * n.para
    
    WAIC1 = -2 * (sum(log(mean.p)) - 2 * sum(log(mean.p) - mean.log.p))
    WAIC2 = -2 * (sum(log(mean.p)) - sum(var.log.p))
    
    if (influence == "FALSE") {
      return(list(CPO = CPO, LPML = LPML, DIC = DIC, EAIC = EAIC, EBIC = EBIC, WAIC1 = WAIC1,
                  WAIC2 = WAIC2))
    } else {
      aa <- exp(log(rowMeans(inner.CPO)) + log(likapprox))
      IL <- rowMeans(log(aa))
      JL <- rowMeans((aa - 1) * log(aa))
      LL <- rowMeans(abs(aa - 1))
      CHL <- rowMeans((aa - 1)^2)
      
      return(list(CPO = CPO, LPML = LPML, DIC = DIC, EAIC = EAIC, EBIC = EBIC, WAIC1 = WAIC1,
                  WAIC2 = WAIC2, KL = IL, JDist = JL, LDist = LL, ChiDist = CHL))
    }
  }
  
  summary = function(out,
                     level = 95,
                     digits = 2,           #round to digits
                     table  = NULL,        #create html or latex or pandoc table
                     cap    = NULL,        #table caption
                     Family = Family,
                     Print = FALSE,
                     time,
                     cens = c("Left", "Right")
  ){
    
    p = dim(out$beta)[2]
    name2 = c()
    for (k in 1:p) {
      name2[k] = paste("beta", k, sep = "")
    }
    names = c(name2, "sigma2", "tau", "nu")
    
    #get summary statistics
    beta = colMeans(out$beta)
    sigma2 = mean(out$sigma2)
    tau = mean(out$tau)
    nu = mean(out$nu)
    
    c.point = c(beta, sigma2, tau, nu)
    
    se.beta = apply(out$beta, 2, sd)
    se.sigma2 = sd(out$sigma2)
    se.tau = sd(out$tau)
    se.nu = sd(out$nu)
    
    c.sd = c(se.beta, se.sigma2, se.tau, se.nu)
    
    HPD.sigma2 = hpd(out$sigma2, alpha = level/100)
    HPD.tau = hpd(out$tau, alpha = level/100)
    HPD.nu = hpd(out$nu, alpha = level/100)
    
    HPD.beta = matrix(0, nrow = p, ncol = 2)
    for (i in 1:p)
    {
      HPD.beta[i,] = hpd(out$beta[, i], alpha = level/100)
    }
    
    HPD.out = rbind(HPD.beta, HPD.sigma2, HPD.tau, HPD.nu)
    
    #check which coefficients are "significant"
    c.sign = ifelse((HPD.out[, 2] > 0 & HPD.out[, 1] > 0) | (HPD.out[, 2] < 0 & HPD.out[, 1] < 0), "*", "")
    
    #make everything into table
    coefs = format(round(cbind(c.point, c.sd, HPD.out), digits), nsmall = digits)
    coefs = cbind(coefs, c.sign)
    
    hpd.upper = paste0("UB HPD(", level, "%)")
    hpd.lower = paste0("LB HPD(", level, "%)")
    
    header = c("Mean", "SD", hpd.lower, hpd.upper, paste0(level, "% CI excl. 0"))
    
    df = data.frame(coefs)
    rownames(df) = names
    colnames(df) = header
    
    if (is.null(table)) {
      result = kable(df, align = c("c","c","c","c","c"), row.names = T)
    } else {
      if (table == "latex") {
        result = kable(df[,-ncol(df)], align = c("c","c","c","c","c"), row.names = T,
                       format = "latex",
                       booktabs = T,
                       linesep = "",
                       caption = cap)
      }
      if (table == "pandoc") {
        result = kable(df[, -ncol(df)], align = c("c","c","c","c","c"), row.names = T,
                       format = "pandoc",
                       caption = cap)
      }
    }
    
    timecat = ifelse(time > 300,
                     paste0(round(time / 60, 2), " CPU minutes."),
                     paste0(round(time, 2), " CPU seconds."))
    
    if (Print) {
      cat(cli::rule(line = 2, line_col = "orange"), "\n")
      cat(cli::rule(
        center = cli::col_green(
          '** Bayesian results of ', cens, ' censored ', Family, ' linear regression model **'
        ),
        line_col = "deeppink3",
        line = 2
      ),
      "\n")
      cat(cli::rule(
        center = cli::col_green(
          '** Analysis based on ', n.iter, ' posterior draws after a burn-in period of ',
          burnin, ' iterations **'
        ),
        line_col = "deeppink3",
        line = 2
      ),
      "\n")
      cat(cli::rule(
        center = cli::col_green(
          '** MCMC sampling took a total of ', timecat, ' **'
        ),
        line_col = "deeppink3",
        line = 2
      ),
      "\n")
      cat(cli::rule(line = 2, line_col = "orange"), "\n")
    }
    
    result = list(result = result, df = df)
    return(result)
  }
  
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
  
  if (influence == T) {
    if (is.null(spacing)) {
      stop(massage("Spacing must be specified if influence or criteria is TRUE"))
    }
  }
  
  fam = c("AL", 
          "GHST", 
          "VG", 
          "H", 
          "Normal", 
          "NMVBS",
          "NMVL",
          "NIG")
  
  if (all(Family != fam)) {
    stop(massage("Family is not specified correctly. Check documentation!"))
  }
  
  if ((burnin >= n.iter) || burnin < 1)
  {
    stop(massage("Invalid burnin number"))
  }
  if (!is.numeric(n.iter) || n.iter < 1)
  {
    stop(massage("Invalid number of iterations"))
  }
  if (n.thin >= n.iter)
  {
    stop(massage("Invalid number of lag (n.thin) for posterior sample"))
  }
  
#-------------------------------------------------------------------------------
#----------------------------Start of the algorithm-----------------------------
#-------------------------------------------------------------------------------
  begin = proc.time()[1]
  
  if (Print) {
    cat(cli::rule(line = 2, line_col = "orange"), "\n")
    cat(cli::rule(
      center = cli::col_green(
        "* Bayesian estimation of ",
        cens, ' censored ', Family, " linear regression model *"
      ),
      line_col = "deeppink3",
      line = 2
    ),
    "\n")
  }
  
  if (prior == "Exp" || prior == "Hierar") {
    eta1 = 0; eta2 = 0
  }
  
  # for example 
  if (prior == "Gamma") {
    eta1 = eta[1]; eta2 = eta[2]
  }
  
  M <- n.iter
  n = dim(x)[1]
  p = dim(x)[2]
  
  out = Gibbs.sampleing(cc,
                        y,
                        x,
                        n.iter = n.iter,
                        n.thin = n.thin,
                        burnin = burnin,
                        cens = cens,
                        hyper_set = hyper_set,
                        prior = prior,
                        hyper = hyper,
                        Family = Family,
                        n.chains = n.chains,
                        Prog.Bar = Prog.Bar,
                        Print = Print)
  
  end = proc.time()[1]
  time = end - begin
  
  Out1 = summary(out = out,
                 level = level,
                 digits = 3,
                 table = table,
                 cap = " ",
                 Family = Family,
                 Print = Print,
                 time,
                 cens = cens)
  
  if (Print) print(Out1$result)
  para.out = as.data.frame(apply(Out1$df[, -5], 2, as.numeric))
  
  if (influence == FALSE) {
    crit = criteria(y, x, espac = spacing, out, Family, influence, cens = cens)
    KL = JDist = LDist = ChiDist = NULL
  }
  if (influence == TRUE) {
    crit = criteria(y, x, espac = spacing, out, Family, influence, cens = cens)
    KL = crit$KL
    JDist = crit$JDist
    LDist = crit$LDist
    ChiDist = crit$ChiDist
  }
  critFin =
    c(crit$LPML,
      crit$DIC,
      crit$EAIC,
      crit$EBIC,
      crit$WAIC1,
      crit$WAIC2)
  critFin = round(t(as.matrix(critFin)), digits = 3)
  dimnames(critFin) = list(c("Value"), c("LPML", "DIC", "EAIC", "EBIC", "WAIC1", "WAIC2"))
  
  if (Print) {
    cat('\r \n')
    cat(cli::rule(line = 2, line_col = "orange"), "\n")
    cat(cli::rule(
      center = cli::col_green("** Model selection criteria **"),
      line_col = "deeppink3",
      line = 2
    ),
    "\n")
    DF = data.frame(critFin)
    rownames(DF) = "Value"
    colnames(DF) = c("LPML", "DIC", "EAIC", "EBIC", "WAIC1", "WAIC2")
    print(kable(DF, align = c("c","c","c","c","c","c"), row.names = T))
    cat('\r \n')
    cat(cli::rule(line = 2, line_col = "orange"), "\n")
  }
  
  if (plot) {
    library(ggplot2)
    library(gridExtra)
    
    # Dummy data
    datasig = data.frame(x = 1:length(out$sigma2), y = out$sigma2 )
    p.sig = ggplot(datasig, aes(x = x, y = y)) +
      geom_line() + ylab("") + ggtitle(expression(sigma^2)) +
      xlab("Iteration")
    
    if (Family !=  "Normal"){
      datalam = data.frame( x = 1:length(out$tau), y = out$tau )
      p.lam = ggplot(datalam, aes(x = x, y = y)) +
        geom_line() + ylab("") + ggtitle(expression(tau)) +
        xlab("Iteration")
    }
    
    if (Family != "AL" && Family !=  "Normal") {
      datanu = data.frame( x = 1:length(out$nu), y = out$nu )
      p.nu = ggplot(datanu, aes(x = x, y = y)) +
        geom_line() + ylab("") + ggtitle(expression(nu)) +
        xlab("Iteration")
    }
    
    if (Family ==  "Normal") {
      plist = list(p.sig)
      ii = 1
    }
    if (Family ==  "AL") {
      plist = list(p.lam, p.sig)
      ii = 2
    }
    if (Family !=  "Normal" && Family !=  "AL") {
      plist = list(p.lam, p.sig, p.nu)
      ii = 3
    }
    
    for (i in 1:dim(x)[2]) {
      dataname = paste("databeta", i, sep = "")
      pname = paste("p.beta", i, sep = "")
      x = paste("beta[", i, "]", sep = "")
      tname = parse(text = x)
      
      assign(dataname, data.frame( x = 1:length(out$beta[, i]), y = out$beta[, i] ))
      plist[[ii + i]] = ggplot(get(paste("databeta", i, sep = "")), aes(x = x, y = y)) +
        geom_line() + ylab("") + ggtitle(tname) + xlab("Iteration")
    }
    
    nCol = floor(sqrt(length(plist))) + 1
    #do.call("grid.arrange", c(plist, ncol = nCol))
    
    if (menu(c("Yes", "No"), title = "Do you want to draw posterior histograms?") == 1) {
      if (Family !=  "Normal") par(mfrow = c(3, ceiling((p + 3)/2)))
      if (Family ==  "Normal") par(mfrow = c(3, ceiling((p + 1)/2)))
      
      hist(out$sigma2, xlab = expression(sigma^2), prob = T, main = " ")
      if (Family !=  "Normal") {
        hist(out$tau, xlab = expression(tau), prob = T, main = " ")
      }
      if (Family !=  "Normal" && Family !=  "AL") {
        hist(out$nu, xlab = expression(nu), prob = T, main = " ")
      }
      for (i in 1:p) {
        xx = paste("beta[", i, "]", sep = "")
        tname = parse(text = xx)
        hist(out$beta[, i], xlab = tname, prob = T, main = " ")
      }
    }
  }
  
  if (sample.out == TRUE) {
    return(
      list(
        para = para.out,
        tau = out$tau,
        sigma2 = out$sigma2,
        beta = out$beta,
        lambda=out$lambda,
        nu = out$nu,
        KL = KL,
        JDist = JDist,
        LDist = LDist,
        ChiDist = ChiDist,
        Model.selection = critFin
      )
    )
  }
  
  if (sample.out == FALSE) {
    return(
      list(
        para = para.out,
        KL = KL,
        JDist = JDist,
        LDist = LDist,
        ChiDist = ChiDist,
        Model.selection = critFin
      )
    )
  }
}

#-------------------------------------------------------------------------------


