# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #
# # # introduction  # # # # # #
# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #

  ## parameter, initial, and max time value
     lambda <- 0.3
	 p0 <- 8.0s
     tmax <- 10
  ## true solution
     plot(function(t) p0*exp(lambda*t), lwd=2, col=2,
     xlim=c(0, tmax), xlab='time', ylab='p')

# # # # # # # # # # # # # # # #

  ## length of step
     h <- 0.1
  ## number of steps
     n <- ceiling(tmax/h)
  ## initialize
     ts <- 0.0
	 ps <- p0

# # # # # # # # # # # # # # # #
  ## method
  for (i in 2:n) {
     ## approximations
     ps[i] <- ps[i-1] + h*(lambda*ps[i-1])
     ts[i] <- ts[i-1] + h
}

# # # # # # # # # # # # # # # #

  ## plotting approximates
     lines(ts, ps, type='l', lty=3, lwd=2, pch=19)
  ## adding a legend
     legend('bottomright', c('exact', 'approx'),
     lwd=2, col=c(2, 1), lty=c(1, 3))
     
# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #

  ## load package
  library(deSolve) ## see text
  ## give model
  POPmod <- function(t, x, parms)  {
     with(as.list(c(parms, x)), {
     dP <- lambda*P
     list(c(dP))
})}

# # # # # # # # # # # # # # # #

  ## list parameter(s)
     parms  <- c(lambda = 0.3)
  ## specify times solution should be returned
     times <- seq(0, 10, by=0.01)
  ## specify initial condition(s)
     xstart <- c(P = 8.0)
  ## apply solver
     out <-  lsoda(xstart, times, POPmod, parms)
     out <-  as.data.frame(out)
     
# # # # # # # # # # # # # # # #

  ## plot solution
     plot(out$time, out$P, type='l', col='red',
     lwd=2, ylim=c(0, max(out$P)), las=1)
     
# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #
# # # section 2 # # #
# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #

  # a Predator-Prey model
  PPmod <- function(t, x, parms) {
     with(as.list(c(parms, x)), {
     dH <- alpha*H - beta*H*L
     dL <- gamma*H*L - delta*L
     list(c(dH, dL))
})}

# # # # # # # # # # # # # # # #

  ## parameters
     (parms <- c(alpha=1.1, beta=0.05, gamma=0.01,
     delta=0.5))
  ## vector of time steps
     (times <- seq(0, 500, by=0.1))
  ## initial conditions
     (xstart <- c(H = 50, L = 8))

# # # # # # # # # # # # # # # #

  out <- lsoda(xstart, times, PPmod, parms)
  out <- as.data.frame(out)

# # # # # # # # # # # # # # # #

  ## plot 1
     matplot(out$time, out[, c('H', 'L')], type='l',
     lwd=2, col=c('blue', 'red'), lty=1,
     ylim=c(0, 150), xlab='Time', ylab='Population
     densities (number/area)',  las=1)
     legend('topleft', c('Prey (H)', 'Predator (L)'),
     col=c('blue', 'red'), lty=1, bg='white')
  ## plot 2
     plot(out$H, out$L, type='l', xlim=c(0, 150),
     ylim=c(0, 50), xlab='Prey (H)',
     ylab='Predator (L)', las=1)
     
# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #
# # # # # Cats  # # # # # # # #
# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #

TOXOmod <- function(t, x, parms)  {
   with(as.list(c(parms, x)), {
## susceptible cats
   dCs <- B - Omega*E*Cs - M*Cs
## infected cats
   dCi <- Omega*E*Cs - Gamma*Ci - Delta.i*M*Ci
## chronically infected cats
   dCc <- Gamma*Ci - Delta.c*M*Cc
## environment
   dE <- (Lambda.i*Ci + Lambda.c*Cc) - eta*E
   list(c(dCs, dCi, dCc, dE))
})}

# # # # # # # # # # # # # # # #

  ## The parameters
     parms <- c(
        B = 60*(1/52), ## birth rate of kittens per week
        Omega = 0.046, ## rate of infection
        M = 0.2*(1/52), ## mortality of cats
        Gamma = 11/15,  ## loss of acute infection
        Delta.i = 1.05, ## added mortality from acute
        Delta.c = 1.0005, ## added mortality from chronic
        eta = 4*(1/52), ## oocyst decay rate
        Lambda.i = (0.07/11), ## acute shedding rate
        Lambda.c = 0) ## chronic shedding rate
  ## vector of timesteps
     times <- seq(0, 50*52, by=1)
  ## initial conditions
     xstart <- c(Cs = 495, Ci = 5, Cc = 0, E = 0)
     
# # # # # # # # # # # # # # # #

  ## run model
     out <- lsoda(xstart, times, TOXOmod, parms)
     out <- as.data.frame(out)
     
# # # # # # # # # # # # # # # #

## plot susceptibles
plot(out$time, out$Cs, type='l', col='#B276B2', lty=1, 
   ylim=c(0, 500), las=1, lwd=3, 
   xlab='Time (in weeks)', ylab='Cats')
## add acute then chronic
lines(out$time, out$Ci, col='#F17CB0', lty=3, lwd=3)
lines(out$time, out$Cc, col='#FAA43A', lty=2, lwd=3)
## legend
legend('topright', c('S', 'I', 'C'), lty=c(1, 3, 2), 
   lwd=2, col=c('#B276B2', '#F17CB0', '#FAA43A'))
## acute infection closeup
par(fig = c(0.25, 0.80, 0.5, 0.95), new = T)
plot(out$time, out$Ci, type='l', col='#F17CB0', lty=3, 
   lwd=3, xlim=c(0, 600), ylim=c(0, 35),
   xlab='Time (in weeks)', ylab='Infectious', las=1)

# # # # # # # # # # # # # # # #

plot(out$time, out$E, type='l', col='#60BD68', las=1,
   lwd=3, ylim=c(0, 2), xlab='Time (in weeks)', 
   ylab='Oocysts in Environment (billions)')
## legend
legend('topright', 'E', lwd=2, col='#60BD68')


# # # # # # # # # # # # # # # #

  ## sets range of parameters
     omegas <- seq(0, 0.15, length=101)
  ## initialize 'vals' for storage
     vals <- NULL
     sols <- NULL
# # # # # # # # # # # # # # # #

  ## steps through omegas
  for(i in 1:length(omegas)){
     parms['Omega'] <- omegas[i]
  ## solve and store solution
      out  <- lsoda(xstart, times, TOXOmod, parms)
      sols[[i]] <- as.data.frame(out)
  ## saves last row, 2nd through last column of solution
      vals <- rbind(vals, out[nrow(out), 2:ncol(out)])
}

# # # # # # # # # # # # # # # #

## panel 1
    matplot(omegas, vals[,1:3], type='l', lty=c(1,3,2),
       lwd=2, col=c('#B276B2', '#F17CB0', '#FAA43A'),
       xlab='Rate of infection',
       ylab='Cat steady-state values', las=1)
## panel 2
    plot(omegas, vals[,4], type='l', lty=1, lwd=2,
       col='#60BD68', las=1, xlab='Rate of infection',
       ylab='Oocyst steady-state values (billions)')

# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #
# # # # #  stochastic # # # # #
# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #

  set.seed(100) ## for reproducibility
  ## parameters
     (r <- 0.4)
     (mu <- 0.1)
  ## transition probabilities
     (a <- c(r, mu))
  ## transitions
     (trans <- c(1, -1))
  ## initialization
     (Xs <- X <- 8)
	 (taus <- 0)

# # # # # # # # # # # # # # # #

  for(i in 2:200){
     if(X==0)break;
     (rn <- runif(2))
     (rate <- a*Xs[i-1])
     (pr <- cumsum(rate)/sum(rate))
     (rxn <- min(which(pr > rn[1])))
     ## print(rxn) (useful to debug)
     (X <- X + trans[rxn])
     (Xs <- c(Xs, X))
     (taus[i] <- - log(rn[2])/sum(rate))
}

# # # # # # # # # # # # # # # #

  ## calculate actual time
     (times <- cumsum(taus))
     (tmax <- ceiling(max(times)))
     (xmax <- ceiling(max(Xs)))

# # # # # # # # # # # # # # # #

  plot(function(t)Xs[1]*exp((r-mu)*t), col=2, lwd=2,
     xlim=c(0, tmax), ylim=c(0, xmax),
     xlab='Time', ylab='Population size', las=1)
  points(times, Xs, pch=19, lwd=2)
  
# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #
# # # stochastic predator # # #
# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #

 parms <- c(alpha=1.1, beta=0.05, gamma=0.01, delta=0.5)
  Xs <- X <- t(c(H=50, L=8))
  taus <- 0
  trans <- rbind(c(1, -1, 0, 0), c(0, 0, 1, -1))
  
# # # # # # # # # # # # # # # #

 for(i in 2:10000){
  ## stops at extinction
     if(sum(X) == 0)break;
  ## same as before
     rn <- runif(2, 0, 1)
  ## rate contains more possibilities
     rate <- with(as.list(parms),
        c(alpha*X[1], beta*X[1]*X[2],
        gamma*X[1]*X[2], delta*X[2]))
  ## same as before
     prop <- cumsum(rate)/sum(rate)
     rxn <- min(which(rn[1] < prop))
  ## apply transition and store state
     X <- X + t(as.matrix(trans[, rxn], 2, 1))
     Xs <- rbind(Xs, X) ## other ways to do this
  ## store time increment
     taus[i] <- -log(rn[2])/sum(rate)
}

# # # # # # # # # # # # # # # #
  ## calculate event times
     times <- cumsum(taus)
  ## organize output
     dat <- as.data.frame(cbind(times, Xs))
     names(dat) <- c("time", "H", "L")

# # # # # # # # # # # # # # # #

  ## a fairly generic plot, one solution
  matplot(dat$time, dat[, c("H", "L")], type='l', lwd=1,
     col=c('blue', 'red'), lty=1, xlab="Time",
     ylab="Population densities (number/area)", las=1)

# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #
# # # # circadian # # # # # # #
# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #

 ## the model
  CLOCKmod <- function(t, x, parms)  {
    with(as.list(c(parms, x)), {
      if(t > t.switch1){
      a <- 1 ## alcohol value switch
      }
      if(t > t.switch2){
      a <- 1
      }
  # conc. of Per2/Cry mRNA
      dY1 <- (a*v1b*(Y7+c))/(k1b*(1+(Y3/k1i)^p+Y7+c))-k1d*Y1
  # conc. of PER2/CRY complex in the cytoplasm
      dY2 <- k2b*(Y1)^q - k2d*Y2 - k2t*Y2 + k3t*Y3
  # conc. of PER2/CRY complex in the nucleus
      dY3 <- k2t*Y2 - k3t*Y3 - k3d*Y3
  # conc. of Bmal1 mRNA
      dY4 <- (v4b*(Y3)^r)/((k4b)^r + (Y3)^r) - k4d*Y4
  # conc. of BMAL1 protein in the cytoplasm
      dY5 <- k5b*Y4 - k5d*Y5 - k5t*Y5 + k6t*Y6
  # conc. of BMAL1 protein in the nucleus
      dY6 <- k5t*Y5 - k6t*Y6 - k6d*Y6 + k7a*Y7 - k6a*Y6
  # conc. of transcriptionally active form BMAL1
      dY7 <- k6a*Y6 - k7a*Y7 - k7d*Y7
  # the output
      list(c(dY1, dY2, dY3, dY4, dY5, dY6, dY7))
})}

# # # # # # # # # # # # # # # #

(parms <- c(v1b = 9, k1b = 1, k1i = 0.56, c = 0.01,
   p = 8, k1d = 0.12, k2b = 0.3, q = 2, k2d = 0.05,
   k2t = 0.24, k3t = 0.02, k3d = 0.12, v4b = 3.6,
   k4b = 2.16, r = 3, k4d = 0.75, k5b = 0.24,
   k5d = 0.06, k5t = 0.45, k6t = 0.06, k6d = 0.12,
   k6a = 0.09, k7a = 0.003, k7d = 0.09,
   a = 1, t.switch1 = 5000, t.switch2 = 5000))
times <- seq(0, 5000, by=0.1)
(xstart <- c(Y1 = 0.25, Y2 = 0.3, Y3 = 1.15, Y4 = 0.85,
   Y5 = 0.72, Y6 = 1.35, Y7 = 1.075))

# # # # # # # # # # # # # # # #

  ## calculate solution
     out  <- lsoda(xstart, times, CLOCKmod, parms)
     out  <- as.data.frame(out)
  ## plots
     plot(out$time, out$Y3, type='l', col='gray',
        lty=1, lwd=2, ylim=c(0, 5), xlim=c(0, 250),
        xlab='Time (hours)',
        ylab='Protein amounts (nM)', las=1)
     lines(out$time, out$Y5 + out$Y6 + out$Y7,
        type='l', col='black', lty=2, lwd=2)
        
# # # # # # # # # # # # # # # #

  ## fix parameter value, repeat as necessary
     parms['a'] <- 0.1
     out2  <- lsoda(xstart, times, CLOCKmod, parms)
     out2  <- as.data.frame(out2)
  ## add corresponding solutions
     lines(out2$time, out2$Y3, type='l', col='red',
        lty=1, lwd=2)
     lines(out2$time, out2$Y5 + out2$Y6 + out2$Y7,
        type='l', col='red4', lty=2, lwd=2)

# # # # # # # # # # # # # # # #

 ## add legend
     legend('topleft', c('PER2/CRY Protein',
        'BMAL1 Protein', 'PER2/CRY Protein (w/ alcohol)',
        'BMAL1 Protein (w/ alcohol)'), col=c('gray',
        'black', 'red', 'red4'), lty=c(1, 2, 1, 2), lwd=2)

# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #
# # # # # influenza # # # # # #
# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #

  FLUmod <- function(t, x, parms)  {
    with(as.list(c(parms, x)), {
  #for all, m is parameter for natural mortality
  #Susceptible
     dS <- p - m*S - beta*S*I - zeta*S + omega*(V + N + R)
  #Vaccination
     dV <- zeta*S - m*V - omega*V
  #Infectious
     dI <- beta*S*I - gamma*I - m*I
  #Recovered Naturally
     dN <- (1 - delta)*(1 - rho)*gamma*I - m*N - omega*N
  #Treated, coded as R to avoid confusion with T for TRUE
     dR <- (1 - delta)*(rho)*gamma*I - m*R - omega*R
  #Mortality
     dD <- delta*gamma*I
  #Added cost equation
     dC <- z*zeta*S + b*beta*S*I + d*delta*gamma*I +
       + n*(1 - delta)*(1 - rho)*gamma*I  +
       + r*(1 - delta)*(rho)*gamma*I
      list(c(dS, dV, dI, dN, dR, dD, dC)) # the output
  })}


# # # # # # # # # # # # # # # #
## parameters
   parms <- c(p = (13)*(1/52), m = (0.013)*(1/52),
      beta = (0.08)*(1/52), zeta = (1.2)*(1/52),
      gamma = (12)*(1/52), omega = (0.8)*(1/52),
      delta = 0.0005, rho = (1/3), z = 52.92,
      b = 3.00, d = 542430.70, n = 110.28, r = 5613.65)
## time points
   times <- seq(0, 52*2, by=0.1)
## initial values
   xstart <- parms['p']/parms['m']*c(S = 0.98, V = 0.01,
      I = 0.01, N = 0, R = 0, D = 0, C = 0)

# # # # # # # # # # # # # # # #
out <- as.data.frame(lsoda(xstart, times, FLUmod, parms))
# # # # # # # # # # # # # # # #

  ## panel 1
  matplot(out$time, out[, c('S', 'V', 'I', 'R','N')],
     type='l', lty=1, lwd=2, col=c("red", "gray", "green", "blue", "magenta"),
     xlab='Time (weeks)', ylab='Population densities (number/area)', las=1)
  legend('topright', c('Susceptible', 'Vaccinated',
     'Infectious', 'Treated', 'Natural'),
     col=c("red", "gray", "green", "blue", "magenta"), lty=1, lwd=2)
  ## panel 2 (cost scaled to millions)
  plot(out$time, out$C/1e6, type='l', lty=1,
     xlab='Time (weeks)', ylab='Cost (millions USD)', las=1)
# # # # # # # # # # # # # # # #

  ## time vector
     times <- seq(0, 52*100, by=0.1)
  ## range of zeta values
     zetas <- c(seq(0.0, 0.2, by=0.01))
  ## object to store model output
     vals <- NULL
  for(i in 1:length(zetas)){
     parms['zeta'] <- zetas[i]
     out  <-  lsoda(xstart, times, FLUmod, parms)
     out  <- as.data.frame(out)
     if(i/10==floor(i/10))print(i) ## progress bar
     vals <- rbind(vals, out[nrow(out), ])
}

# # # # # # # # # # # # # # # #

  matplot(zetas, vals[, c('S', 'V', 'I', 'R', 'N')],
     type='l', lty=1, lwd=2, col=c("red", "gray", "green", "blue", "magenta"),
     xlim=c(0, 0.2), ylim=c(0, 1000), xlab='Vaccination
     rate', ylab='Population steady-state values', las=1)
  legend('topleft', c('Susceptible', 'Vaccinated',
     'Infectious', 'Treated', 'Natural'),
     col=c("red", "gray", "green", "blue", "magenta"), lty=1, lwd=2)
     
# # # # # # # # # # # # # # # #

