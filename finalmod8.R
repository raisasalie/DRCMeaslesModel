rm(list = ls(all=TRUE))
setwd("~/Documents/Masters 2020:2021/MMID/Ass1")
####################################
### data
### to be used in in S V E I R (from WHO, 2019)
# Dec2019 
cases <- 3550
#how many in Kinshasa (estimate)
cases <- cases*(14.3/85)
#Vaccination rate - 2008 est
prop.vacc <- 0.92
prop.S <- 1-prop.vacc

### population distribution among patches - randomly allocated 
### spatial layer
n.vills <- 3 
kinshasa.pop <- 14.3*(10^6) - cases # 2020
capital.pop <- 8.9*(10^6)*(1.0433)^5 # incorporate growth from (2015)

vills.pop <- kinshasa.pop - capital.pop

villages <- rep(0, n.vills)
set.seed(2020)
for (i in 1:vills.pop){
  v <- sample(c(1:n.vills),1)
  villages[v] <- villages[v] + 1
}

pop <- c(capital.pop, villages)
names(pop) <- c("cap", "vill1", "vill2", "vill3")
#proportions belonging to each
props <-pop/sum(pop)

### social layer
# proportion children (0-14 years) - from DRC factbook, CIA 
p <- 0.4125
p <- c(p, 1-p)
ps <- matrix(NA, ncol = 4, nrow = 2)
colnames(ps) <- names(pop)
rownames(ps) <- c("children", "adults")
for(j in 1:ncol(ps)){
  for(i in 1:nrow(ps)){
    ps[i,j] <- p[i]*pop[j]
  }
}
ps <- t(ps)
# cases matrix
case.m <- ps/sum(ps)*cases

################################
#### Parameters 
# birth and rate per month (from crude rate)
mu <- 40.6/(1000*12)
# natural death rate per month (from crude)
alph <- 9.3/(1000*12)
# death due to infection
phi <- 0.125/12 #0.125 per year
# vaccine effectiveness 
e <- 0.95
# rate of immunity 
kapp <- 0.14286/12 # 0.14286 per year
# beta
b <- 0.09091
# incubation period (rate of infectiousness)
i <- 12 
eps <- 1/i *30
# recovery rate 
r <- 21-i #days, recovery time  
gamm <- 1/r * 30 
# vaccination coverage in locality i
etas <- c(rep(0.95, 20))
############################ 
# social contact matrix
# all locations, Congo
library(readxl)
contmat <- read_excel("MUestimates_all_locations_1.xlsx", sheet = "Congo")
# upper bound of ages (naming cols and rows)
ages <- seq(from = 5, by = 5, length.out = ncol(contmat))
colnames(contmat) <- ages
#rownames(contmat) <- ages
heatmap(as.matrix(contmat), keep.dendro = F)
# set up final social contact matrix 
sum11 <- sum(contmat[c(1:2),c(1:2)])
sum22 <- sum(contmat[c(3:16),c(3:16)])
sum12 <- sum(contmat[c(1:2),c(3:16)])
sum21 <- sum(contmat[c(3:16),c(1:2)])
# change to monthly
thetasoc <- 30*matrix(c(sum11, sum12, sum21, sum22), ncol = 2, byrow = T)
colnames(thetasoc) <- c("0-15", "16+")
rownames(thetasoc) <- colnames(thetasoc)
##########################################
#spatial contact matrix - simulated
thetaspa <- diag(rep(1,4))
set.seed(2020)
e <- 30*runif(1)
# all connect to 1 - capital
thetaspa[1, 2:4] <- e
thetaspa[2:4, 1] <- e
# 2 and 3 connected 
thetaspa[2,3] <- e
thetaspa[3,2] <- e
# 4 and 3 connected4
thetaspa[4,3] <- e
thetaspa[3,4] <- e
#######################################################
#model
patch.i <- function(t, x, parms){
  ##############################
  
  #####################################
  with(as.list(c(parms, x)),{
    # population of each patch, within each age group
    N11 <- S11+E11+I11+R11+V11
    N12 <- S12+E12+I12+R12+V12
    
    N21 <- S21+E21+I21+R21+V21
    N22 <- S22+E22+I22+R22+V22  
    
    N31 <- S31+E31+I31+R31+V31
    N32 <- S32+E32+I32+R32+V32
    
    N41 <- S41+E41+I41+R41+V41
    N42 <- S42+E42+I42+R42+V42
    
    # proportion which are children 
    #pr <- (N11+N21+N31+N41)/(N11+N21+N12+N22+N31+N32+N41+N42)
    # for both age groups
    #pr <- c(pr, 1-pr)
    
    # divide into components
    N <- matrix(c(N11,N12,N21,N22,N31,N32, N41, N42), ncol = 2, nrow = 4, byrow = T)
    S <- matrix(c(S11,S12,S21,S22,S31,S32, S41, S42), ncol = 2, nrow = 4, byrow = T)
    V <- matrix(c(V11,V12,V21,V22,V31,V32, V41, V42), ncol = 2, nrow = 4, byrow = T)
    E <- matrix(c(E11,E12,E21,E22,E31,E32, E41, E42), ncol = 2, nrow = 4, byrow = T)
    I <- matrix(c(I11,I12,I21,I22,I31,I32, I41, I42), ncol = 2, nrow = 4, byrow = T)
    R <- matrix(c(R11,R12,R21,R22,R31,R32, R41, R42), ncol = 2, nrow = 4, byrow = T)
    
    out <- c()
    #mat <- matrix(NA)
    # for each patch i
    for (i in 1:4){
      # for each social group j
      for (j in 1:2){
        #set up indicator var
        for (j in 1:2){
          if(j==1){ind = 1}else{ind=0}
        }
        
        dS <- ind*mu*(1-eta[i])*N[i,j] - b*I[i,j]*S[i,j]*sum(thetasoc[,j])/N[i,j] + sum(S[,j]*thetaspa[,i]) - S[i,j]*sum(thetaspa[i,])
        dV <- ind*mu*eta[i]*N[i,j] - alph*V[i,j] - kapp*V[i,j] + sum(V[,j]*thetaspa[,i]) - V[i,j]*sum(thetaspa[i,])
        dE <- b*I[i,j]*S[i,j]*sum(thetasoc[,j])/N[i,j] - (eps+alph)*E[i,j] + sum(E[,j]*thetaspa[,i]) - E[i,j]*sum(thetaspa[i,])
        dI <- eps*E[i,j] - (gamm+alph+phi)*I[i,j] + sum(I[,j]*thetaspa[,i]) - I[i,j]*sum(thetaspa[i,])
        dR <- gamm*I[i,j] + kapp*V[i,j] - alph*R[i,j] + sum(R[,j]*thetaspa[,i]) - R[i,j]*sum(thetaspa[i,])
        
        # order of out: SEIRV (as with start)
        out <- c(out, dS, dE, dI, dR, dV)
        
        # set negs to zero
        for (k in 1:length(out)){
          if (out[k]>0){}else{out[k]<-0}
        }
      }
    }
    list(out)
  })
}

eta <- rep(0.95, 4)
parms <- c(mu=mu, alph=alph, phi=phi, kapp=kapp, b=b, 
           eps=eps, gamm=gamm, eta=eta, 
           thetasoc=thetasoc, thetaspa=thetaspa)
# double ratio for children, and patch 1
r <- 1/14
start <- c(S11=prop.S*ps[1,1], E11=0, I11=cases*r*4, R11=prop.vacc*ps[1,1], V11=0,
           S12=prop.S*ps[1,2], E12=0, I12=cases*r, R12=prop.vacc*ps[1,2], V12=0,
           S21=prop.S*ps[2,1], E21=0, I21=cases*r*2, R21=prop.vacc*ps[2,1], V21=0,
           S22=prop.S*ps[2,2], E22=0, I22=cases*r, R22=prop.vacc*ps[2,2], V22=0,
           S31=prop.S*ps[3,1], E31=0, I31=cases*r*2, R31=prop.vacc*ps[3,1], V31=0,
           S32=prop.S*ps[3,2], E32=0, I32=cases*r, R32=prop.vacc*ps[3,2], V32=0,
           S41=prop.S*ps[4,1], E41=0, I41=cases*r*2, R41=prop.vacc*ps[4,1], V41=0,
           S42=prop.S*ps[4,2], E42=0, I42=cases*r, R42=prop.vacc*ps[4,2], V42=0)
# set seq of times
# monthly for next three years 
times <- seq(0, 24, 1)
# solve DE's
library(deSolve)
run <- ode(y = start, times = times, func = patch.i, parms=parms)
# divide into patches 
patches <- list(NA, NA, NA, NA) 
names(patches) <- c("P1", "P2", "P3", "P4")
for(i in 1:4){ 
  # for each patch
  res <- matrix(NA, ncol = (ncol(run))/4+1, nrow = nrow(run))
  res[,c(1:11)] <- run[,c(1,c((10*i-8):(10*i+1)))]
  colnames(res) <- rep(NA, 11)
  colnames(res)[c(1:11)] <- colnames(run)[c(1,c((10*i-8):(10*i+1)))]
  
  # add totals 
  Ni <- rowSums(res[,-1])
  Ni1 <- rowSums(res[, c(2:6)])
  Ni2 <- Ni-Ni1
  res <- cbind(res, Ni1, Ni2, Ni)
  colnames(res)[(ncol(res)-2):ncol(res)] <- paste("N", i, c(1,2,"all"), sep="")
  
  # add to list
  patches[[i]] <- res 
}

# plot 
cols <- c("red", "green", "purple", "pink", 
          "orange", "black", "darkblue", "darkred", "darkgreen")

for (k in 1:4){ #k=patch
  #
  #pdf(file = paste("patch",k,".pdf", sep=""), width = 12, compress = F)
  soln <- patches[[k]]
  
  # in tens of thousands
  soln <- soln/(10^4)
  
  par(mfrow=c(1,3))
  for (i in 1:2){#i=age group
    run <- soln[,c((5*i-3):(5*i+1))][,-c(1,4,5)]
    #
    # plot SEI - others not of interest
    plot(run[,1], type = "l", lwd = 2, col = cols[1], 
         ylim=c(0, max(run)), xlim=c(0,24),
         ylab ="Number of people ('0 000)", xlab="Month")
    for (j in 2:(ncol(run))){
      lines(run[,j], col = cols[j], lwd = 2)
    }
    legend("bottomright", legend = colnames(run), col = cols, lwd = 2)
    
  }
  
  #plot N's and R 
  run <- soln[ ,c(5,10,(ncol(soln)-2):ncol(soln))]
  plot(run[,1], type = "l", lwd = 2, col = cols[1], 
       ylim=c(0, max(run)), xlim=c(0,24),
       ylab ="Number of people ('0 000)", xlab="Month")
  for (j in 2:ncol(run)){
    lines(run[,j], col = cols[j], lwd = 2)
  }
  legend("bottomright", legend = colnames(run), col = cols, lwd = 2)
  #dev.off()
}

########################################
# NLP
#matrix(0.95, ncol =4, nrow=1)
runmod <- function(eta){
  #eta = array of 4
  parms <- c(mu=mu, alph=alph, phi=phi, kapp=kapp, b=b, 
             eps=eps, gamm=gamm, eta=eta, 
             thetasoc=thetasoc, thetaspa=thetaspa)
  run <- ode(y = start, times = times, func = patch.i, parms=parms)
  
  return(run)
}

#objective function
eval_f <- function(x){
  Is <- runmod(eta=x)[, c(4,9,14,19,24,29,34,39)]
  
  return(list("objective" = max(rowSums(Is)), "gradient"=rep(0,4)))
  #
}

set.seed(2020)
x0 <- runif(4)

#lower/upper bounds
lb <- rep(0,4)
ub <- rep(1,4)

# available vaccines
vmax <- 8*(10^5)


#contraints
eval_g_ineq <- function(x){
  constr <- -mu*24*t(as.matrix(pop))%*%as.matrix(x)+vmax
  
  grad <- c(x[2]+x[3]+x[4],
            x[1]+x[3]+x[4],
            x[1]+x[2]+x[4],
            x[1]+x[2]+x[3])
  
  return(list("constraints"=constr, "jacobian" = grad))
}
#,

library(nloptr)

#local_opts <- list( "algorithm" = "NLOPT_LD_LBFGS", "xtol_rel"  = 1.0e-7 )
opts <- list( "algorithm" = "NLOPT_LN_COBYLA",
              "xtol_rel"  = 1.0e-2,
              "maxeval"   = 500)
set.seed(2020)
res <- nloptr(x0=x0,
              eval_f=eval_f,
              lb=lb,
              ub=ub,
              eval_g_ineq=eval_g_ineq,
              opts=opts)

# solutions
eta.opt <- res$solution
names(eta.opt) <- c("Patch 1 (capital city)", paste("Patch", c(2:4), sep=" "))
library(xtable)
t <- as.matrix(eta.opt)
colnames(t) <- c("$eta_i$")
xtable(t(t), digits = 4)
res$objective
res$num_constraints_ineq

########################################################
######## plots from optimal eta
parms <- c(mu=mu, alph=alph, phi=phi, kapp=kapp, b=b, 
           eps=eps, gamm=gamm, eta=eta.opt, 
           thetasoc=thetasoc, thetaspa=thetaspa)

# double ratio for children, and patch 1
r <- 1/14
start <- c(S11=prop.S*ps[1,1], E11=0, I11=cases*r*4, R11=prop.vacc*ps[1,1], V11=0,
           S12=prop.S*ps[1,2], E12=0, I12=cases*r, R12=prop.vacc*ps[1,2], V12=0,
           S21=prop.S*ps[2,1], E21=0, I21=cases*r*2, R21=prop.vacc*ps[2,1], V21=0,
           S22=prop.S*ps[2,2], E22=0, I22=cases*r, R22=prop.vacc*ps[2,2], V22=0,
           S31=prop.S*ps[3,1], E31=0, I31=cases*r*2, R31=prop.vacc*ps[3,1], V31=0,
           S32=prop.S*ps[3,2], E32=0, I32=cases*r, R32=prop.vacc*ps[3,2], V32=0,
           S41=prop.S*ps[4,1], E41=0, I41=cases*r*2, R41=prop.vacc*ps[4,1], V41=0,
           S42=prop.S*ps[4,2], E42=0, I42=cases*r, R42=prop.vacc*ps[4,2], V42=0)

# set seq of times
# monthly for next two years 
times <- seq(0, 24, 1)
# solve DE's
library(deSolve)
sim <- ode(y = start, times = times, func = patch.i, parms=parms)

######################################
# divide into patches 
to.patch <- function(sim){
  patches <- list(NA, NA, NA, NA) 
  names(patches) <- c("P1", "P2", "P3", "P4")
  for(i in 1:4){ 
    # for each patch
    res <- matrix(NA, ncol = (ncol(sim))/4+1, nrow = nrow(sim))
    res[,c(1:11)] <- sim[,c(1,c((10*i-8):(10*i+1)))]
    colnames(res) <- rep(NA, 11)
    colnames(res)[c(1:11)] <- colnames(sim)[c(1,c((10*i-8):(10*i+1)))]
    
    # add totals 
    Ni <- rowSums(res[,-1])
    Ni1 <- rowSums(res[, c(2:6)])
    Ni2 <- Ni-Ni1
    res <- cbind(res, Ni1, Ni2, Ni)
    colnames(res)[(ncol(res)-2):ncol(res)] <- paste("N", i, c(1,2,"all"), sep="")
    
    # add to list
    patches[[i]] <- res 
  }
  
  return(patches)
}


# plot 
cols <- c("red", "green", "purple", "pink", 
          "orange", "black", "darkblue", "darkred", "darkgreen")

plot.patch <- function(soln){
  soln <- soln/(10^4)
  
  # plot EIRV
  #par(mfrow=c(1,3))
  for (i in 1:2){#i=age group
    run <- soln[,c((5*i-3):(5*i+1))][,-c(1,4)]
    #
    # plot SEI - others not of interest
    plot(run[,1], type = "l", lwd = 2, col = cols[1], 
         ylim=c(0, max(run)), xlim=c(0,24),
         ylab ="Number of people ('0 000)", xlab="Month", 
         main=paste("j = ", i, sep=" "))
    for (j in 2:(ncol(run))){
      lines(run[,j], col = cols[j], lwd = 2)
    }
    legend("right", legend = colnames(run), col = cols, lwd = 1, 
           cex=.75, box.lwd = .5, 
           xjust = 0, x.intersp = .5, y.intersp = .75)
    
  }
  
  #plot N's and R 
  run <- soln[ ,c(2,7,5,10,(ncol(soln)-2):ncol(soln))]
  plot(run[,1], type = "l", lwd = 2, col = cols[1], 
       ylim=c(0, max(run)), xlim=c(0,24),
       ylab ="Number of people ('0 000)", xlab="Month")
  for (j in 2:ncol(run)){
    lines(run[,j], col = cols[j], lwd = 2)
  }
  legend("topleft", legend = colnames(run), col = cols, lwd = 1, 
         cex=0.75, xjust = 0, x.intersp = .5, y.intersp = .75)
  #dev.off()
}

for (k in 1:4){ #k=patch
  pdf(file = paste("patch",k,".pdf", sep=""), width = 12, height = 4, compress = F)
  soln <- to.patch(sim)[[k]]
  
  par(mfrow=c(1,3))
  plot.patch(soln = soln)
  
  dev.off()
}

#################################
# Sensitivity - b
#function for running ODE's, plotting
run.sim <-function(parms){
  # run
  sim <- ode(y = start, times = times, func = patch.i, parms=parms)
  
  #plot each patch
  for (k in 1:4){
    plot.patch(soln = to.patch(sim)[[k]])
  }
}

#bseq <- seq(from=0, to=0.9, length.out = 10)
bseq <- c(0.1,0.9)
for (v in 1:length(bseq)){
  # define new parms 
  parms1 <- c(mu=mu, alph=alph, phi=phi, kapp=kapp, b=bseq[v], 
              eps=eps, gamm=gamm, eta=eta.opt, 
              thetasoc=thetasoc, thetaspa=thetaspa)
  
  pdf(file = paste("beta",v, ".pdf", sep = ""), compress = F, width = 7)
  #dev.new(width=7, height=6, unit="in")
  par(mfrow=c(4,3), cex=.5)
  run.sim(parms1)
  dev.off()
  Sys.sleep(0.5)
}

##############################
# Sensitivity  - gamm, eps
# incubation period (rate of infectiousness)
iseq <- c(10:12) 
epsseq <- rep(0,3)
for(i in 1:length(iseq)){
  epsseq[i] <- 1/iseq[i] *30
}
# recovery rate 
rseq <- rep(21,3)-epsseq #days, recovery time  
gammseq <- rep(0,3)
for(i in 1:length(gammseq)){
  gammseq[i] <- 1/rseq[i] * 30 
}  

### plots for each
for (v in 1:length(gammseq)){
  # define new parms 
  parms1 <- c(mu=mu, alph=alph, phi=phi, kapp=kapp, b=b, 
              eps=epsseq[v], gamm=gammseq[v], eta=eta.opt, 
              thetasoc=thetasoc, thetaspa=thetaspa)
  
  pdf(file = paste("gammeps",v, ".pdf", sep = ""), compress = F, width = 7)
  #dev.new(width=7, height=6, unit="in")
  par(mfrow=c(4,3), cex=.5)
  run.sim(parms1)
  dev.off()
  Sys.sleep(0.5)
}

