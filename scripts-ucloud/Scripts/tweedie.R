library(tweedie)
library(statmod)

#setwd("")
dat = read.table(file = "input_110-2010-step500_10sims.txt", header = T, sep = "\t")

# set gen to either of 110,510,1010,1510,2010
gen = 510

# Test generation age for different simulations (1-10) of 'gen'
for (sim in 1:10){
    arm.lengths = dat[dat$Age==gen & dat$sim_num==sim, "arm.lengths"]
    length(arm.lengths) # 100 samples, 200 observations (both sides of mut)
    
    # Simple method
    print(paste("Simple method:", 100/mean(arm.lengths)))

    # out <- tweedie.profile(arm.lengths ~ 1, p.vec=seq(1.5, 2.5, by=0.2))
    # out$p.max
    # out$ci
    # 
    # summary(glm(arm.lengths ~ 1, family = tweedie(var.power = 2, link.power = 0)))
    
    ## Fitting Tweedie (and gamma) distributions
    ll.tweedie <- function(par, y) {
        print(par)
        #out <- sum(log(dtweedie(y = y, mu = par[1], phi = exp(par[2]), power = par[3])))
        out <- sum(log(dtweedie(y = y, mu = 2/par[1], phi = 1/2, power = par[2])))
        #out <- sum(log(tweedie.dev(y = y, mu = 2/par[1], power = par[2])))
        #out <- sum(log(dtweedie(y = y, mu = 2/par[1], phi = exp(par[2]), power = par[3])))
        return(out)
    }
    ll.gamma <- function(par, y) {
        print(par)
        out <- sum(log(dtweedie(y = y, mu = 2/par[1], phi = exp(par[2]), power = 2)))
        #out <- sum(log(dtweedie(y = y, mu = 2/par[1], phi = 1/2, power = 2)))
        return(out)
    }
    
    ll.tweedie(par = c(2/mean(arm.lengths), 1.5), y = arm.lengths)
    tt.tweedie = optim(par = c(2/mean(arm.lengths), 1.5), fn = ll.tweedie, y = arm.lengths, control = list(fnscale = -1), hessian = TRUE,
                       lower = c(0.01,1.01), upper = c(100,1.99), method = "L-BFGS-B" )
    tt.gamma = optim(par = c(2/mean(arm.lengths),1/2), fn = ll.gamma, y = arm.lengths, control = list(fnscale = -1), hessian = TRUE,
                     lower = c(0.01,0.01), upper = c(1000, 10), method = "L-BFGS-B" ) ##method = "Brent", 
    
    std.error.tw = sqrt(diag(solve(-tt.tweedie$hessian)))
    std.error.gm = sqrt(diag(solve(-tt.gamma$hessian)))
    
    ll.tweedie(par = tt.tweedie$par, y = arm.lengths) ### Tweedie model
    ll.gamma(par = tt.gamma$par, y = arm.lengths) ### Gamma model
    -2*(ll.tweedie(par = tt.tweedie$par, y = arm.lengths) - ll.gamma(par = tt.gamma$par, y = arm.lengths)) # Check for difference
    
    # Print model estimates and CI
    data.frame("Ic_Min"  = tt.tweedie$par - qnorm(0.975) *std.error.tw, "Estimates" = tt.tweedie$par, "Ic_Max" = tt.tweedie$par + qnorm(0.975) *std.error.tw)
    data.frame("Ic_Min"  = tt.gamma$par - qnorm(0.975) *std.error.gm, "Estimates" = tt.gamma$par, "Ic_Max" = tt.gamma$par + qnorm(0.975) *std.error.gm)

    # Print Tweedie age
    print(data.frame("Ic_Min" = 100/(2/(tt.tweedie$par[1] - qnorm(0.975) *std.error.tw[1])),
                    "Age_estimate" = 100/(2/tt.tweedie$par[1]),
                    "Ic_Max" = 100/(2/(tt.tweedie$par[1] + qnorm(0.975) *std.error.tw[1])))
    )
    
    # Print Gamma age
    #print(paste("Gamma:", 100*length(arm.lengths)/sum(arm.lengths))) ## Gamma age from Gandolfo paper
    print(paste("Gamma:", 100/(2/tt.gamma$par[1]))) ## Gamma age
}