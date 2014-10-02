jagsparallel <- function(data, inits, parameters.to.save,
                         model.file = "model.bug", n.chains = 2, n.iter = 2000, 
                         n.burnin = floor(n.iter/2),
                         n.thin = max(1, floor((n.iter - n.burnin)/1000)), 
                         n.cluster = n.chains, DIC = TRUE, working.directory = 
                         NULL, jags.seed = 123, digits = 5,
                         RNGname = c("Wichmann-Hill", "Marsaglia-Multicarry", 
                                     "Super-Duper", "Mersenne-Twister"), 
                         jags.module = c("glm", "dic")) {
    library(R2jags)
    library(foreach)
    library(doParallel)
    library(doRNG)
    library(abind)
    jags.inits <- if (missing(inits)) {
        NULL
    } else {
        inits
    }
    cl <- makeCluster(n.cluster)
    registerDoParallel(cl)
    jags.model <- model.file
    set.seed(jags.seed)
    res <- foreach (n=1:n.chains, .inorder=FALSE, .packages=c("R2jags")) %dorng% {
        jags(data, inits, parameters.to.save=parameters.to.save, 
             model.file=model.file, n.chains=1, n.iter = n.iter,
             n.burnin = n.burnin, n.thin = n.thin, DIC = DIC,
             working.directory = working.directory, jags.seed = jags.seed, 
             progress.bar = "none", digits = digits, RNGname = RNGname, 
             jags.module = jags.module)
    }
    stopCluster(cl)
    result <- NULL
    model <- NULL
    for (ch in 1:n.chains) {
        result <- abind(result, res[[ch]]$BUGSoutput$sims.array, along = 2)
        model[[ch]] <- res[[ch]]$model
    }
    if (is.function(model.file)) {
        model.file <- substitute(model.file)
    }
    result <- as.bugs.array2(result, model.file = model.file, program = "jags", 
                             DIC = DIC, n.iter = n.iter, n.burnin = n.burnin, 
                             n.thin = n.thin)
    out <- list(model = model, BUGSoutput = result, parameters.to.save = 
                parameters.to.save, model.file = model.file, n.iter = n.iter, 
                DIC = DIC)
    class(out) <- c("rjags.parallel", "rjags")
    return(out)
}
