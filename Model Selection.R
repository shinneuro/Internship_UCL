# Analysis ----------------------------------------------------------------


env_models = new.env()
lapply("MinhoModels_DeltaSMf.RData", load, envir = env_models)
models = as.list(env_models)
model_name <- names(models) #sapply(seq_along(models), function(idx) names(models)[idx])

# models <- lapply("MinhoModels_DeltaSMf.Rdata",load,.GlobalEnv)

N <- length(models[[1]]) # number of participants
M <- length(models) # number of models

participant <- vector("list",length=length(models[[1]]))
para <- vector("list",length=5)
AIC <- vector("double",length=length(models))
delta_AIC <- vector("double",length=length(models))
w_AIC <- vector("double",length=length(models))
N_AIC <- vector("integer",length=length(models))
BIC <- vector("double",length=length(models))
delta_BIC <- vector("double",length=length(models))
w_BIC <- vector("double",length=length(models))
N_BIC <- vector("integer",length=length(models))


for (i in 1:N) { # number of participants
  for (j in 1:M) { # number of models
    para[[j]] <- models[[j]][[i]]$par
    AIC[j] <- models[[j]][[i]]$value + 2*length(para[[j]])
    BIC[j] <- models[[j]][[i]]$value + log(200)*length(para[[j]]) # number of choice
  }
  # delta
  for (j in 1:M) {
    delta_AIC[j] <- AIC[j]-min(AIC)
    delta_BIC[j] <- BIC[j]-min(BIC)
  }
  # weights
  for (j in 1:M) {
    w_AIC[j] <- exp(-0.5*delta_AIC[j]) / sum(exp(-0.5*delta_AIC))
    w_BIC[j] <- exp(-0.5*delta_BIC[j]) / sum(exp(-0.5*delta_BIC))
    
  }
  
  selected_AIC <- model_name[which.max(w_AIC)]
  selected_BIC <- model_name[which.max(w_BIC)]
  
  N_AIC[which.max(w_AIC)] <- N_AIC[which.max(w_AIC)] + 1
  N_BIC[which.max(w_BIC)] <- N_BIC[which.max(w_BIC)] + 1
  
  participant[[i]] <- list(para, AIC, delta_AIC, w_AIC, selected_AIC, BIC, delta_BIC, w_BIC, selected_BIC)
  names( participant[[i]]) <- c("para","AIC", "delta_AIC", "w_AIC", "selected_AIC", "BIC", "delta_BIC", "w_BIC", "selected_BIC")
}



save(participant, N_AIC, N_BIC, file="participants_analysis.RData")

