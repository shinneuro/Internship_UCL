# Analysis ----------------------------------------------------------------


env_models = new.env()
lapply("final_prospect_only.RData", load, envir = env_models)
models = as.list(env_models)
model_name <- names(models) #sapply(seq_along(models), function(idx) names(models)[idx])

# models <- lapply("MinhoModels_DeltaSMf.Rdata",load,.GlobalEnv)

N <- length(models[[1]]) # number of participants
M <- length(models) # number of models

participant <- vector("list",length=length(models[[1]]))
para <- vector("list",length=length(models))
names(para) <- model_name
AIC <- vector("double",length=length(models))
delta_AIC <- vector("double",length=length(models))
w_AIC <- vector("double",length=length(models))
#N_AIC <- matrix(0,ncol=length(models), nrow = 1)
#colnames(N_AIC) <- model_name
BIC <- vector("double",length=length(models))
delta_BIC <- vector("double",length=length(models))
w_BIC <- vector("double",length=length(models))
#N_BIC <- matrix(0,ncol=length(models), nrow = 1)
#colnames(N_BIC) <- model_name
N_IC <-  matrix(0,ncol=length(models), nrow = 2)
rownames(N_IC) <- c("N_AIC","N_BIC")
colnames(N_IC) <- model_name


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
  
  N_IC[1,which.max(w_AIC)] <- N_IC[1,which.max(w_AIC)] + 1 # N_AIC
  N_IC[2,which.max(w_BIC)] <- N_IC[2,which.max(w_BIC)] + 1 # N_BIC
  
  participant[[i]] <- list(para, AIC, delta_AIC, w_AIC, selected_AIC, BIC, delta_BIC, w_BIC, selected_BIC)
  names( participant[[i]]) <- c("para","AIC", "delta_AIC", "w_AIC", "selected_AIC", "BIC", "delta_BIC", "w_BIC", "selected_BIC")
}


selected <- matrix(0,nrow=79,ncol=2)
colnames(selected) <- c("selected_AIC","selected_BIC")

for (i in 1:79) {
   selected[i,1] <- participant[[i]]$selected_AIC
   selected[i,2] <- participant[[i]]$selected_BIC
}


save(participant, N_IC, selected, file="participants_analysis_final_prospect_only.RData")

