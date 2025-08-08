## Data simulation using MIDAsim (More info on https://github.com/mengyu-he/MIDASim)

# Load packages

library(MIDASim)
library(nleqslv)
library(survival)
library(microbiome)

# General setup
# Which settings to explore

alln <- c(50, 200) # total sample size
allb <- c(1.5, 3, 5)  # effect size
rnt_perc <- c(0.1,0.2) # number of causal taxa

eg <- expand.grid(n = alln, beta = allb, rnt_perc = rnt_perc, 
                  stringsAsFactors = F)

#------------DIETSWAP-------------
path2 = "path/to/save/simulated counts/for /Dietswap/data"

data("dietswap")
dietswap2<-subset_samples(dietswap, timepoint==2)
count.dietswap <- t(as.matrix(dietswap2@otu_table))
count.dietswap.setup <- MIDASim.setup(count.dietswap, mode = 'parametric')

dim(count.dietswap)
for(iscen in 1:nrow(eg)){
  iscen=6
  n.taxa = count.dietswap.setup$n.taxa

  condition <- eg[iscen,]
  n.data <- condition$n
  beta1 <- condition$beta # effect of the variable of interest
  rnt <- round(condition$rnt_perc * n.taxa)

  fn <- sprintf("Setting%d(n_%d_b_%s_rnt_%d)",iscen,n.data,beta1,rnt)

  dir.create(file.path(path2, fn), showWarnings = FALSE)
  
  
  # Top 60 of most abundant taxa are selected for possible modification
  m = 60
  top.index = order(count.dietswap.setup$mean.rel.abund, decreasing =T)[1:m]
  x.index = sample(top.index, rnt, replace =F) # get causal taxa set: randomly select rnt=10/20 taxa from the most abundant 50 taxa
  
  # Group indicator (added to the data in the end, i.e. defines the two groups)
  x.true = c(rep(0, n.data/2), rep(1, n.data/2))
  x.true = x.true - mean(x.true)
  
  # Generate 100 datasets per scenario under consideration
  
  data.all = ra.all = matrix(NA, n.data, n.taxa)
 
  new.mean.rel.abund = old.mean.rel.abund= count.dietswap.setup$mean.rel.abund
  
   fraction1 = 0.7
   fraction1_selected = x.index[sample(1:length(x.index),round(fraction1*length(x.index)),replace=FALSE)]
   fraction2_selected = x.index[!(x.index %in%fraction1_selected)]

   new.mean.rel.abund[fraction1_selected] =   old.mean.rel.abund[fraction1_selected] * 1/beta1
   new.mean.rel.abund[fraction2_selected] =   old.mean.rel.abund[fraction2_selected] * beta1
   new.mean.rel.abund <- new.mean.rel.abund/sum(new.mean.rel.abund)
  
  new.lib.size = sample(min(count.dietswap.setup$lib.size):max(count.dietswap.setup$lib.size), size=n.data/2, replace=TRUE)
  
  count.dietswap.modified.1 <- try(MIDASim.modify(count.dietswap.setup,
                                             mean.rel.abund = old.mean.rel.abund,
                                             lib.size = new.lib.size))
  
  count.dietswap.modified.2 <- try(MIDASim.modify(count.dietswap.setup,
                                             mean.rel.abund = new.mean.rel.abund,
                                             lib.size = new.lib.size))
  
  
  if ("try-error" %in% class(count.dietswap.modified.2)){
    break()
  }else{
    for(sim in 1:100){
      temp1 <- MIDASim(count.dietswap.modified.1)
      temp2 <- MIDASim(count.dietswap.modified.2)
      
      data.all <- rbind(temp1$sim_count,temp2$sim_count)
      ra.all<- rbind(temp1$sim_rel,temp2$sim_rel)
      
      if (any(is.na(rowSums(data.all)))){
        next()
      }
      
      cens.prob <- colMeans(data.all == 0) #the proportion of samples where the count of that taxon is zero.
      
      trueDA <- rep(0, n.taxa)
      trueDA[x.index] <- 1
      colnames(data.all) <- paste0("otu", 1:ncol(data.all))
      midas.data <- list(trueDA, data.all, cens.prob, x.true)
      fn2 <- file.path(file.path(path2, fn),sprintf("sim%d.Rdata",sim))
      save(midas.data, file = fn2)
    }}}





#-----------IBD-----------
path = "path/to/save/IBD/data"
data("count.ibd")
count.ibd.setup <- MIDASim.setup(count.ibd, mode = 'parametric')

for(iscen in 1:nrow(eg)){
  #prop_DA<-0.1
  n.taxa = count.ibd.setup$n.taxa
  
  condition <- eg[iscen,]
  n.data <- condition$n
  beta1 <- condition$beta # effect of the variable of interest
  rnt <- round(condition$rnt_perc * n.taxa)
  
  fn <- sprintf("Setting%d(n_%d_b_%s_rnt_%d)",iscen,n.data,beta1,rnt)
  dir.create(file.path(path3, fn), showWarnings = FALSE)
  
  
  # Top 60 of most abundant taxa are selected for possible modification
  m = 60
  top.index = order(count.ibd.setup$mean.rel.abund, decreasing =T)[1:m]
  x.index = sample(top.index, rnt, replace = T) # get causal taxa set: randomly select rnt=10/20 taxa from the most abundant 50 taxa
  
  # Group indicator (added to the data in the end, i.e. defines the two groups)
  x.true = c(rep(0, n.data/2), rep(1, n.data/2))
  x.true = x.true - mean(x.true)
  
  # Generate 100 datasets per scenario under consideration
  data.all = ra.all = matrix(NA, n.data, n.taxa)
  
  new.mean.rel.abund = old.mean.rel.abund= count.ibd.setup$mean.rel.abund
  
  #Compositional data (Upregulated and down regulated DA taxa)
  fraction1 = 0.7
  fraction1_selected = x.index[sample(1:length(x.index),round(fraction1*length(x.index)),replace=FALSE)]
  fraction2_selected = x.index[!(x.index %in%fraction1_selected)]

  new.mean.rel.abund[fraction1_selected] =   old.mean.rel.abund[fraction1_selected] * 1/beta1
  new.mean.rel.abund[fraction2_selected] =   old.mean.rel.abund[fraction2_selected] * beta1
 
  new.mean.rel.abund <- new.mean.rel.abund/sum(new.mean.rel.abund)
  
  new.lib.size = sample(min(count.ibd.setup$lib.size):max(count.ibd.setup$lib.size), size=n.data/2, replace=TRUE)
  
  count.ibd.modified.1 <- try(MIDASim.modify(count.ibd.setup,
                                                  mean.rel.abund = old.mean.rel.abund,
                                                  lib.size = new.lib.size))
  
  count.ibd.modified.2 <- try(MIDASim.modify(count.ibd.setup,
                                                  mean.rel.abund = new.mean.rel.abund,
                                                  lib.size = new.lib.size))
  
  
  if ("try-error" %in% class(count.ibd.modified.2)){
    break()
  }else{
    for(sim in 1:100){
      temp1 <- MIDASim(count.ibd.modified.1)
      temp2 <- MIDASim(count.ibd.modified.2)
      
      data.all <- rbind(temp1$sim_count,temp2$sim_count)
      ra.all<- rbind(temp1$sim_rel,temp2$sim_rel)
      
      if (any(is.na(rowSums(data.all)))){
        next()
      }
      
      cens.prob <- colMeans(data.all == 0) #the proportion of samples where the count of that taxon is zero.
      
      trueDA <- rep(0, n.taxa)
      trueDA[x.index] <- 1
      colnames(data.all) <- paste0("otu", 1:ncol(data.all))
      midas.data <- list(trueDA, data.all, cens.prob, x.true)
      fn2 <- file.path(file.path(path3, fn),sprintf("sim%d.Rdata",sim))
      save(midas.data, file = fn2)
    }}}

#------------DIETSWAP WITHOUT COMPOSITIONALITY----------
# Which settings to explore

alln <- c(50, 200) # total sample size
allb <- c(1, 1.5, 3, 5)  # effect size
rnt <- c(20,50) # number of causal taxa

eg <- expand.grid(n = alln, beta = allb, rnt = rnt, 
                  stringsAsFactors = F)

for(iscen in 1:nrow(eg)){
  
  condition <- eg[iscen,]
  n.data <- condition$n
  beta1 <- condition$beta # effect of the variable of interest
  rnt <- condition$rnt
  
  n.taxa = count.dietswap.setup$n.taxa
  
  fn <- sprintf("Setting%d(highcomp)",iscen)
  
  dir.create(file.path("./Sim_datasets/Dietswap", fn), showWarnings = FALSE)
  
  
  # Top 60 of most abundant taxa are selected for possible modification
  m = 60
  top.index = order(count.dietswap.setup$mean.rel.abund, decreasing =T)[1:m]
  x.index = sample(top.index, rnt, replace =F) # get causal taxa set: randomly select rnt=10/20 taxa from the most abundant 50 taxa
  
  # Group indicator
  x.true = c(rep(0, n.data/2), rep(1, n.data/2))
  x.true = x.true - mean(x.true)
  
  # Generate 100 datasets per scenario under consideration
  
  data.all = ra.all = matrix(NA, n.data, n.taxa)
  
  # effect1 <- (x.true * beta1)
  # new.mean.rel.abund = old.mean.rel.abund= count.dietswap.setup$mean.rel.abund
  # new.mean.rel.abund[x.index[1:(rnt/2)]] =   new.mean.rel.abund[x.index[1:(rnt/2)]] *exp(0.5*beta1)
  # new.mean.rel.abund[x.index[(rnt/2+1):(rnt)]] =  new.mean.rel.abund[x.index[(rnt/2+1):(rnt)]] *exp(-0.5*beta1)
  # 
  # new.mean.rel.abund <- new.mean.rel.abund/sum(new.mean.rel.abund)
  
  
  new.mean.rel.abund = old.mean.rel.abund= count.dietswap.setup$mean.rel.abund
  
  
  ## Following Thas (personal communication 7/2/2025)
  fraction1 = 0.9
  fraction1_selected = x.index[sample(1:length(x.index),round(fraction1*length(x.index)),replace=FALSE)]
  fraction2_selected = x.index[!(x.index %in%fraction1_selected)]
  
  new.mean.rel.abund[fraction1_selected] =   old.mean.rel.abund[fraction1_selected] * 1/beta1
  new.mean.rel.abund[fraction2_selected] =   old.mean.rel.abund[fraction2_selected] * 1.2
  
  sum0 = sum(old.mean.rel.abund[-x.index])
  sum1 = sum(new.mean.rel.abund[fraction1_selected])
  sum2 = sum(new.mean.rel.abund[fraction2_selected])
  table(new.mean.rel.abund/old.mean.rel.abund)
  gam_d = (1-sum0-sum1)/sum2
  new.mean.rel.abund[fraction2_selected] =   old.mean.rel.abund[fraction2_selected] *  1.2 *gam_d
  
  
  
  ## High compositionality setting (no adjustment, so all LFC of all taxa change due to compositionality, setting 16 used)
  
  #new.mean.rel.abund[x.index] =   new.mean.rel.abund[x.index] * beta1
  #new.mean.rel.abund <- new.mean.rel.abund/sum(new.mean.rel.abund)
  
  new.lib.size = sample(min(count.dietswap.setup$lib.size):max(count.dietswap.setup$lib.size), size=n.data/2, replace=TRUE)
  
  count.dietswap.modified.1 <- try(MIDASim.modify(count.dietswap.setup,
                                                  mean.rel.abund = old.mean.rel.abund,
                                                  lib.size = new.lib.size))
  
  count.dietswap.modified.2 <- try(MIDASim.modify(count.dietswap.setup,
                                                  mean.rel.abund = new.mean.rel.abund,
                                                  lib.size = new.lib.size))
  
  
  if ("try-error" %in% class(count.dietswap.modified.2)){
    break()
  }else{
    for(sim in 1:100){
      temp1 <- MIDASim(count.dietswap.modified.1)
      temp2 <- MIDASim(count.dietswap.modified.2)
      
      data.all <- rbind(temp1$sim_count,temp2$sim_count)
      ra.all<- rbind(temp1$sim_rel,temp2$sim_rel)
      
      if (any(is.na(rowSums(data.all)))){
        next()
      }
      
      cens.prob <- colMeans(data.all == 0) #the proportion of samples where the count of that taxon is zero.
      
      trueDA <- rep(0, n.taxa)
      trueDA[x.index] <- 1
      colnames(data.all) <- paste0("otu", 1:ncol(data.all))
      midas.data <- list(trueDA, data.all, cens.prob, x.true)
      fn2 <- file.path(file.path("./Sim_datasets/Dietswap", fn),sprintf("sim%d.Rdata",sim))
      save(midas.data, file = fn2)
    }}}

#-------IBD WITHOUT COMPOSITIONALITY-------------
## Data simulation using MIDAsim
data("count.ibd")
count.ibd.setup <- MIDASim.setup(count.ibd, mode = 'parametric')

# Which settings to explore

alln <- c(50, 200) # total sample size
allb <- c(1, 1.5, 3, 5)  # effect size
rnt <- c(20,50) # number of causal taxa

eg <- expand.grid(n = alln, beta = allb, rnt = rnt, 
                  stringsAsFactors = F)
iscen = 16
for(iscen in 1:nrow(eg)){
  
  condition <- eg[iscen,]
  n.data <- condition$n
  beta1 <- condition$beta # effect of the variable of interest
  rnt <- condition$rnt
  
  n.taxa = count.ibd.setup$n.taxa
  
  fn <- sprintf("Setting%d(highcomp)",iscen)
  dir.create(file.path("./Sim_datasets", fn), showWarnings = FALSE)
  
  
  # Top 60 of most abundant taxa are selected for possible modification
  m = 60
  top.index = order(count.ibd.setup$mean.rel.abund, decreasing =T)[1:m]
  x.index = sample(top.index, rnt, replace =F) # get causal taxa set: randomly select rnt=10/20 taxa from the most abundant 50 taxa
  
  # Group indicator
  x.true = c(rep(0, n.data/2), rep(1, n.data/2))
  x.true = x.true - mean(x.true)
  
  # Generate 100 datasets per scenario under consideration
  
  data.all = ra.all = matrix(NA, n.data, n.taxa)
  
  effect1 <- (x.true * beta1)
  new.mean.rel.abund = old.mean.rel.abund= count.ibd.setup$mean.rel.abund
  new.mean.rel.abund[x.index[1:(rnt/2)]] =   new.mean.rel.abund[x.index[1:(rnt/2)]] *exp(0.5*beta1)
  new.mean.rel.abund[x.index[(rnt/2+1):(rnt)]] =  new.mean.rel.abund[x.index[(rnt/2+1):(rnt)]] *exp(-0.5*beta1)
  
  new.mean.rel.abund <- new.mean.rel.abund/sum(new.mean.rel.abund)
  
  
  new.mean.rel.abund = old.mean.rel.abund= count.ibd.setup$mean.rel.abund
  
  
  # Following Thas (personal communication 7/2/2025) -> used for all settings to make sure only the selected taxa differ
  fraction1 = 0.9
  fraction1_selected = x.index[sample(1:length(x.index),round(fraction1*length(x.index)),replace=FALSE)]
  fraction2_selected = x.index[!(x.index %in%fraction1_selected)]
  
  new.mean.rel.abund[fraction1_selected] =   old.mean.rel.abund[fraction1_selected] * 1/beta1
  new.mean.rel.abund[fraction2_selected] =   old.mean.rel.abund[fraction2_selected] * 1.2
  
  sum0 = sum(old.mean.rel.abund[-x.index])
  sum1 = sum(new.mean.rel.abund[fraction1_selected])
  sum2 = sum(new.mean.rel.abund[fraction2_selected])
  
  gam_d = (1-sum0-sum1)/sum2
  new.mean.rel.abund[fraction2_selected] =   old.mean.rel.abund[fraction2_selected] *  1.2 *gam_d
  
  
  
  new.lib.size = sample(min(count.ibd.setup$lib.size):max(count.ibd.setup$lib.size), size=n.data/2, replace=TRUE)
  
  count.ibd.modified.1 <- try(MIDASim.modify(count.ibd.setup,
                                             mean.rel.abund = old.mean.rel.abund,
                                             lib.size = new.lib.size))
  
  count.ibd.modified.2 <- try(MIDASim.modify(count.ibd.setup,
                                             mean.rel.abund = new.mean.rel.abund,
                                             lib.size = new.lib.size))
  
  
  if ("try-error" %in% class(count.ibd.modified.2)){
    break()
  }else{
    for(sim in 1:100){
      temp1 <- MIDASim(count.ibd.modified.1)
      temp2 <- MIDASim(count.ibd.modified.2)
      
      data.all <- rbind(temp1$sim_count,temp2$sim_count)
      ra.all<- rbind(temp1$sim_rel,temp2$sim_rel)
      
      if (any(is.na(rowSums(data.all)))){
        next()
      }
      
      cens.prob <- colMeans(data.all == 0) #the proportion of samples where the count of that taxon is zero.
      
      trueDA <- rep(0, n.taxa)
      trueDA[x.index] <- 1
      colnames(data.all) <- paste0("otu", 1:ncol(data.all))
      midas.data <- list(trueDA, data.all, cens.prob, x.true)
      fn2 <- file.path(file.path("./Sim_datasets", fn),sprintf("sim%d.Rdata",sim))
      save(midas.data, file = fn2)
    }}}


#--------Creating phyloseq objects from the simulated counts------
create_phyloseq <- function(midas_data) {
  # otu table
  midas_otutable <- midas_data[[2]]
  rownames(midas_otutable) <- paste0("sample", 1:nrow(midas_otutable))  # Set row names for samples
  
  # taxa table
  n.taxa <- ncol(midas_otutable)
  midas_taxtable <- matrix(nrow = n.taxa, ncol = 2)
  midas_taxtable[, 1:2] <- c(midas_data[[1]], midas_data[[3]])
  colnames(midas_taxtable) <- c("DA_ind", "cens.prob")
  midas_taxtable <- as.data.frame(midas_taxtable) 
  midas_taxtable$isDA <- ifelse(midas_taxtable$DA_ind == 1, "TRUE", "FALSE")
  rownames(midas_taxtable) <- colnames(midas_otutable)
  midas_taxtable <- as.matrix(midas_taxtable)
  
  # sample data
  midas_samdata <- as.data.frame(midas_data[[4]])
  rownames(midas_samdata) <- rownames(midas_otutable)
  colnames(midas_samdata) <- "x.index"
  midas_samdata$group <- ifelse(midas_samdata$x.index == 0.5, 1, 0)  # Create a 'group' column for the groupings
  
  # Create phyloseq object
  midas_phyloseq <- phyloseq(otu_table(midas_otutable, taxa_are_rows = FALSE),
                             sample_data(midas_samdata),
                             tax_table(midas_taxtable))
  
  return(midas_phyloseq)
}

#### Put datasets in a list ####
# Transforming the counts into a phyloseq object
#### Put the datasets in each setting in a list ####
data_files <- list.files(path2 = path,pattern = "*.Rdata", full.names=TRUE)
#per setting: Example setting 3
Setting6_n_200_b_5_rnt_13 <- list()
for (i in seq_along(data_files)) {
  temp_env <- new.env()
  load(data_files[i], envir = temp_env)
  Setting6_n_200_b_5_rnt_13[[i]] <- temp_env[[ls(temp_env)[1]]]
}

# Loop through the list and create phyloseq object for each dataset
Setting6_n200_b5_rnt13 <- list()
for (i in seq_along(Setting6_n_200_b_5_rnt_13)) {
  Setting6_n200_b5_rnt13[[i]] <- create_phyloseq(Setting6_n_200_b_5_rnt_13[[i]])
}

setwd("folder to save the data lists")
save(Setting6_n200_b5_rnt13, file = 'Setting6_n200_b5_rnt13.RData')



