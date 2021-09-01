## Repeatability in avian migratory timings

#devtools::install_github("itchyshin/orchard_plot", subdir = "orchaRd", force = TRUE, build_vignettes = TRUE)

# Load packages #######################################
pacman::p_load(SciViews,
               tidyverse,
               metafor, # package for meta-analysis
               rotl, # package for phylogeny
               ape, # package for phylogeny
               cowplot, # combining multiple plots
               here, # making folder path usable for all
               clubSandwich, # package to assist metafor
               orchaRd,
               MuMIn,
               stringr
)

# Functions needed #######################################

# Function for Fisher's Z transformation (Zr) for correlation-based repeatabilities (r) and ICC (Holtmann et al. 2017, Table 1)
Zr_transformation <- function(r,K,Est) {
  
  if(Est == "ICC") {Zr <- 0.5*ln(((1+(K-1)*r)/(1-r)))}
  if(Est == "r") {Zr <- 0.5*ln(((1+r)/(1-r)))} # TODO I fixed this- but I guess for r k is 2 so probably does not matter
  Zr
}

# Function for sampling variance for correlation-based repeatabilities (r) and ICC (Holtmann et al. 2017, Table 1)
Calc_SV <- function(K,N,Est){
  
  if(Est == "ICC") {VZr <- K/(2*((N-2)*(K-1)))}
  if(Est == "r") {VZr <- 1/(N-3)}
  VZr
}

# Function to obtain I^2 total and separate I2 from multilevel-meta-analytic model (written by Shinichi Nakagawa)
I2 <- function(model, method = c("Wolfgang", "Shinichi")){
  ## evaluate choices
  method <- match.arg(method)
  
  # Wolfgang's method
  if(method == "Wolfgang"){
    W <- solve(model$V) 
    X <- model.matrix(model)
    P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
    I2_total <- sum(model$sigma2) / (sum(model$sigma2) + (model$k - model$p) / sum(diag(P)))
    I2_each  <- model$sigma2 / (sum(model$sigma2) + (model$k - model$p) / sum(diag(P)))
    names(I2_each) = paste0("I2_", model$s.names)
    
    # putting all together
    I2s <- c(I2_total = I2_total, I2_each)
    
    # or my way
  } else {
    # sigma2_v = typical sampling error variance
    sigma2_v <- sum(1/model$vi) * (model$k-1) / (sum(1/model$vi)^2 - sum((1/model$vi)^2)) 
    I2_total <- sum(model$sigma2) / (sum(model$sigma2) + sigma2_v) #s^2_t = total variance
    I2_each  <- model$sigma2 / (sum(model$sigma2) + sigma2_v)
    names(I2_each) = paste0("I2_", model$s.names)
    
    # putting all together
    I2s <- c(I2_total = I2_total, I2_each)
  }
  return(I2s)
}

# Function to get estimates from rma objects (metafor) (Hayward et al. 2021)
get_est <- function(model, mod = " ") {
  
  name <- as.factor(str_replace(row.names(model$beta), mod, ""))
  estimate <- as.numeric(model$beta)
  lowerCL <- model$ci.lb
  upperCL <- model$ci.ub
  
  table <- tibble(name = name, estimate = estimate, lowerCL = lowerCL, upperCL = upperCL)
}

# Function function to get prediction intervals (crediblity intervals) from rma objects (metafor) (Hayward et al. 2021)
get_pred <- function(model, mod = " ") {
  name <- as.factor(str_replace(row.names(model$beta), mod, ""))
  len <- length(name)
  
  if (len != 1) {
    newdata <- matrix(NA, ncol = len, nrow = len)
    for (i in 1:len) {
      # getting the position of unique case from X (design matrix)
      pos <- which(model$X[, i] == 1)[[1]]
      newdata[, i] <- model$X[pos, ]
    }
    pred <- predict.rma(model, newmods = newdata)
  } else {
    pred <- predict.rma(model)
  }
  lowerPR <- pred$cr.lb
  upperPR <- pred$cr.ub
  
  table <- tibble(name = name, lowerPR = lowerPR, upperPR = upperPR)
}


R2 <- function(model){
  warning("Conditional R2 is not meaningful and the same as marginal R2\n")
  
  # fixed effect variance
  fix <- var(as.numeric(as.vector(model$b) %*% t(as.matrix(model$X))))
  
  # marginal
  R2m <- fix / (fix + sum(model$sigma2))
  R2
  #Rm <- round(100*R2m, 3)
  
  # conditional
  R2c <- (fix + sum(model$sigma2) - model$sigma2[length(model$sigma2)]) / 
    (fix + sum(model$sigma2))
  
  R2s <- c(R2_marginal = R2m, R2_coditional = R2c)
  return(R2s)
}


# custom function for extracting mean and CI from each metafor model (Publication bias paper)
estimates.CI <- function(model){
  db.mf <- data.frame(model$b,row.names = 1:nrow(model$b))
  db.mf <- cbind(db.mf,model$ci.lb,model$ci.ub,row.names(model$b))
  names(db.mf) <- c("mean","lower","upper","estimate")
  return(db.mf[,c("estimate","mean","lower","upper")])
}


# Read in and sort data #######################################

df <- read_csv(here("data", "Meta-analysis_data.csv")) # I had to use read_csv ratther than read.csv

# Descriptive 
df %>% summarise(mean = mean(n), median = median(n), min = min(n), max = max(n), mean_k = mean(k), median_k = median(k), min_k = min(k), max_k = max(k))
df %>% group_by(method) %>% summarise(median = median(n), mean = mean(n), min = min(n), max = max(n))
df %>% group_by(method, taxa) %>% summarise(median = median(n), mean = mean(n), min = min(n), max = max(n), cohort= n_distinct(cohort_ID), paper=n_distinct(paper_ID), es=n_distinct(es_ID))
df %>% group_by(method, taxa) %>% summarise(median = median(k), mean = mean(k), min = min(k), max = max(k))
df %>% group_by(taxa) %>% summarise(n_distinct(species_ID))
no.spp <- df %>% group_by(species_ID) %>% summarise(n_distinct(paper_ID))

# Calculate effect sizes and their sampling variances #######################################
# Creating new columns for effect sizes (Zr) and sampling variance (VZr)
df$VZr <- df$Zr <- NA

# Calculate effect sizes using function above
for (i in as.numeric(rownames(df))) {
  r <- df$R[i]
  K <- df$k[i]
  Est <- df$est[i]
  
  Zr <- Zr_transformation(r, K, Est)
  
  df$Zr[i] <- Zr
  
}
# Calculate sampling variances using function above
for (i in as.numeric(rownames(df))) {
  K <- df$k[i]
  N <- df$n[i]
  Est <- df$est[i]
  
  VZr <- Calc_SV(K,N,Est)
  
  df$VZr[i] <- VZr
  
}



# Phylogeny #######################################

# Using Jetz data (Holtmann et al. 2017) #######################################
species_list <- unique(df$species_latin) # use unique() as some names are repeated
species_list<-gsub(" ", "_", species_list) # replace spaces with underscore
str(species_list)

# Load bird supertree based on Hackett's backbone 
bird_tree_hackett <- read.tree("data/Hackett.tre") # tree provided by Benedikt Holtmann

# Prune phylogentic tree for meta-analysis

# Check the tree
bird_tree_hackett # 9993 tips = species
str(bird_tree_hackett) # has edge (branch) lengths
bird_tree_hackett <- collapse.singles(bird_tree_hackett)

# Get only bird species from the supertree that are also included in collected data
# Extract list of all species from the supertree
bird_tree_species <- as.character(bird_tree_hackett$tip.label)

# Check the overlap of species names between collected data file and the supertree
# All species should be present. If not, species names may not match with names in the supertree
intersect(bird_tree_species, species_list) 
# character(39) - should be 47, need to check which species do/do not match

# gives list of species names which are not matched
species_list[!species_list %in% bird_tree_species]
#"Limosa_lapponica_baueri"
#"Cygnus_columbianus_bewickii"
#"Limosa_limosa_limosa"
#"Limosa_limosa_islandica"
#"Catharacta_antarctica_lonnbergi"
#"Pterodroma_deserta"
#"Chen_canagicus"
#"Anser_caerulescens_atlanticus"

species_list_hackett <- species_list

# Changing subspecies to species, or old genus names from papers to new 
# Only pterodroma deserta: 3 species in the Pterodroma feae/madeira/desertae complex were once believed to be subspecies of a single species: Pterodroma mollis
species_list_hackett <- replace(species_list_hackett, species_list_hackett=="Limosa_lapponica_baueri", "Limosa_lapponica")
species_list_hackett <- replace(species_list_hackett, species_list_hackett=="Limosa_limosa_limosa" | species_list_hackett=="Limosa_limosa_islandica" , "Limosa_limosa")
species_list_hackett <- replace(species_list_hackett, species_list_hackett=="Cygnus_columbianus_bewickii", "Cygnus_columbianus")
species_list_hackett <- replace(species_list_hackett, species_list_hackett=="Catharacta_antarctica_lonnbergi", "Catharacta_antarctica")
species_list_hackett <- replace(species_list_hackett, species_list_hackett=="Pterodroma_deserta", "Pterodroma_mollis")
species_list_hackett <- replace(species_list_hackett, species_list_hackett=="Chen_canagicus", "Chen_canagica")
species_list_hackett <- replace(species_list_hackett, species_list_hackett=="Anser_caerulescens_atlanticus", "Chen_caerulescens")

# Now check and see if all species are present in supertree
intersect(bird_tree_species, species_list_hackett) # = 46 (not 47, because L. l. limosa and L. i. islandica both changed to L. limosa)

# At the moment, 'species_ID' column is based on species rather than subspecies, so goes up to s046 not s047

# Prune supertree to the list of taxa included in the data
pruned_birds_stree <- drop.tip(bird_tree_hackett, bird_tree_hackett$tip.label[-match(species_list_hackett, bird_tree_hackett$tip.label)])

# Check if tree is binary and ultrametric
is.binary(pruned_birds_stree) # TRUE
is.ultrametric(pruned_birds_stree) # TRUE

# Save pruned tree to be use in the meta-analysis
write.tree(pruned_birds_stree,
           file = "data/birds_meta-analysis_tree.tre", append = FALSE,
           digits = 10, tree.names = FALSE
)

# Plot tree to see how it looks like
plot(pruned_birds_stree, label.offset = 2, cex = 0.8, main = "'Hackett tree'", cex.main = 1, line = 0.5) # with branch lengths

# Make phylogenetic correlation matrix

# Using metafor - correlation matrix for species in tree
varcor <- vcv(pruned_birds_stree, corr = TRUE)


# Using rotl package (O'Dea et al. 2019) #######################################
# matching names from open tree taxonomy
taxa <- tnrs_match_names(names = levels(df$species_latin), context_name = "Animals")
# which names return more than 1 match?
inspect(taxa, ott_id = taxa$ott_id[taxa$number_matches != 1])
# fixing names with more than 1 match
taxa[taxa$number_matches != 1, ] <- inspect(taxa, ott_id = taxa$ott_id[taxa$number_matches != 1])[2, ]
# Return the induced subtree on the synthetic tree that relates a list of nodes
tr <- tol_induced_subtree(ott_id(taxa), label="name")
plot(tr)

is.binary(tr) # TRUE


# Compute branch lengths  

# correlation matrix to fit to the model
tr$tip.label <- as.factor(tr$tip.label) 
levels(tr$tip.label) <- levels(df$species_latin) # making sure names match
tr$tip.label <- as.character(tr$tip.label) # converting names back to character

# compute branch lengths of tree
phylo_branch <- compute.brlen(tr, method = "Grafen", power = 1)

# check tree is ultrametric
is.ultrametric(phylo_branch) # TRUE

# saving phylogeneic matrix
phylo_cor <- vcv(phylo_branch, varcor = T)



# Meta-analysis #######################################
# Make sure to have calculated effect sizes and sampling variances above

# TODO - need to create a variance-covariance matrix at the cohort level (there is an alternative for this method)
VCV <- impute_covariance_matrix(vi = df$VZr, cluster = df$cohort_ID, r = 0.5)
# put this in for V = 
# Why 0.5 - see https://onlinelibrary.wiley.com/doi/full/10.1111/mec.14031

# TODO - you need to create a column called phylogeny, which matches your tree
# copied column to edit species names to match tree (as some sp. names are from old papers, or are subspecies)
df$species_latin_hackett <- df$species_latin 
df$species_latin_hackett <-replace(df$species_latin_hackett, df$species_latin_hackett == "Limosa lapponica baueri", "Limosa lapponica")
df$species_latin_hackett <- replace(df$species_latin_hackett, df$species_latin_hackett=="Limosa limosa limosa" | df$species_latin_hackett=="Limosa limosa islandica" , "Limosa limosa")
df$species_latin_hackett <-replace(df$species_latin_hackett, df$species_latin_hackett == "Cygnus columbianus bewickii", "Cygnus columbianus")
df$species_latin_hackett <-replace(df$species_latin_hackett, df$species_latin_hackett == "Catharacta antarctica lonnbergi", "Catharacta antarctica")
df$species_latin_hackett <-replace(df$species_latin_hackett, df$species_latin_hackett == "Pterodroma deserta", "Pterodroma mollis")
df$species_latin_hackett <-replace(df$species_latin_hackett, df$species_latin_hackett == "Chen canagicus", "Chen canagica")
df$species_latin_hackett <-replace(df$species_latin_hackett, df$species_latin_hackett == "Anser caerulescens atlanticus", "Chen caerulescens")


df$phylogeny<-gsub(" ", "_", df$species_latin_hackett)

#meta-analytic modelling
ma_test1 <- rma.mv(yi = Zr, V = VCV, # instead of VZr 
                   random = list(~1 | es_ID, 
                                 ~1 | paper_ID, 
                                 ~1 | cohort_ID, 
                                 ~1 | species_ID), 
                   data = df)

summary(ma_test1)

## Calculating I^2
round(i2_ml(ma_test1)*100,1) 

# TODO - see this paper - for reason why we need to put both species_ID and 
#
meta_model1 <- rma.mv(yi = Zr, V = VCV, 
                      random = list(~1 | es_ID, 
                                    ~1 | paper_ID, 
                                    ~1 | cohort_ID, 
                                    ~1 | species_ID, # non-phylo effect
                                    ~1 | phylogeny), # phylo effect
  R = list(phylogeny = varcor), # phylogenetic relatedness
  data = df)

summary(meta_model1)
tanh(0.5186); tanh(0.0687) # effect size estimate into r
# mean k 2.69
## Calculating I^2
round(i2_ml(meta_model1)*100,1) 

# I2_total      I2_es_ID   I2_paper_ID  I2_cohort_ID I2_species_ID  I2_phylogeny 
# 85.9          51.1           0.0           0.0          27.3           7.5 

# TODO you could keep it all but you could take cohort ID  + paper 
# or accoding to % explaning - we can take out a couple of random effects
# TODO - I forgot about -  method="ML")
meta_model2 <- rma.mv(yi = Zr, V = VCV, 
                      random = list(~1 | es_ID, 
                                    ~1 | paper_ID, 
                                   # ~1 | cohort_ID, 
                                    ~1 | species_ID, 
                                    ~1 | phylogeny
                                   ),
                      R = list(phylogeny = varcor), # added in phylogeny
                      data = df,
                      method="ML")

meta_model3 <- rma.mv(yi = Zr, V = VCV, 
                      random = list(~1 | es_ID, 
                                    #~1 | paper_ID, 
                                    # ~1 | cohort_ID, 
                                    ~1 | species_ID, 
                                    ~1 | phylogeny
                                    ),
                      R = list(phylogeny = varcor), # added in phylogeny
                      data = df, 
                      method="ML")

# we cannot take more random effects
# to do multilevel models we need at least one higher level
meta_model4 <- rma.mv(yi = Zr, V = VCV, 
                      random = list(~1 | es_ID, 
                                    #~1 | paper_ID, 
                                    # ~1 | cohort_ID, 
                                    ~1 | species_ID 
                                    #~1 | phylogeny
                                    ),
                      R = list(phylogeny = varcor), # added in phylogney
                      data = df,
                      method="ML")

# TODO 
# comparing AIC is one way or you can just delete ones which does not account very much or you can leave everything
aic1 <- AIC(meta_model1)
aic2 <- AIC(meta_model2) 
aic3 <- AIC(meta_model3) 
aic4 <- AIC(meta_model4) 
 
aic1; aic2; aic3; aic4

# or you can do liklihood ratio test
# these are all not all different s you could use model 4
# but you may want ot use model 1 - that his fine
# TODO, we have done a simulation study - here - https://ecoevorxiv.org/su4zv/ (almost accepted in Methods in Ecol Evol - minor revision)
anova(meta_model1, meta_model2)
anova(meta_model2, meta_model3)
anova(meta_model3, meta_model4)

# TODO - visualise with orchaRd plot

orchard_plot(meta_model1, xlab = "Correlation coefficient (r)", transfm = "tanh", cb="TRUE")



# TODO - please do meta-regression
# Method of tracking

meta_regression1 <- rma.mv(yi = Zr, V = VCV, 
                      mods = ~ method,
                      random = list(~1 | es_ID, 
                                    ~1 | paper_ID, 
                                    ~1 | cohort_ID, 
                                    ~1 | species_ID, 
                                    ~1 | phylogeny),
                      R = list(phylogeny = varcor), # added in phylogney
                      data = df)

summary(meta_regression1)

# reordering
df$method <- factor(df$method, levels = c("Conventional", "GLS", "Satellite"))

meta_regression1b <- rma.mv(yi = Zr, V = VCV, 
                            mods = ~ relevel(method, ref = "GLS"),
                            random = list(~1 | es_ID, 
                                          ~1 | paper_ID, 
                                          ~1 | cohort_ID, 
                                          ~1 | species_ID, 
                                          ~1 | phylogeny),
                            R = list(phylogeny = varcor), # added in phylogney
                            data = df)

summary(meta_regression1b)

r2_ml(meta_regression1)

res_meta_regression1b <- get_est(meta_regression1b, mod = "method")

tanh(res_meta_regression1b$estimate)

# TODO for orchard plot - need meta-regression without intercept (see vignette)
meta_regression1c <- rma.mv(yi = Zr, V = VCV, 
                            mods = ~ method -1,
                            random = list(~1 | es_ID, 
                                          ~1 | paper_ID, 
                                          ~1 | cohort_ID, 
                                          ~1 | species_ID, 
                                          ~1 | phylogeny),
                            R = list(phylogeny = varcor), # added in phylogney
                            data = df)

# orchard_plot(meta_regression1b, mod = "method", xlab = "Zr (effect size)")
# added transformation 'tanh' to convert back to the raw correlation coefficient scale (r) - but this is different to ICC ?
# e.g. tanh(0.4838) = 0.4993 (r) but ICC = 0.378 (which takes into account mean k of dataset, according to equations in Holtmann et al. 2017)
orchard_plot(meta_regression1c, mod = "method", xlab = "Correlation coefficient (r)", transfm = "tanh", cb="TRUE")


# orchard plot for method but with ICC ####
# function
to_ICC <- function(x,k){ 
  (exp(2*x)-1)/(exp(2*x)+k-1)
}

# get model results

mr1 <- mod_results(meta_regression1c, mod="method")
mr1_mod_table <- mr1$mod_table

mr1_data <- mr1$data

# back-transform model results to ICC
# calculate k for each method separately

df %>% group_by(method) %>% summarise(mean(k))


mr1_data <- mr1_data %>% mutate(k = case_when(moderator == "GLS" ~ 2.20, 
                                              moderator == "Conventional" ~ 3.13,
                                              moderator == "Satellite" ~ 3.28))

mr1_data$yi_ICC <- to_ICC(mr1_data$yi, mr1_data$k)

mr1_mod_table$k <- c(3.13, 2.20, 3.28)

for(i in names(mr1_mod_table)[2:6]){
  
  mr1_mod_table[i] <- to_ICC(mr1_mod_table[i], mr1_mod_table$k)
  
}


mr1_data$moderator <- factor(mr1_data$moderator, levels = mr1_mod_table$name, labels = mr1_mod_table$name)
mr1_data$scale <- (1/sqrt(mr1_data[,"vi"]))
mr1_mod_table$K <- as.vector(by(mr1_data, mr1_data[,"moderator"], function(x) length(x[,"yi"])))
mr1_group_no <- nrow(mr1_mod_table)


(fig_ma <- ggplot(data = mr1_mod_table, aes(x = estimate, y = name)) + 
    #scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.2)) + 
    ggbeeswarm::geom_quasirandom(data = mr1_data, aes(x = yi_ICC, y = moderator, size = scale, colour=moderator), groupOnX=FALSE, alpha = 0.5) +
    geom_errorbarh(aes(xmin = lowerPR, xmax = upperPR), height = 0, show.legend = F, size = 0.5, alpha = 0.6) + # CI
    geom_errorbarh(aes(xmin = lowerCL, xmax = upperCL), height = 0, show.legend = F, size = 1.2) + 
    geom_vline(xintercept = 0, linetype = 2, colour = "black", alpha = 0.5) + # creating dots and different size (bee-swarm and bubbles)
    geom_point(aes(fill=name), size = 3, shape = 21) + 
    ggplot2::annotate("text", x = max(mr1_data$yi_ICC) + (max(mr1_data$yi_ICC)*0.10), y = (seq(1, mr1_group_no, 1)+0.3), 
                      label = paste("italic(k)==", mr1_mod_table$K), parse = TRUE, hjust = "right", size = 3.5) +
    ggplot2::theme_bw() +
    ggplot2::guides(fill="none", colour="none") +
    ggplot2::theme(legend.position = c(1,0), legend.justification = c(1,0), legend.title = element_text(size=9), 
                   legend.direction = "horizontal", legend.background = element_blank(), axis.text.y = element_text(size=10, 
                                                                                                                    colour="black", hjust=0.5, angle=90)) +
    ggplot2::labs(x=paste("Effect size (ICC)"), y="", size=paste("Precision (1/SE)")))


# Annual event

# reordering
df$annual_event <- factor(df$annual_event, levels = c("Arrival_breed", "Depart_breed", "Nonbreed_arrival", "Nonbreed_depart"))

meta_regression2 <- rma.mv(yi = Zr, V = VCV, 
                           mods = ~ annual_event,
                           random = list(~1 | es_ID, 
                                         ~1 | paper_ID, 
                                         ~1 | cohort_ID, 
                                         ~1 | species_ID, 
                                         ~1 | phylogeny),
                           R = list(phylogeny = varcor), # added in phylogney
                           data = df)

summary(meta_regression2)


meta_regression2b <- rma.mv(yi = Zr, V = VCV, 
                            mods = ~ relevel(annual_event, ref = "Depart_breed"),
                            random = list(~1 | es_ID, 
                                          ~1 | paper_ID, 
                                          ~1 | cohort_ID, 
                                          ~1 | species_ID, 
                                          ~1 | phylogeny),
                            R = list(phylogeny = varcor), # added in phylogney
                            data = df)

summary(meta_regression2b)

meta_regression2c <- rma.mv(yi = Zr, V = VCV, 
                            mods = ~ relevel(annual_event, ref = "Nonbreed_arrival"),
                            random = list(~1 | es_ID, 
                                          ~1 | paper_ID, 
                                          ~1 | cohort_ID, 
                                          ~1 | species_ID, 
                                          ~1 | phylogeny),
                            R = list(phylogeny = varcor), # added in phylogney
                            data = df)

summary(meta_regression2c)

r2_ml(meta_regression2)

# TODO for orchard plot - need meta-regression without intercept (see vignette)
meta_regression2d <- rma.mv(yi = Zr, V = VCV, 
                            mods = ~ annual_event -1,
                            random = list(~1 | es_ID, 
                                          ~1 | paper_ID, 
                                          ~1 | cohort_ID, 
                                          ~1 | species_ID, 
                                          ~1 | phylogeny),
                            R = list(phylogeny = varcor), # added in phylogney
                            data = df)

# orchard_plot(meta_regression2b, mod = "annual_event", xlab = "Zr (effect size)")
orchard_plot(meta_regression2d, mod = "annual_event", xlab = "Correlation coefficient (r)", transfm = "tanh", cb="TRUE")

# Ecological group
# Does controlling for ecological group when having species_ID and phylogeny in as random effects make sense?

meta_regression3 <- rma.mv(yi = Zr, V = VCV, 
                           mods = ~ taxa,
                           random = list(~1 | es_ID, 
                                         ~1 | paper_ID, 
                                         ~1 | cohort_ID, 
                                         ~1 | species_ID, 
                                         ~1 | phylogeny),
                           R = list(phylogeny = varcor), # added in phylogney
                           data = df)

summary(meta_regression3)

r2_ml(meta_regression3)

# reordering
df$taxa <- factor(df$taxa, levels = c("Waterbird", "Seabird", "Landbird"))

meta_regression3b <- rma.mv(yi = Zr, V = VCV, 
                           mods = ~ relevel(taxa, ref = "Seabird"),
                           random = list(~1 | es_ID, 
                                         ~1 | paper_ID, 
                                         ~1 | cohort_ID, 
                                         ~1 | species_ID, 
                                         ~1 | phylogeny),
                           R = list(phylogeny = varcor), # added in phylogney
                           data = df)

summary(meta_regression3b)

# TODO for orchard plot - need meta-regression without intercept (see vignette)
meta_regression3c <- rma.mv(yi = Zr, V = VCV, 
                            mods = ~ taxa -1,
                            random = list(~1 | es_ID, 
                                          ~1 | paper_ID, 
                                          ~1 | cohort_ID, 
                                          ~1 | species_ID, 
                                          ~1 | phylogeny),
                            R = list(phylogeny = varcor), # added in phylogney
                            data = df)

# orchard_plot(meta_regression3b, mod = "taxa", xlab = "Zr (effect size)")
orchard_plot(meta_regression3c, mod = "taxa", xlab = "Correlation coefficient (r)", transfm = "tanh", cb="TRUE")

# Sex
# reordering
df$sex <- factor(df$sex, levels = c("F", "M", "B"))

meta_regression4 <- rma.mv(yi = Zr, V = VCV, 
                           mods = ~ sex,
                           random = list(~1 | es_ID, 
                                         ~1 | paper_ID, 
                                         ~1 | cohort_ID, 
                                         ~1 | species_ID, 
                                         ~1 | phylogeny),
                           R = list(phylogeny = varcor), # added in phylogney
                           data = df)

summary(meta_regression4)

meta_regression4b <- rma.mv(yi = Zr, V = VCV, 
                            mods = ~ relevel(sex, ref = "F"),
                            random = list(~1 | es_ID, 
                                          ~1 | paper_ID, 
                                          ~1 | cohort_ID, 
                                          ~1 | species_ID, 
                                          ~1 | phylogeny),
                            R = list(phylogeny = varcor), # added in phylogney
                            data = df)

summary(meta_regression4b)

r2_ml(meta_regression4)

# TODO for orchard plot - need meta-regression without intercept (see vignette)
meta_regression4c <- rma.mv(yi = Zr, V = VCV, 
                            mods = ~ sex -1,
                            random = list(~1 | es_ID, 
                                          ~1 | paper_ID, 
                                          ~1 | cohort_ID, 
                                          ~1 | species_ID, 
                                          ~1 | phylogeny),
                            R = list(phylogeny = varcor), # added in phylogney
                            data = df)

# orchard_plot(meta_regression4b, mod = "sex", xlab = "Zr (effect size)")
orchard_plot(meta_regression4c, mod = "sex", xlab = "Correlation coefficient (r)", transfm = "tanh", cb="TRUE")


# Can I look at if repeatability decreases with the number of observations per individual ?
# Using k ?

meta_regression5 <- rma.mv(yi = Zr, V = VCV, 
                            mods = ~ k,
                            random = list(~1 | es_ID, 
                                          ~1 | paper_ID, 
                                          ~1 | cohort_ID, 
                                          ~1 | species_ID, 
                                          ~1 | phylogeny),
                            R = list(phylogeny = varcor), # added in phylogney
                            data = df)

# do I need to centre? Doesn't seem to make a difference
df$k.c <- as.vector(scale(df$k, scale = F))

meta_regression5 <- rma.mv(yi = Zr, V = VCV, 
                           mods = ~ k.c,
                           random = list(~1 | es_ID, 
                                         ~1 | paper_ID, 
                                         ~1 | cohort_ID, 
                                         ~1 | species_ID, 
                                         ~1 | phylogeny),
                           R = list(phylogeny = varcor), # added in phylogney
                           data = df)

summary(meta_regression5)

# Plotting using code from Hayward et al. 2021 
# but using Zr, not ICC or R
pred_meta_regression5 <- predict.rma(meta_regression5)

fig_meta_regression5 <-  df %>% 
  mutate(ymin = pred_meta_regression5$ci.lb, 
         ymax = pred_meta_regression5$ci.ub,
         ymin2 = pred_meta_regression5$cr.lb,
         ymax2 = pred_meta_regression5$cr.ub,
         pred = pred_meta_regression5$pred) %>% 
  ggplot(aes(x = (k), y = Zr, size = (1/VZr) + 3, )) +
  geom_point(shape = 21, fill = "grey90") +
  geom_smooth(aes(y = ymin2), method =  "loess", se = FALSE, lty =  "dotted", lwd = 0.25, colour = "#0072B2") +
  geom_smooth(aes(y = ymax2), method =  "loess", se = FALSE, lty = "dotted", lwd = 0.25, colour = "#0072B2") +
  geom_smooth(aes(y = ymin), method =  "loess", se = FALSE,lty = "dotted", lwd = 0.25, colour ="#D55E00") +
  geom_smooth(aes(y = ymax), method =  "loess", se = FALSE, lty ="dotted", lwd = 0.25, colour ="#D55E00") + 
  geom_smooth(aes(y = pred), method =  "loess", se = FALSE, lty ="dashed", lwd = 0.5, colour ="black")

# Model selection (multi-predictor model) ####

# creates a new function to run in MuMIn
updated.rma.mv <- updateable(rma.mv)
# updated.rma.mv

# testing the new function use method = 'ML' so that we can compare AIC
mr_full <- updated.rma.mv(yi = Zr, V = VCV, 
                          mods = ~ method + taxa + sex + annual_event,
                          random = list(~1 | es_ID, 
                                        ~1 | paper_ID, 
                                        ~1 | cohort_ID, 
                                        ~1 | species_ID, 
                                        ~1 | phylogeny),
                          R = list(phylogeny = varcor), # phylogenetic matrix
                          method = "ML", 
                          data = df)


### additional methods for "rma.mv" class (made by Kamil Barton)
### we need this to run model selection with rma.mv in MuMIn

formula.rma.mv <- function (x, ...) return(eval(getCall(x)$mods))

makeArgs.rma.mv <-
  function (obj, termNames, comb, opt, ...) {
    ret <- MuMIn:::makeArgs.default(obj, termNames, comb, opt)
    names(ret)[1L] <- "mods"
    ret
  }

nobs.rma.mv <-
  function (object, ...)
    attr(logLik(object), "nall")

coefTable.rma.mv <- function (model, ...)
  MuMIn:::.makeCoefTable(model$b, model$se, coefNames = rownames(model$b))

# testing dredge dredge(full.model, evaluate=F) # show all candidate models n = 16
candidates <- dredge(mr_full)

# displays delta AICc <2
candidates_aic2 <- subset(candidates, delta < 2)

# model averaging it seems like models are using z values rather than t values
mr_averaged_aic2 <- summary(model.avg(candidates, delta < 2))

# relative importance of each predictor
importance <- importance(candidates)

# use REML if not for model comparision
model1 <- rma.mv(yi = Zr, V = VCV, 
                 mods = ~  annual_event + sex,
                 random = list(~1 | es_ID, 
                               ~1 | paper_ID, 
                               ~1 | cohort_ID, 
                               ~1 | species_ID, 
                               ~1 | phylogeny),
                 R = list(phylogeny = varcor), # phylogenetic matrix
                 method="REML", 
                 data = df)
model2 <- rma.mv(yi = Zr, V = VCV, 
                 mods = ~  annual_event,
                 random = list(~1 | es_ID, 
                               ~1 | paper_ID, 
                               ~1 | cohort_ID, 
                               ~1 | species_ID, 
                               ~1 | phylogeny),
                 R = list(phylogeny = varcor), # phylogenetic matrix
                 method="REML", 
                 data = df)
model3 <- rma.mv(yi = Zr, V = VCV, 
                 mods = ~  annual_event + taxa,
                 random = list(~1 | es_ID, 
                               ~1 | paper_ID, 
                               ~1 | cohort_ID, 
                               ~1 | species_ID, 
                               ~1 | phylogeny),
                 R = list(phylogeny = varcor), # phylogenetic matrix
                 method="REML", 
                 data = df)
model4 <- rma.mv(yi = Zr, V = VCV, 
                 mods = ~  annual_event + method,
                 random = list(~1 | es_ID, 
                               ~1 | paper_ID, 
                               ~1 | cohort_ID, 
                               ~1 | species_ID, 
                               ~1 | phylogeny),
                 R = list(phylogeny = varcor), # phylogenetic matrix
                 method="REML", 
                 data = df)
model5 <- rma.mv(yi = Zr, V = VCV, 
                 mods = ~  annual_event + sex + taxa,
                 random = list(~1 | es_ID, 
                               ~1 | paper_ID, 
                               ~1 | cohort_ID, 
                               ~1 | species_ID, 
                               ~1 | phylogeny),
                 R = list(phylogeny = varcor), # phylogenetic matrix
                 method="REML", 
                 data = df)

# getting averaged R2 and variance components not provided by the MuMIn package - cannot get this to work ??
# average_sigma2 <- weighted.mean(x = c(model1$sigma2, model2$sigma2, model3$sigma2, model4$sigma2, model5$sigma2), w = candidates_aic2$weight)
# average_R2 <- weighted.mean(x = c(R2(model1)[1], R2(model2)[1], R2(model3)[1], R2(model4)[1], R2(model5)[1]) , w = candidates_aic2$weight)


# TODO - please check for publication bias and time lag bias
# read this paper - https://ecoevorxiv.org/k7pmz
# here is the associated code - https://github.com/itchyshin/publication_bias

# Publication bias 
# To test for publication bias, first fit a phylogenetic multilevel meta-regression to explore whether there is some evidence of small-study effects in meta-analytic dataset. Fit a uni-moderator phylogenetic multilevel meta-regression including the effect sizes' standard errors (SE or sei) as the only moderator (see Equation 21 from main text). This meta-regression will provide some information about the existence of small-study effects (i.e. asymmetry in the distribution of effect sizes)

# creating a variable for the standard error of each effect size (i.e. the square root of the sampling variance)
df$sei <- sqrt(df$VZr)

# Application of Equation 21 from the main text
publication.bias.model.r.se <- rma.mv(yi = Zr, V = VCV,
                                      mod = ~1 + sei,
                                      random = list(~1 | es_ID, 
                                                    ~1 | paper_ID, 
                                                    ~1 | cohort_ID, 
                                                    ~1 | species_ID, 
                                                    ~1 | phylogeny),
                                      R = list(phylogeny = varcor), # added in phylogney
                                      data=df)

print(publication.bias.model.r.se,digits=3)
#Slope = 0.182. Not significant. Means that effect sizes with larger SE (more uncertain effect sizes) DO NOT tend to be larger
estimates.publication.bias.model.r.se <- estimates.CI(publication.bias.model.r.se)
round(estimates.publication.bias.model.r.se[2,2],2)
round(orchaRd::r2_ml(publication.bias.model.r.se)[[1]]*100,1)


# Time-lag bias
# To test for time-lag bias (also called decline effects) we can first fit a uni-moderator phylogenetic multilevel meta-regression including the year of publication (mean-centred) as the only moderator (see Equation 23 from main text). The estimated slope for year of publication will provide some evidence on whether effect sizes have changed linearly over time since the first effect size was published

df$pub_year <- as.numeric(df$pub_year)
df$pub_year.c <- as.vector(scale(df$pub_year, scale = F)) # two rows from my unpublished paper on RI petrels

# Application of Equation 23 from the main manuscript
publication.bias.model.r.timelag <- rma.mv(yi = Zr, V = VCV,
                                           mods= ~1 + pub_year.c,
                                           random = list(~1 | es_ID, 
                                                         ~1 | paper_ID, 
                                                         ~1 | cohort_ID, 
                                                         ~1 | species_ID, 
                                                         ~1 | phylogeny),
                                           R = list(phylogeny = varcor), # added in phylogney
                                           data=df)

summary(publication.bias.model.r.timelag) # estimate of pub year v close to zero
round(orchaRd::r2_ml(publication.bias.model.r.timelag)[[1]]*100,1)


# All-in publication bias test (multi-moderator)  
# once happy with what moderators to include in final model ?

publication.bias.model.r.all.se <- rma.mv(yi = Zr, V = VCV,
                                         mods= ~1 + # -1 removes the intercept
                                           sei +
                                           pub_year.c +
                                         annual_event + method + taxa + sex, 
                                         random = list(~1 | es_ID, 
                                                       ~1 | paper_ID, 
                                                       ~1 | cohort_ID, 
                                                       ~1 | species_ID, 
                                                       ~1 | phylogeny),
                                         R = list(phylogeny = varcor), # added in phylogney
                                         data=df)

summary(publication.bias.model.r.all.se)
# sei & pub_year still not significant - no small study effect or publication bias. Same as uni-variate models.

round(orchaRd::r2_ml(publication.bias.model.r.all.se)[[1]]*100,1)





# Make map of locations of repeatability studies #######################################

worldmap <- map_data('world')
#add point for RI
#RI <- data.frame(long= 57.30, lat = -20.11, stringsAsFactors = FALSE) #make point for RI
names(df)
(world <- ggplot() + geom_polygon(data = worldmap, aes(x=long, y = lat, group = group), fill = "grey", 
                                  colour = "grey39", size = 0.5) + 
    coord_fixed(1.3) +
    theme_minimal() +
    scale_colour_manual(values=c("#F9766E", "#619DFF", "#00BB38")) +
    geom_point(data=df, aes(x= long, y =lat, colour=taxa, shape=method), size =6) +
    geom_rect(aes(xmin= -12, xmax= 23, ymin= 32, ymax= 64), fill = NA, colour = "black", size = 1.2) +
    #geom_point(data=RI, aes(x=long, y=lat), shape=24, fill= "yellow", size = 12) +
    theme(panel.grid = element_line(colour = "white"), axis.title = element_blank(), axis.text = element_blank(), 
          axis.ticks = element_blank(), legend.position = "none"))

##' Zoomed in plot on Europe
(zoom <- ggplot() + geom_polygon(data = worldmap, aes(x=long, y = lat, group = group), fill = "grey", 
                                 colour = "grey39", size = 0.5) + 
    coord_sf(xlim= c(-11,23), ylim= c(32,64)) +
    theme_minimal() +
    theme(legend.text=element_text(size=13), legend.title=element_text(size=15, face="bold")) +
    theme(panel.grid = element_line(colour= "white"),
          axis.title = element_blank(), axis.text = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA), 
          axis.ticks = element_blank()) +
    scale_colour_discrete(name = "Ecological Group", labels=c("Landbird", "Waterbird", "Seabird")) +
    scale_shape_discrete(name = "Method") +
    geom_point(data=df, aes(x= long, y =lat, shape=method, colour=taxa),size =6))

plot_grid(world, zoom, nrow=1, rel_widths = c(1.5,1))

ggsave("figs/Map_of_tagging_locations.jpg")

