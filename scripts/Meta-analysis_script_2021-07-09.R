## Repeatability in avian migratory timings

devtools::install_github("itchyshin/orchard_plot", subdir = "orchaRd", force = TRUE, build_vignettes = TRUE)


# Load packages #######################################
pacman::p_load(SciViews,
               tidyverse,
               metafor, # package for meta-analysis
               rotl, # package for phylogeny
               ape, # package for phylogeny
               cowplot, # combining multiple plots
               here, # making folder path usable for all
               clubSandwich, # package to assist metafor
               orchaRd
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







# Read in and sort data #######################################

df <- read_csv(here("data", "Meta-analysis_data.csv")) # I had to use read_csv ratther than read.csv
str(df)
df

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
                                    ~1 | species_ID, 
                                    ~1 | phylogeny),
  R = list(phylogeny = varcor), # added in phylogney
  data = df)

summary(meta_model1)

## Calculating I^2
round(i2_ml(meta_model1)*100,1) 

# I2_total      I2_es_ID   I2_paper_ID  I2_cohort_ID I2_species_ID  I2_phylogeny 
# 85.9          51.1           0.0           0.0          27.3           7.5 

# TODO you could keep it all but you could take cohort ID  + paper 
# or accoding to % explaning - we can tak out a couple of random effects

meta_model2 <- rma.mv(yi = Zr, V = VCV, 
                      random = list(~1 | es_ID, 
                                    ~1 | paper_ID, 
                                   # ~1 | cohort_ID, 
                                    ~1 | species_ID, 
                                    ~1 | phylogeny
                                   ),
                      R = list(phylogeny = varcor), # added in phylogney
                      data = df)

meta_model3 <- rma.mv(yi = Zr, V = VCV, 
                      random = list(~1 | es_ID, 
                                    #~1 | paper_ID, 
                                    # ~1 | cohort_ID, 
                                    ~1 | species_ID, 
                                    ~1 | phylogeny
                                    ),
                      R = list(phylogeny = varcor), # added in phylogney
                      data = df)

# we cannot take more random effects
# to do mulileve models we need at least one higher level
meta_model4 <- rma.mv(yi = Zr, V = VCV, 
                      random = list(~1 | es_ID, 
                                    #~1 | paper_ID, 
                                    # ~1 | cohort_ID, 
                                    ~1 | species_ID 
                                    #~1 | phylogeny
                                    ),
                      R = list(phylogeny = varcor), # added in phylogney
                      data = df)

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

orchard_plot(meta_model1, xlab = "Zr (effect size)")


# TODO - please do meta-regression
# one example

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

summary(meta_regression1)

r2_ml(meta_regression1)

# TODO for orchard plot (read vignette)
meta_regression1b <- rma.mv(yi = Zr, V = VCV, 
                            mods = ~ method -1,
                            random = list(~1 | es_ID, 
                                          ~1 | paper_ID, 
                                          ~1 | cohort_ID, 
                                          ~1 | species_ID, 
                                          ~1 | phylogeny),
                            R = list(phylogeny = varcor), # added in phylogney
                            data = df)

# TODO - please check for publication bias and time lag bias
# read this paper - https://ecoevorxiv.org/k7pmz
# here is the associated code - https://github.com/itchyshin/publication_bias

orchard_plot(meta_regression1b, mod = "method", xlab = "Zr (effect size)")


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

