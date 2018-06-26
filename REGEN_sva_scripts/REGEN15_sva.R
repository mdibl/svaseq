#----------------------------------#
#           REGEN Analysis         #
#          Kristoph Naggert        #
#             06/21/2018           #
#----------------------------------#

#-------------------#
#       Setup       #
#-------------------#

# Load in libraries for analysis
library(sva)
library(pamr)
library(limma)

# Read in the data with header set to TRUE (First row becomes column names)
data15 <- read.table("/Users/knaggert/Desktop/summer_2018/regene_for_sva/REGENE15_full_datanorm.txt", header = TRUE)
design15 <- read.table("/Users/knaggert/Desktop/summer_2018/regene_for_sva/REGENE15_design.dat", header = TRUE)

# Move the gene_name column so that it becomes row_names
data15.2 <- data15[,-1]
rownames(data15.2) <- data15[,1]

#-------------------------#
#       Filter Data       #
#-------------------------#

# Reverse the log2 transformation (Loop through each entry in the dataframe and preform the transformation)
for (i in 1:ncol(data15.2)) {
  for (j in 1:nrow(data15.2)) {
    n <- 2^data15.2[j,i]
    data15.2[j,i] <- n
  }
}

# Filter so that we are working with a dataset of between 15,000 and 20,000 elements
# The filter is such that there must be a least 2 sample with at least 4 reads for a given gene
filter_15 <- apply(data15.2, 1, function(x) length(x[x>4])>=2)
filtered_15 <- data15.2[filter_15,]

# Make the dataframe into a matrix so that we can run sva functions on it
dat0_15 <- as.matrix(filtered_15)


#--------------------------#
#       Applying SVA       #
#--------------------------#

# Since the outcome of the sva function is varaible we set a seed so that we always get the same output
set.seed(18)

# Create the full model matrix - including both the adjustment variables and the variable of interest
mod_15 = model.matrix(~as.factor(TERM_1), data = design15)

# The null model contains only the adjustment variables. Since we are not ad- justing for any other variables 
# in this analysis, only an intercept is included in the model.
mod0_15 = model.matrix(~1, data = design15)

# The function returns a list with four components, sv, pprob.gam, pprob.b, n.sv.
# sv is a matrix whose columns correspond to the estimated surrogate variables.
# pprob.gam is the posterior probability that each gene is associated with one or more latent variables
# pprob.b is the posterior probability that each gene is associated with the variables of interest
# n.sv is the number of surrogate variables estimated by the sva.

# sva 5 times and compile results for interpretation
svobj_15 = sva(dat0_15, mod_15, mod0_15)
svobj2_15 = sva(dat0_15, mod_15, mod0_15)
svobj3_15 = sva(dat0_15, mod_15, mod0_15)
svobj4_15 = sva(dat0_15, mod_15, mod0_15)
svobj5_15 = sva(dat0_15, mod_15, mod0_15)


#---------------------------------#
#       Dump svovj to table       #
#---------------------------------#

# Turn the probabilty doubles into matrices to combine into one larger matrix
a15 <- as.matrix(svobj_15$pprob.gam)
b15 <- as.matrix(svobj2_15$pprob.gam)
c15 <- as.matrix(svobj3_15$pprob.gam)
d15 <- as.matrix(svobj4_15$pprob.gam)
e15 <- as.matrix(svobj5_15$pprob.gam)

# Combine the matrices together
outmatrix_15 <- cbind(a15, b15, c15, d15, e15)

# Add a column with the gene names
outmatrix_15 <- cbind(outmatrix_15, rownames(filtered_15))

# Write the matrix to a csv file
write.csv(outmatrix_15, file = "/Users/knaggert/Desktop/summer_2018/regene_for_sva/sva_pprob.gam_dump_REGEN15.csv")

#-----------------------------#
#       Applying SVASeq       #
#-----------------------------#

# Load in color library and set a color scheme
library(RColorBrewer)
colors <- brewer.pal(3, "Set2")

# The group variable in this case consists of our two different treatments
group = as.factor(rep(c("Ctl", "Trt"), each=3))

# Set null and alternative models
mod1_15 = model.matrix(~group)
mod0_15 = cbind(mod1[, 1])

# Now we can apply svaseq to estimate the latent factor. We don't specify n.sv, because
# svaseq will calculate it for us if left undefined
svseq_15 = svaseq(dat0_15, mod1_15, mod0_15)

# Plot the sv and the weights for each sample
plot(svseq_15, pch=19, col = colors[group])





