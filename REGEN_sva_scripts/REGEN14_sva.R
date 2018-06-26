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
data14 <- read.table("/Users/knaggert/Desktop/summer_2018/regene_for_sva/REGENE14_full_datanorm.txt", header = TRUE)
design14 <- read.table("/Users/knaggert/Desktop/summer_2018/regene_for_sva/REGENE14_design.dat", header = TRUE)

# Move the gene_name column so that it becomes row_names
data14.2 <- data14[,-1]
rownames(data14.2) <- data14[,1]

#-------------------------#
#       Filter Data       #
#-------------------------#

# Rearrange columns so that the three controls are together and the three treated samples are together
data14.3 <- data14.2[,c("REGENS.248", "REGENS.249", "REGENS.253", "REGENS.250", "REGENS.251", "REGENS.252")]

# Reverse the log2 transformation (Loop through each entry in the dataframe and preform the transformation)
for (i in 1:ncol(data14.3)) {
  for (j in 1:nrow(data14.3)) {
    n <- 2^data14.3[j,i]
    data14.3[j,i] <- n
  }
}

# Filter so that we are working with a dataset of between 15,000 and 20,000 elements
# The filter is such that there must be a least 2 sample with at least 4 reads for a given gene
filter_14 <- apply(data14.3, 1, function(x) length(x[x>4])>=2)
filtered_14 <- data14.3[filter_14,]

# Make the dataframe into a matrix so that we can run sva functions on it
dat0_14 <- as.matrix(filtered_14)

#--------------------------#
#       Applying SVA       #
#--------------------------#

# Since the outcome of the sva function is varaible we set a seed so that we always get the same output
set.seed(444)

# Create the full model matrix - including both the adjustment variables and the variable of interest
mod_14 = model.matrix(~as.factor(TERM_1), data = design14)

# The null model contains only the adjustment variables. Since we are not ad- justing for any other variables 
# in this analysis, only an intercept is included in the model.
mod0_14 = model.matrix(~1, data = design14)

# The function returns a list with four components, sv, pprob.gam, pprob.b, n.sv.
# sv is a matrix whose columns correspond to the estimated surrogate variables.
# pprob.gam is the posterior probability that each gene is associated with one or more latent variables
# pprob.b is the posterior probability that each gene is associated with the variables of interest
# n.sv is the number of surrogate variables estimated by the sva.

# sva 5 times and compile results for interpretation
svobj_14 = sva(dat0_14, mod_14, mod0_14)

svobj2_14 = sva(dat0_14, mod_14, mod0_14)

svobj3_14 = sva(dat0_14, mod_14, mod0_14)

svobj4_14 = sva(dat0_14, mod_14, mod0_14)

svobj5_14 = sva(dat0_14, mod_14, mod0_14)

svobj6_14 = sva(dat0_14, mod_14, mod0_14)

svobj7_14 = sva(dat0_14, mod_14, mod0_14)

#---------------------------------#
#       Dump svovj to table       #
#---------------------------------#

# 2/7 svobj test returned 0 significant SVs (Tests 1 and 4)

# Turn the probabilty doubles into matrices to combine into one larger matrix
a14 <- as.matrix(svobj2_14$pprob.gam)
b14 <- as.matrix(svobj3_14$pprob.gam)
c14 <- as.matrix(svobj5_14$pprob.gam)
d14 <- as.matrix(svobj6_14$pprob.gam)
e14 <- as.matrix(svobj7_14$pprob.gam)

# Combine the matrices together
outmatrix_14 <- cbind(a14, b14, c14, d14, e14)

# Add a column with the gene names
outmatrix_14 <- cbind(outmatrix_14, rownames(filtered_14))

# Write the matrix to a csv file
write.csv(outmatrix_14, file = "/Users/knaggert/Desktop/summer_2018/regene_for_sva/sva_pprob.gam_dump_REGEN14.csv")

#-----------------------------#
#       Applying SVASeq       #
#-----------------------------#

# Load in color library and set a color scheme
library(RColorBrewer)
colors <- brewer.pal(3, "Set2")

# The group variable in this case consists of our two different treatments (DMSO and ZF143)
group = as.factor(rep(c("Ctl", "Trt"), each = 3))

# Set null and alternative models
mod1_14 = model.matrix(~group)
mod0_14 = cbind(mod1[, 1])

# Now we can apply svaseq to estimate the latent factor. We don't specify n.sv, because
# svaseq will calculate it for us if left undefined
svseq_14 = svaseq(dat0_14, mod1_14, mod0_14)

# Plot the sv and the weights for each sample
plot(svseq_14, pch=19, col = colors[group])

