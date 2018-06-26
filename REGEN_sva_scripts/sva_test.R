#----------------------------------#
#           REGEN Analysis         #
#          Kristoph Naggert        #
#             06/25/2018           #
#----------------------------------#

#-------------------#
#       Setup       #
#-------------------#

# Load in libraries used for analysis
library(sva)
library(pamr)
library(limma)

# Read in the data with header set to TRUE (First row becomes column names)
dataT <- read.table("/Users/knaggert/Desktop/summer_2018/regene_for_sva/REGENE15_full_datanorm.txt", header = TRUE)
designT <- read.table("/Users/knaggert/Desktop/summer_2018/regene_for_sva/REGENE15_design.dat", header = TRUE)

dataT14 <- read.table("/Users/knaggert/Desktop/summer_2018/regene_for_sva/REGENE14_full_datanorm.txt", header = TRUE)
designT14 <- read.table("/Users/knaggert/Desktop/summer_2018/regene_for_sva/REGENE14_design.dat", header = TRUE)

# Move the gene_name column so that it becomes row_names for both sets
dataT2 <- dataT[,-1]
rownames(dataT2) <- dataT[,1]

dataT14.2 <- dataT14[,-1]
rownames(dataT14.2) <- dataT14[,1]

# Reverse the log2 transformation for both the 14 and 15 datasets
for (i in 1:ncol(dataT2)) {
  for (j in 1:nrow(dataT2)) {
    n <- 2^dataT2[j,i]
    dataT2[j,i] <- n
  }
}

for (i in 1:ncol(dataT14.2)) {
  for (j in 1:nrow(dataT14.2)) {
    n <- 2^dataT14.2[j,i]
    dataT14.2[j,i] <- n
  }
}

# Filter so that we are working with datasets of between 15,000 and 20,000 elements
filter_T <- apply(dataT2, 1, function(x) length(x[x>4])>=2)
filtered_T <- dataT2[filter_T,]

filter_T14 <- apply(dataT14.2, 1, function(x) length(x[x>4])>=2)
filtered_T14 <- dataT14.2[filter_T,]

# Make the dataframes into a matrix
dat0_test <- as.matrix(filtered_T)

dat0_14test <- as.matrix(filtered_T14)

# The group variable in this case consists of our two different treatments
group = as.factor(rep(c("Ctl", "Trt"), each=3))

# Create the null and alternative models for 15
mod1_test = model.matrix(~group)
mod0_test = cbind(mod1_test[, 1])

# Create the null and alternative models for 14
mod1_test14 = model.matrix(~group)
mod0_test14 = cbind(mod1_test14[, 1])

# Run svaseq to determine the number of SV for each dataset (15 and 14 respectively)
svseq_test = svaseq(dat0_test, mod1_test, mod0_test)

svseq_test14 = svaseq(dat0_14test, mod1_test14, mod0_test14)

# Bind the SV to the alternative model so we can compare it to the alternative model that does not include SVs
modSv_test = cbind(mod1_test, svseq_test$sv)

modSv_test14 = cbind(mod1_test14, svseq_test14$sv)

# The f.pvalue function can be used to calculate parametric F-test p-values for each row of a data matrix
# The F-test compares the models modSv and mod1. They must be nested models, so all of the variables in mod1 must appear in modSv.
# p.adjust simply gives us the adjusted p-values for each row of the data matrix
wacky_p = f.pvalue(dat0_test, modSv_test, mod1_test)
wacky_q = p.adjust(wacky_p, method = "BH")

whacky_p = f.pvalue(dat0_14test, modSv_test14, mod1_test14)
whacky_q = p.adjust(whacky_p, method = "BH")

# Save the adj. p-values as matrices
qh <- as.matrix(whacky_q)
qq <- as.matrix(wacky_q)

# Combine the two matrices of adj. p-values into one larger matrix
outmatrix_q <- cbind(qh, qq)

# Write that combined matrix to .csv file
write.csv(outmatrix_q, file = "/Users/knaggert/Desktop/summer_2018/regene_for_sva/REGENE14_15_adjP.csv")
