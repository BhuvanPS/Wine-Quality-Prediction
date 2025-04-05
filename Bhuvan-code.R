################################# 
# You can use this template to draft the script for your Assessment 2 of SIT718.
# More clarification and related resources can be found at
# https://d2l.deakin.edu.au/d2l/le/content/1527306/viewContent/7891667/View
#################################

#############################################################################################
# save your code as "name-code.R" (where ''name'' is replaced with your surname or first name).
#############################################################################################

##################################
#Question 1 - Understand the Data
##################################

data.raw <- as.matrix(read.table("RedWine.txt"))
ID<-224113776
set.seed(ID) # using your student ID number for reproducible sampling with the seed function


data.subset <- data.raw[sample(1:1599, 460), c(1:6)]

data.variable.names <- c("citric acid", "chlorides", "total sulfur dioxide", "pH", "alcohol","quality")
colnames(data.subset)<-data.variable.names
min(data.subset[,5])
min(data.subset[,6])
max(data.subset[,6])
library(graphics) # Base R graphics



library(ggplot2)

# First, let's explore the distributions of our variables
par(mfrow=c(2,3))
# Create 5 scatterplots function (for each X variable against the variable of interest Y) 

for(i in 1:5){
  plot(data.subset[,i],data.subset[,6],xlab = data.variable.names[i],ylab = "Quality",main=paste(data.variable.names[i],"vs Quality"))
}
# Create 6 histograms for each X variable and Y
par(mfrow=c(1,1))
for(i in 1:6){
  hist(data.subset[,i],xlab = data.variable.names[i],main = paste("Histogram of ",data.variable.names[i]))
}


################################
#Question 2 - Transform the Data
################################


I<-c("citric acid", "chlorides", "total sulfur dioxide", "alcohol","quality") # Choose any four X variables and Y
ncol(data.subset)
length(I)
variables_for_transform <- data.subset[,I]  # obtain a 460 by 5 matrix

# for each variable, you need to figure out a good data transformation method, 
# such as Polynomial, log and negation transformation. The k-S test and Skewness 
# calculation may be helpful to select the transformation method
library(moments)

analyze_normality <- function(x, var_name) {
  mean_val <- mean(x)
  median_val <- median(x)
  mode_val <- as.numeric(names(sort(table(round(x, 2)), decreasing = TRUE)[1])) # Approximate mode
  
  if (length(mode_val) == 0 || is.na(mode_val) || is.nan(mode_val)) {
    mode_val <- NA # Or some other default value
  }
  
  skewness_val <- skewness(x)
  
  return(c(mean_val, median_val, mode_val, skewness_val))
}
normality_results_before <- data.frame(
  "citric acid" = analyze_normality(variables_for_transform[, 1], "citric acid"),
    "chlorides" = analyze_normality(variables_for_transform[, 2], "chlorides"),
    "total sulfur dioxide" = analyze_normality(variables_for_transform[, 3], "total sulfur dioxide"),
    "alcohol" = analyze_normality(variables_for_transform[, 4], "alcohol"),
    "quality" = analyze_normality(variables_for_transform[, 5], "quality")
)
rownames(normality_results_before) <- c("Mean", "Median", "Mode", "Skewness")
print(normality_results_before)

minmax <- function(x) {
    (x - min(x)) / (max(x) - min(x))
  }
transform <- function(data) {
  transformed_data <- matrix(0, nrow = nrow(data), ncol = 5)

  # chlorides(X2)
  transformed_data[, 2] <- log1p(log1p(data[, 2]^0.001))

  # total sulfur dioxide (X3) - log1p transformation
  transformed_data[, 3] <- log1p(data[, 3])

  # citric acid (X1)
  transformed_data[, 1] <- data[, 1]^0.6
print(data[8,4])
  # alcohol (X4) - log1p transformation
  transformed_data[, 4] <- log1p(log1p(data[, 4]))

  # quality (Y) - No transformation (ordinal)
  transformed_data[, 5] <- data[, 5]
write.table(transformed_data, "unscaled-transformed.txt")
  # Min-Max Scaling
  

  return(transformed_data)
}

data.transformed <- transform(variables_for_transform)
  data.transformed[, 1:5] <- apply(data.transformed[, 1:5], 2, minmax)


# Normality analysis after scaling
normality_results_after_scaling <- data.frame(
  "citric acid" = analyze_normality(data.transformed[, 1], "citric acid"),
    "chlorides" = analyze_normality(data.transformed[, 2], "chlorides"),
    "total sulfur dioxide" = analyze_normality(data.transformed[, 3], "total sulfur dioxide"),
    "alcohol" = analyze_normality(data.transformed[, 4], "alcohol"),
    "quality" = analyze_normality(data.transformed[, 5], "quality")
)
rownames(normality_results_after_scaling) <- c("Mean", "Median", "Mode", "Skewness")
print("Normality Analysis of Scaled Transformed Data:")
print(normality_results_after_scaling)



# Save this transformed data to a text file
write.table(data.transformed, "name-transformed.txt")  # replace ??name?? with either your surname or first name.


##########################################
#Question 3 - Build models and investigate
##########################################

source("AggWaFit718.R")

data.transformed_copy <- as.matrix(read.table("name-transformed.txt"))  # import your saved data
cor(data.transformed_copy)
# Get weights for Weighted Arithmetic Mean with fit.QAM() 
wam.result <- fit.QAM(cbind(data.transformed_copy[,1:4],data.transformed_copy[,5]), stats.1 = "wam_stats.txt")

# Get weights for Power Mean p=0.5 with fit.QAM()

result.wpm05 <- fit.QAM(cbind(data.transformed_copy[,1:4],data.transformed_copy[,5]), g = PM05, g.inv = invPM05, stats.1 = "wpm05_stats.txt")
# Get weights for Power Mean p=2 with fit.QAM()
result.wpm2 <- fit.QAM(cbind(data.transformed_copy[,1:4], data.transformed_copy[,5]),g=QM ,stats.1 = "wpm2_stats.txt")

# Get weights for Ordered Weighted Average with fit.OWA()
OWA_weights <- fit.OWA(cbind(data.transformed_copy[,1:4], data.transformed_copy[,5]),stats.1 = "owa_stats.txt")
print(OWA_weights)
#choquet
choquet<-fit.choquet(cbind(data.transformed_copy[,1:4], data.transformed_copy[,5]),stats.1 = "choquet.txt")


# Read the output file
# Function to extract stats from the output file
extract_stats <- function(file_path, weights_start_index, weights_end_index) {
  output_data <- readLines(file_path)

  rmse <- as.numeric(gsub("RMSE ", "", output_data[1]))
  avg_abs_error <- as.numeric(gsub("Av. abs error ", "", output_data[2]))
  pearson_corr <- as.numeric(gsub("Pearson correlation ", "", output_data[3]))
  spearman_corr <- as.numeric(gsub("Spearman correlation ", "", output_data[4]))

  weights_lines <- output_data[weights_start_index:weights_end_index]
  weights <- sapply(weights_lines, function(line) {
    as.numeric(strsplit(line, " ")[[1]][2])
  })

  weights_df <- data.frame(
    i = 1:length(weights),
    w_i = weights
  )
rownames(weights_df) <- c("citric acid","chlorides","total sulfur dioxide","alcohol")

  print(paste("RMSE:", rmse))
  print(paste("Av. abs error:", avg_abs_error))
  print(paste("Pearson correlation:", pearson_corr))
  print(paste("Spearman correlation:", spearman_corr))
  print(weights_df)
}

# Call the function for wam_stats.txt, with weights starting from line 6 and ending at line 9.
extract_stats("wam_stats.txt", 6, 9)
extract_stats("wpm05_stats.txt", 6, 9)
extract_stats("wpm2_stats.txt", 6, 9)
#######################################
#Question 4 - Use Model for Prediction
#######################################

new_input <- c(0.8, 0.63, 37, 2.51, 7.0) 

load_wam_weights <- function(file_path) {
  output_data <- readLines(file_path)
  weights_lines <- output_data[6:9]  # Assuming weights start from line 6 to 9
  weights <- sapply(weights_lines, function(line) {
    as.numeric(strsplit(line, " ")[[1]][2])
  })
  return(weights)
}

transformed_new_input_no_scale <- data.frame(
  X1 = as.numeric(new_input[1])^0.6,
  X2 = log1p(log1p(as.numeric(new_input[2])^0.001)),
  X3 = log1p(as.numeric(new_input[3])),
  X4 = log1p(log1p(as.numeric(new_input[5])))
)
print(transformed_new_input_no_scale)
transformed_train_data <- read.table("unscaled-transformed.txt")
train_min <- apply(transformed_train_data[, 1:4], 2, min)
train_max <- apply(transformed_train_data[, 1:4], 2, max)
print(train_min)
print(train_max)
# Scale the new input using the training data's min and max.
scaled_new_input <- (transformed_new_input_no_scale - train_min) / (train_max - train_min)
print(scaled_new_input)
# Load the weights
wam_weights <- load_wam_weights("wam_stats.txt")

# Calculate the predicted quality using the Weighted Arithmetic Mean
predicted_quality <- sum(as.numeric(scaled_new_input[1, 1:4]) * wam_weights)
print(predicted_quality) 
print(predicted_quality*5 +3)
scaled_new_input <- apply(transformed_new_input, 2, minmax)
print(transformed_new_input)
# Make predictions using the loaded models
# Assuming AggWaFit718.R provides prediction functions or weights
# Adapt this part based on the actual output of AggWaFit718.R
x<-c(0.290573580053069,0.194675814651332,0.0432655706497103,0.471485034645866)

predicted_wam <- sum(transformed_new_input*x)
print(predicted_wam*5 +3)
print(max(data.subset[,6]))
print(min(data.subset[,6]))

predicted_wpm_0.5 <- weighted.mean(scaled_new_input, wpm.weights.0.5)
predicted_wpm_2 <- weighted.mean(scaled_new_input, wpm.weights.2)
predicted_owa <- owa(scaled_new_input, owa.weights) 

# transforming the four variables in the same way as in question 2 



# applying the transformed variables to the best model selected from Q3 for Y prediction



# Reverse the transformation to convert back the predicted Y to the original scale and then round it to integer




#############################################################################################
# References 
# Following Harvard style: https://www.deakin.edu.au/students/studying/study-support/referencing/harvard
#############################################################################################

# You must cite all the datasets and packages you used for this assessment. 
#
#