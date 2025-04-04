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


I <- c("citric acid","chlorides","total sulfur dioxide","alcohol","quality") # Choose any four X variables and Y
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
  analyze_normality(variables_for_transform[, 1], "citric acid"),
    analyze_normality(variables_for_transform[, 2], "chlorides"),
    analyze_normality(variables_for_transform[, 3], "total sulfur dioxide"),
    analyze_normality(variables_for_transform[, 4], "alcohol"),
    analyze_normality(variables_for_transform[, 5], "quality")
)
rownames(normality_results_before) <- c("Mean", "Median", "Mode", "Skewness")
print("Normality Analysis of Original Data:")
print(normality_results_before)

data.transformed <- matrix(0, nrow = nrow(variables_for_transform), ncol = 5)

#chlorides(X2)
data.transformed[, 2] <- log1p(log1p(variables_for_transform[, 2]^0.001))

# total sulfur dioxide (X3) - Box-Cox transformation

data.transformed[, 3] <- log1p(variables_for_transform[, 3])

# pH (X4) - No transformation
data.transformed[, 1] <- variables_for_transform[, 1]

# alcohol (X5) - log1p transformation
data.transformed[, 4] <- log1p(log1p(variables_for_transform[, 4]))

# quality (Y) - No transformation (ordinal)
data.transformed[, 5] <- variables_for_transform[, 5]

# Min-Max Scaling (applied before z-score) - EXCLUDING quality
minmax <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}
data.transformed[, 1:5] <- apply(data.transformed[, 1:5], 2, minmax) # Apply minmax to columns 1-4

# Print the scaled data
print("Scaled Transformed Data:")
print(data.transformed)

# Normality analysis after scaling
normality_results_after_scaling <- data.frame(
  analyze_normality(variables_for_transform[, 1], "citric acid"),
    analyze_normality(variables_for_transform[, 2], "chlorides"),
    analyze_normality(variables_for_transform[, 3], "total sulfur dioxide"),
    analyze_normality(variables_for_transform[, 4], "alcohol"),
    analyze_normality(variables_for_transform[, 5], "quality")
)
rownames(normality_results_after_scaling) <- c("Mean", "Median", "Mode", "Skewness")
print("Normality Analysis of Scaled Transformed Data:")
print(normality_results_after_scaling)

# Save the scaled transformed data to a text file
write.table(data.transformed, "scaled-transformed.txt")

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


#######################################
#Question 4 - Use Model for Prediction
#######################################

new_input <- c(0.8, 0.63, 37, 2.51, 7.0) 

new_input_for_transform <- new_input[c("index")] # choose the same four X variables as in Q2 

new_X1 <- 0.
new_X2 <- 37
new_X4 <- 2.51
new_X5 <- 7

new_input <- data.frame(X1 = new_X1, X2 = new_X2, X4 = new_X4, X5 = new_X5)

# Apply the same transformations and scaling
transformed_new_input <- data.frame(
  X1 = log1p(log1p(new_input$X1^0.001)),
  X2 = log1p(new_input$X2^0.01),
  X4 = new_input$X4^0.05,
  X5 = log1p(log1p(new_input$X5^0.001))^0.1
)

scaled_new_input <- apply(transformed_new_input, 2, minmax)
print(transformed_new_input)
# Make predictions using the loaded models
# Assuming AggWaFit718.R provides prediction functions or weights
# Adapt this part based on the actual output of AggWaFit718.R
x<-c(0.396366041843691, 0.186332944275711,0, 0.417301013880598)

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