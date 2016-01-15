
# function to print the param estimates (i.e. error rater matrices, true label matrix, etc.)
print_estimates <- function(param_list, max_param_index){
    max_param <- param_list[["max_params_per_run"]][[paste("run", max_param_index, sep = "_")]]
    # print the error rate for every rater
    for(i in 1:length(max_param[["error_rate_list"]])){
        print(paste("error rate for rater", i, sep = "_"))
        print(round(max_param[["error_rate_list"]][[i]], 3))
        print(rep("-",10))
    }
    
    # print the probability of the categories/label types
    print("probability of categories/label types:")
    print(round(max_param[["prob_label_types"]], 3))
    print(rep("-",10))
    
    # print the estimate true label probability (i.e. consensus ratings)
    print("probability of true label for each rated object:")
    print(round(max_param[["truelabel_mat"]], 4))
    print(rep("-",10))
    
    # plot the loglikelihood diagnostic to check if it was increasing at each iteration
    plot(param_list[["loglikelihood_diagnostic"]][[paste("run", max_param_index, sep = "_")]])
}

#import the implementation of Dawid et al. (1979) algorithm
working_dir <- getwd()

source(paste(working_dir, "em_dawid79.R", sep = "/"))

convergence_threshold <- 0.000001

# data labels for 5 hypothetical raters to rate if the url belongs to a search engine company or not
# labels coding 1: NO (it is not a search engine company)
#               2: YES (it is a search engine company)
url <- c("google.com","bing.com","youtube.com","baidu.com","yale.edu")
r1 <- c(1,1,1,1,1)
r2 <- c(2,1,1,2,1)
r3 <- c(2,2,1,2,1)
r4 <- c(2,2,1,2,1)
r5 <- c(1,1,2,1,2)
test_data <- data.frame(url,r1,r2,r3,r4,r5)
rating_categories <- 1:2
# get the consensus ratings using the data to calculate the initial estimates
# of the probability of the true labels T_{ij}
consensus_syntheticdata <- find_consensus_rating(test_data[,-1],
                                                 rating_categories,
                                                 num_runs = 1,
                                                 num_iter = 100,
                                                 convergence_threshold)
# get the index of the parameters that have highest loglikelihood
max_param_index <- consensus_syntheticdata[["max_param_index"]]
# if muliple parameters have the highest loglikilihood, get the first one
# it is more reasonable in this scenario to check each of the parameter estimates to decide which to choose
if(length(max_param_index) > 1){
    max_param_index <- max_param_index[1]
}

# print the output of the model
print_estimates(consensus_syntheticdata, max_param_index)

#---------------------------------------------------------
# Dawid et al. (1979) data set
dataset_dir <- paste(working_dir, "dataset", sep = "/")
dawid_em_data <- read.table(file = paste(dataset_dir, "dawid_em_data.csv", sep = "/"),
                            header = TRUE,
                            sep = ",",
                            quote = "",
                            stringsAsFactors = FALSE)
rating_categories <- 1:4
consensus_dawid_data<- find_consensus_rating(dawid_em_data[,-1], 
                                             rating_categories,
                                             num_runs = 1,
                                             num_iter = 100,
                                             convergence_threshold,
                                             multiple_per_rater = TRUE)

# get the index of the parameters that have highest loglikelihood
max_param_index <- consensus_dawid_data[["max_param_index"]]
# if muliple parameters have the highest loglikilihood, get the first one
# it is more reasonable in this scenario to check each of the parameter estimates to decide which to choose
if(length(max_param_index) > 1){
    max_param_index <- max_param_index[1]
}

# print the output of the model
print_estimates(consensus_dawid_data, max_param_index)