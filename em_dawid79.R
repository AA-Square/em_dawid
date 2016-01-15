# initialize randomly the hidden/missing variables T_{ij} for every i in [1,I] and j in [1,J]
init_truelabel_random <- function(num_objects, rating_categories){
    num_categories <- length(rating_categories)
    obj_truelabel_mat <- matrix(-1, num_objects, num_categories)
    for (i in 1:num_objects){
        obj_truelabel_mat[i, ] <- prop.table(runif(length(rating_categories)))
    }
    return(obj_truelabel_mat)
}

# initialize the hidden/missing variables T_{ij} for every i in [1,I] and j in [1,J] using the counts
# from the responses of every rater. See eq(3.1) in Dawid et al.
init_truelabel_majvote <- function(raters_ratings){
    
    num_raters <- length(raters_ratings)
    num_objects <- nrow(raters_ratings[[1]])
    num_categories <- ncol(raters_ratings[[1]])
    total_counts <- matrix(0, num_objects, num_categories)
    
    for(i in 1:num_raters){
        total_counts <- total_counts + raters_ratings[[i]]
    }
    
    normalized_total_counts <- total_counts / apply(total_counts, 1, FUN = sum)
    
    return(normalized_total_counts)
    
}

# parse multiple ratings provided by the same rater for the same object using regex
parse_concat_ratings <- function(concat_rating){
    
    # parameters:
    #   concat_rating: concatenated ratings by a rater for the same object. i.e. "1111"
    # return:
    #   rating_vec: an integer vector of the ratings i.e. [1,1,1]
    pattern <- "[0-9]"
    rating_vec <- strapply(concat_rating, pattern, c)[[1]]
    return(as.integer(rating_vec))
}

# parse_concat_num <- function(x){
#     flag <- TRUE
#     parsed_numbers <- -1
#     counter <- 1
#     while(flag){
#         right_most_digit <- x%%10
#         parsed_numbers <- c(parsed_numbers, right_most_digit)
#         temp <- x/10
#         x <- temp - (right_most_digit*0.1)
#         if(x == 0){
#             flag <- FALSE
#         }
#     }
#     return(rev(parsed_numbers[-1]))
# }



# count the number of each rating type/category in the rater's provided rating for one object
count_rating_types <- function(rating, types){
    counts <- rep(-1, length(types))
    for(i in 1:length(types)){
        counts[i] <- length(which(rating == types[i]))
    }
    return(counts)
}

count_rating_types_perrater <- function(ratings, rating_categories, multiple_per_rater){
    num_raters <- ncol(ratings)
    num_categories <- length(rating_categories)
    rater_object_count_types <- list()
    num_objects <- nrow(ratings)
    
    for(i in 1:num_raters){
        rater_responses <- ratings[, i]
        count_matrix <- matrix(-1, num_objects, num_categories)
        for(j in 1:num_objects){
            if(multiple_per_rater){
                count_matrix[j, ] <- t(count_rating_types(parse_concat_ratings(rater_responses[j]), rating_categories))
            }else{
                count_matrix[j, ] <- t(count_rating_types(rater_responses[j], rating_categories))
            }
        }
        rater_object_count_types[[i]] <- count_matrix
    }
    
    return(rater_object_count_types)
}

# estimation of the error rates of every rater see eq (2.3) in Dawid et al. 
estimate_error_rate <- function(truelabel_mat, raters_ratings){
    num_raters <- length(raters_ratings)
    #getting the dimension of the error matrix per rater
    num_col <- ncol(truelabel_mat)
    num_row <- num_col
    # list to hold the error rater matrix for every rater
    error_rate_list <- list()
    
    for(i in 1:num_raters){
        error_rate_mat <- matrix(-1, num_row, num_col)
        # current rater ratings matrix of dimension num_objects x num_rating_categories
        rater_mat <- raters_ratings[[i]]
        for(j in 1:num_row){
            current_truelabels <- truelabel_mat[, j]
            for(k in 1:num_col){
                count_rater_labeltype <- rater_mat[, k]
                error_rate_mat[j, k] <- sum(current_truelabels * count_rater_labeltype)
            }
            # normalize 
            error_rate_mat[j, ] <- error_rate_mat[j, ] / sum(error_rate_mat[j, ])
        }
        # save the error rate matrix of every rater in a list
        error_rate_list[[i]] <- error_rate_mat
    }
    return(error_rate_list)
    
}


# estimate the probability of a rating/label of type j where j \in [1,J]
estimate_pj <- function(truelabel_mat){
    p_j <- apply(truelabel_mat, 2, FUN = sum)
    return(p_j/sum(p_j)) 
}




# update the probability of T_{ij} (i.e. the probability of the true label for every object i)
# it uses the log domain to estimate the probabilites i.e. log(p(T_{ij})) in order not to risk underflow
update_truelabel <- function(error_rate_list, raters_ratings, prob_categ){
    
    # computing the log sum of probabilitis 
    log_sum_exp <- function(logjointdata){
        C_max <- max(logjointdata)
        return(C_max + log(sum(exp(logjointdata - C_max))))
    }
    
    num_objects <- nrow(raters_ratings[[1]])
    num_categories <- length(prob_categ)
    num_raters <- length(error_rate_list)
    truelabel_mat <- matrix(-1, num_objects, num_categories)
    # probability of the data/responses observed for single object
    p_data_object <- rep(0, num_objects)
    
    for(i in 1:num_objects){
        for(j in 1:num_categories){
            a_j <- 0
            for(k in 1:num_raters){
                rater_err <- error_rate_list[[k]][j, ]
                rater_label_counttypes <- raters_ratings[[k]][i, ]
                index_nonzero_count <- which(rater_label_counttypes != 0)
                a_j <- a_j + sum(log(rater_err[index_nonzero_count]) * rater_label_counttypes[index_nonzero_count])
                
            }
            truelabel_mat[i, j] <- a_j + log(prob_categ[j])
        }
        p_data_object[i] <- log_sum_exp(truelabel_mat[i, ])
        truelabel_mat[i, ] <- exp(truelabel_mat[i, ] - p_data_object[i])
    }
    return(list(truelabel_mat, p_data_object))
    
}

# stop iterating if the difference in the likelihood of the data from iteration t and t+1 is below the threshold
check_convergence_likelihood <- function(prev_likelihood, curr_likelihood, conv_threshold){
    
    diff <- abs(prev_likelihood - curr_likelihood)
    if (diff > conv_threshold){
        return(FALSE)
    }
    
    return(TRUE)
}

# get the index of the parameters that have the maximum log likelihood across all runs
get_max_params_index <- function(max_params_per_run){
    loglikelihood <- rep(-1,length(max_params_per_run))
    for(i in 1:length(max_params_per_run)){
        loglikelihood[i] <- max_params_per_run[[i]][[4]]
    }
    max_index <- which(max(loglikelihood) == loglikelihood)
    return(max_index)
} 

# check if the dependency package is installed
# it is mainly used to install gsubfn package to parse multiple ratings from a single rater for the same object
check_package_dependency <- function(required_packages){
    cond <- required_packages %in% installed.packages()[, "Package"]
    pack_toinstall <- required_packages[!cond]
    if(length(pack_toinstall)){
        install.packages(pack_toinstall)
        load(pack_toinstall)
    }
}

find_consensus_rating <- function(data, rating_categories, num_runs, num_iter, convergence_threshold, multiple_per_rater = FALSE){
    # install gsubfn package in case multiple ratings from a single rater for the same object is provided
    if(multiple_per_rater){
        check_package_dependency(c("gsubfn"))
    }
    raters_ratings <- count_rating_types_perrater(data, rating_categories, multiple_per_rater)
    num_raters <- ncol(data)
    num_objects <- nrow(data)
    max_params_per_run <- list()   
    loglikelihood_diagnostic <- list()
    
    for(i in 1:num_runs){
        print(i)
        counter <- 0
        converged <- FALSE
        converged_to_threshold <- FALSE
        # initialize the true label matrix using the ratings provided by raters
        # we use both methods; the first run uses the data to initialize and for subsequent
        # runs we use random initalization
        if(i == 1){
            truelabel_mat <- init_truelabel_majvote(raters_ratings)
        }else{
            truelabel_mat <- init_truelabel_random(num_objects, rating_categories)
        }
        
        while(!converged){
            counter <- counter + 1
            #       print(paste("we are in iteration number", counter, sep=" "))
            error_rate_list <- estimate_error_rate(truelabel_mat, raters_ratings)
            prob_label_types <- estimate_pj(truelabel_mat)
            result_list <- update_truelabel(error_rate_list, raters_ratings, prob_label_types)
            # matrix containing the probability distribution on the label types for every object
            truelabel_mat <- result_list[[1]]
            # vector containg the log likelihood of the data on every object
            # summing over all objects would result in the full log likelihood of the data given the current values
            # of the parameters
            loglikelihood <- sum(result_list[[2]])
            if(counter == 1){
                loglikelihood_diagnostic[[paste("run", i, sep = "_")]] <- loglikelihood
                track_params <- list(error_rate_list = error_rate_list,
                                     prob_label_types = prob_label_types,
                                     truelabel_mat = truelabel_mat,
                                     loglikelihood = loglikelihood,
                                     converged_to_threshold = converged_to_threshold)
                max_iter_params <- track_params
            }else{
                loglikelihood_diagnostic[[paste("run", i, sep = "_")]] <- c(loglikelihood_diagnostic[[paste("run", i, sep = "_")]],
                                                                            loglikelihood)
                converged <- check_convergence_likelihood(track_params[["loglikelihood"]], loglikelihood, convergence_threshold)
                
                if(converged){
                    converged_to_threshold <- TRUE
                }else if(counter == num_iter & !converged){
                    converged <- TRUE
                    converged_to_threshold <- FALSE
                }
                
                track_params <- list(error_rate_list = error_rate_list,
                                     prob_label_types = prob_label_types,
                                     truelabel_mat = truelabel_mat,
                                     loglikelihood = loglikelihood,
                                     converged_to_threshold = converged_to_threshold)
                # checking which parameters lead to higher log likelihood 
                if(max_iter_params[["loglikelihood"]] < track_params[["loglikelihood"]]){
                    max_iter_params <- track_params
                }
            } 
        }
        max_params_per_run[[paste("run", i, sep = "_")]] <- max_iter_params
    }
    
    max_index <- get_max_params_index(max_params_per_run)    
    return(list(max_param_index = max_index, 
                max_params_per_run = max_params_per_run, 
                loglikelihood_diagnostic = loglikelihood_diagnostic))
}