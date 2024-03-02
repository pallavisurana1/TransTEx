#' Perform Stratified Sampling, Hypothesis Testing
#'
#' This function integrates the process of calculating sample sizes, performing stratified sampling and hypothesis testing,
#' and finally calculating empirical p-values for tissue-specific analysis.
#'
#' @param data A data frame containing the dataset for analysis with 'tissue' and 'ENST' prefixed columns for transcripts.
#' @param tissue_of_interest The tissue of interest for comparison.
#' @param size The desired sample size or proportion of the total, passed directly to `sample_size_calculator`.
#' @param iterations The number of iterations to perform the subsampling and testing.
#' @param threshold The p-value threshold used for empirical p-value calculation.
#' @return A data frame with transcript IDs, empirical p-values, and other metadata if specified.
#' @export
subsampling_empiricalP <- function(data, tissue_of_interest, size, iterations = 100, threshold = 0.05) {
  # Define transcript columns
  transcript_columns <- names(data)[startsWith(names(data), "ENST")]

  # Function to compute sample size
  sample_size_calculator <- function(population, size = NULL) {
    if (is.null(size)) {
      cochran_n <- round(((1.96)^2 * 0.5 * 0.5) / 0.02^2)
      n <- round(cochran_n / (1 + ((cochran_n - 1) / population)))
    } else if (size >= 0 && size < 1) {
      n <- round(population * size)
    } else if (size < 0) {
      stop("Parameter 'size' must be an integer or a proportion between 0 and 0.99.")
    } else if (size >= 1) {
      n <- size
    }
    return(n)
  }

  # Run stratified sampling and hypothesis testing
  stratified_sample_and_test <- function(data, tissue_of_interest, size, iterations) {
    sample_size <- sample_size_calculator(nrow(data), size)
    p_values <- matrix(NA, nrow = iterations, ncol = length(transcript_columns))

    ## run iter number of times or 100 = default
    for (i in 1:iterations) {

      # Perform stratified sampling
      sampled_data <- data %>% sample_n(size = sample_size, replace = FALSE)

      # Perform hypothesis testing for each transcript
      for (j in 1:length(transcript_columns)) {
        group1_data <- sampled_data[[transcript_columns[j]]][sampled_data$sample_class == tissue_of_interest]
        group2_data <- sampled_data[[transcript_columns[j]]][sampled_data$sample_class != tissue_of_interest]

        # Check if both groups have enough observations
        if (length(group1_data) >= 2 && length(group2_data) >= 2) {
          wilcox_1_smaller_result <- wilcox.test(
            group1_data,
            group2_data,
            alternative = "less"  # Perform one-tailed Wilcoxon test for smaller values
          )

          p_values[i, j] <- wilcox_1_smaller_result$p.value
        } else {
          # If either group has less than 2 observations, set p-value to NA
          p_values[i, j] <- NA
        }
      }
    }

    # Create a data frame with transcript names and p-values
    p_values_df <- data.frame(
      Transcript_ID = transcript_columns,
      P_Value = apply(p_values, 2, function(x) paste(x, collapse = ", "))
    )

    return(p_values_df)
  }

  p_values_df <- stratified_sample_and_test(data, tissue_of_interest, size, iterations)
  return(p_values_df)
}

#' Calculate Empirical P-values
#'
#' This function integrates the process of calculating sample sizes, performing stratified sampling and hypothesis testing,
#' and finally calculating empirical p-values for tissue-specific analysis.
#'
#' @param data A data frame containing the dataset for analysis with 'tissue' and 'ENST' prefixed columns for transcripts.
#' @param threshold The p-value threshold used to determine significance (default at 0.05).
#' @return A data frame with transcript IDs, empirical p-values, and other metadata if specified.
#' @export
calculate_empirical_p <- function(data, threshold = 0.05) {
  data$Empirical_P <- (rowSums(data$P_Value >= threshold, na.rm = TRUE) + 1) / (nrow(data) + 1)
  return(data)
}



