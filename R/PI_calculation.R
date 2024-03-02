#' Calculate Proportion Index (PI) Matrix
#'
#' This function calculates the Proportion Index matrix by first determining the count of values greater
#' than or equal to a threshold (iter) for each transcript (ENST columns) grouped by sample class,
#' then merging with tissue data and converting counts to proportions.
#'
#' @param data A data frame containing sample and ENST columns with values to be tested against iter.
#' @param iter A numeric threshold for determining whether the value in ENST columns is counted (default is 0.5).
#' @return A data frame with proportions for each tissue and transcript (ENST columns).
#' @export
PI_matrix <- function(data, iter = 0.5) {
  if (!all(c("sample", "tissue") %in% names(data))) {
    stop("data must contain 'sample' and 'tissue' columns")
  }
  if (!is.numeric(iter)) {
    stop("'iter' must be numeric")
  }

  # Get number of tissues
  tiss_df = data$tissue %>% table %>% as.data.frame; colnames(tiss_df) = c("tissue", "num")

  # Summarize data to find counts greater than or equal to 'iter'
  greater_df <- data %>%
    group_by(tissue) %>%
    summarize(across(starts_with("ENST"), ~sum(. >= iter, na.rm = TRUE)), .groups = "drop")

  # Merge with tiss_df_1 and calculate proportions
  prop_df <- merge(tiss_df, greater_df, by = "tissue")
  prop_df <- prop_df %>%
    mutate(across(starts_with("ENST"), ~ round(./num , 2)))

  return(prop_df)
}

##--------------------------------------------------------------------------------------------------------

#' Calculate Inflection Point (IPi) for PI matrix
#'
#' This function calculates the Inflection point which is a point on the PDF where the double-differential = 0.
#' than or equal to a threshold (iter) for each transcript (ENST columns) grouped by sample class,
#' then merging with tissue data and converting counts to proportions.
#'
#' @param data A data frame containing sample_class and ENST columns with output from function 'PI_matrix'.
#' @param threshold A numeric threshold for determining whether the value in ENST columns is counted (default is 0.5).
#' @return A data frame with proportions for each tissue and transcript (ENST columns).
#' @export
calculate_inflection_points <- function(data, threshold = 0.5) {
  # Ensure the tissue column is not in the data frame
  data$num <- NULL

  # Initialize the cutoff data frame
  cutoff <- data.frame(tissue = rownames(data), value = rep(NA, nrow(data)))

  for(i in seq_len(nrow(data))) {
    # Extract numeric values for the current tissue
    df1 <- unlist(data[i, ], use.names = FALSE)
    df1 <- df1[df1 > threshold]

    # Calculate the density and find inflection points
    if (length(df1) > 1) {
      df_density <- density(df1)
      x <- df_density$x
      y <- df_density$y

      # Compute first and second derivative
      xprime <- head(x, -1) + diff(x) / 2
      yprime <- diff(y) / diff(x)
      ydoubleprime <- diff(yprime) / diff(xprime)

      # Identify inflection points
      inflection_points <- which(diff(sign(ydoubleprime)) != 0)

      if (length(inflection_points) > 0) {
        # Find last positive inflection point less than 1
        positive_inflections <- subset(data.frame(x = xprime[inflection_points],
                                                  y = yprime[inflection_points]),
                                       x > 0 & y > 0 & x < 1)
        if (nrow(positive_inflections) > 0) {
          cutoff$value[i] <- tail(positive_inflections$x, 1)
        }
      }
    }
  }
  cutoff$value = round(cutoff$value, 2)
  return(cutoff)
}

##--------------------------------------------------------------------------------------------------------

#' Calculate Inflection Point (IPi) helper function to determine transcripts above IPi
#' @param data_prop  PI matrix.
#' @param up_lim  cutoff threshold.
#' @return A list containing the filtered values above threshold.
#' @export
prop_greater = function(data_prop, up_lim){
  df = data_prop
  df = df %>% t() %>% as.data.frame

  df_filtered <- df[apply(df, 1, function(x) any(x >= up_lim)), ]

  ##-----

  df_filtered = df_filtered %>% as.data.frame %>% t() %>% as.data.frame %>% rownames_to_column("tissue")
  filtered_df_lis = split(df_filtered, f = df_filtered$tissue)

  for (i in 1:length(filtered_df_lis)){
    rownames(filtered_df_lis[[i]]) = NULL
    filtered_df_lis[[i]] = filtered_df_lis[[i]] %>% as.data.frame %>% column_to_rownames("tissue") %>% as.data.frame %>% t()
    filtered_df_lis[[i]] = filtered_df_lis[[i]] %>% as.data.frame %>% rownames_to_column("transcript")
    colnames(filtered_df_lis[[i]]) = c("transcript", "prop")
    filtered_df_lis[[i]] = filtered_df_lis[[i]] %>% filter(prop > up_lim)
    filtered_df_lis[[i]]$prop %>% range
  }
  return(filtered_df_lis)
}

##--------------------------------------------------------------------------------------------------------

#' Process Data into Expression Classes
#'
#' This function classifies transcripts into different expression categories such as tissue-specific,
#' tissue-enhanced, widespread based on their expression in a number of tissues.
#'
#' @param data  A list which is the output after applying Inflection point cutoff.
#' @return A list containing the categorized data frames and a count table.
#' @export
divide_transtex_classes <- function(data) {
  # Check if columns needed are present
  required_cols <- c("tissue", "transcript", "prop")
  if (!all(required_cols %in% colnames(data))) {
    missing_cols <- required_cols[!required_cols %in% colnames(data)]
    stop(sprintf("Data must contain '%s' columns. Missing columns: %s",
                 paste(required_cols, collapse = "', '"),
                 paste(missing_cols, collapse = ", ")))
  }

  # Summarize data by transcript
  count_df <- data %>%
    group_by(transcript) %>%
    summarize(
      tissue = paste(tissue, collapse = ", "),
      count = n(),
      mean_prop = mean(prop, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate_all(as.character) %>%
    distinct() %>%
    mutate(count = as.integer(count))

  # Create a table with count frequencies
  count_table <- count_df$count %>% table() %>% as.data.frame()

  # Calculate the total number of tissues
  total_tissues <- length(unique(data$tissue))

  # Categorize transcripts based on count
  transtex <- list(
    TspTrans = count_df %>% filter(count == 1),
    TenhTrans = count_df %>% filter(count >= 2 & count <= floor(0.5 * total_tissues)),
    WideTrans = count_df %>% filter(count > floor(0.5 * total_tissues) & count <= total_tissues)
  )

  # Compile results
  result <- list(count_df = count_df, count_table = count_table, transtex = transtex)
  return(result)
}

##--------------------------------------------------------------------------------------------------------
