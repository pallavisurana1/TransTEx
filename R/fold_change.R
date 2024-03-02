#' Calculate Fold Change for Tissue-Specific
#'
#' @param data A data frame containing TPM values, sample classes, and other relevant metadata.
#' @param category A character string specifying the analysis category: 'tissue-specific' or 'tissue-enhanced'.
#' @param tissue_interest A character vector of tissue names of interest for tissue-specific calculations.
#' @param df2 A data frame containing PI tissue groupings for tissue-enhanced calculations.
#' @return A data frame with log2 fold change calculations.
#' @export
fold_change_TSp <- function(data, category = "tissue-specific", tissue_interest) {
  # Filter and prepare data
  data <- data %>%
    select(starts_with("ENST") | any_of(c("tissue", "sample")))

    # mean TPM in tissue of interest
    mean_tpm_df <- data %>%
      filter(tissue %in% tissue_interest) %>%
      group_by(tissue) %>%
      summarise(across(starts_with("ENST"), mean, na.rm = TRUE))

    # mean TPM not in tissue of interest
    not_mean_tpm_df <- bind_rows(lapply(unique(data$tissue), function(g) {
      df_not_in <- data %>%
        filter(tissue != g) %>%
        summarise(across(starts_with("ENST"), mean, na.rm = TRUE), .groups = 'drop') %>%
        mutate(tissue = g)
    }))

    # Calculate fold change
    fc <- log2(mean_tpm_df[, -which(names(mean_tpm_df) == "tissue")] / not_mean_tpm_df[, -which(names(not_mean_tpm_df) == "tissue")])
    return(data.frame(tissue = mean_tpm_df$tissue, fc))
}

##--------------------------------------------------------------------------------------------------------

#' Calculate Fold Change for Tissue-Enhanced Transcripts
#'
#' @param df1 Data frame with mean TPM values (mTPM) for each sample class, and additional metadata.
#' @param df2 Data frame with PI tissue groupings for tissue-enhanced calculations.
#' @return A modified data frame and log2 fold change calculations.
#' @export
fold_change_TEn <- function(df1, df2) {
  # Ensure column names in df1 start with "ENST"
  transcript_cols <- names(df1)[grepl("^ENST", names(df1))]

  # Initialize vectors to store results
  combined_means_numerator <- numeric(length(transcript_cols))
  combined_means_denominator <- numeric(length(transcript_cols))

  # Process each transcript column
  for (transcript_col in transcript_cols) {
    df2 <- df2 %>%
      rowwise() %>%
      mutate(
        numerator = sum(df1[match(unlist(strsplit(tissue, ", ")), df1$tissue), transcript_col] * df1$num[match(unlist(strsplit(tissue, ", ")), df1$tissue)], na.rm = TRUE) / sum(df1$num[match(unlist(strsplit(tissue, ", ")), df1$tissue)], na.rm = TRUE),
        denominator = sum(df1[which(!df1$tissue %in% unlist(strsplit(tissue, ", "))), transcript_col] * df1$num[which(!df1$tissue %in% unlist(strsplit(tissue, ", ")))], na.rm = TRUE) / sum(df1$num[which(!df1$tissue %in% unlist(strsplit(tissue, ", ")))], na.rm = TRUE)
      ) %>%
      ungroup()

    combined_means_numerator <- c(combined_means_numerator, df2$numerator)
    combined_means_denominator <- c(combined_means_denominator, df2$denominator)
  }

  # Calculate log2 fold change
  df2$log2foldchange <- log2(df2$numerator / df2$denominator)

  # Drop intermediate columns if necessary
  df2 <- df2 %>% select(-numerator, -denominator)

  return(df2)
}

##--------------------------------------------------------------------------------------------------------

#' Calculate Fold Change Based on Maximum Value for Tissue-Specific Transcripts
#'
#' This function calculates the fold change for each transcript in each tissue
#' compared to its maximum expression across the rest of tissues tissues, then applies a log2 transformation.
#'
#' @param data A matrix with tissues, samples and transcript/genes as columns (expression matrix)
#' @return A data frame with tissues and log2-transformed fold changes.
#' @export
fold_change_max_TSp <- function(data) {

  # mean TPM in tissue of interest
  data1 <- data %>%
    group_by(tissue) %>%
    summarise(across(starts_with("ENST"), mean, na.rm = TRUE)) %>% as.data.frame
  data1 = data1 %>% column_to_rownames("tissue")

  # Initialize a matrix to store the log2-transformed fold change results
  result <- matrix(0, nrow = nrow(data1), ncol = ncol(data1))
  rownames(result) <- rownames(data1)
  colnames(result) <- colnames(data1)

  # Iterate through each cell in the input matrix to calculate fold change
  for (i in 1:nrow(data1)) {
    for (j in 1:ncol(data1)) {
      # Current value
      val = data1[i, j]
      # Maximum value in the column, excluding the current row
      max_val = max(data1[-i, j], na.rm = TRUE)
      # Calculate fold change, apply log2 transformation, and handle division by zero or max_val being zero
      result[i, j] = ifelse(max_val > 0, log2(val / max_val + 1e-9), NA) # Adding a small value to avoid log2(0)
    }
  }

  return(result)
}

##--------------------------------------------------------------------------------------------------------

#' Calculate Fold Change for Tissue-Enhanced Transcripts Based on Maximum Value Outside Tissue Group
#'
#' This function computes the mean TPM for specified groups of tissues, compares it to the maximum mean TPM
#' outside those groups, and calculates the log2 fold change max.
#'
#' @param df1 A data frame with TPM values, sample class, and possibly other metadata.
#'            It should include a 'num' column for weighting in mean calculations.
#' @param df2 A data frame with PI tissue groupings for tissue-enhanced calculations.
#' @return A modified df2 with additional columns for combined mean, maximum outside mean, and log2 fold change.
#' @export
TenhTrans_FC_max <- function(df1, df2) {
  # Ensure df1 has 'tissue' and 'num' columns
  if (!("tissue" %in% names(df1)) | !("num" %in% names(df1))) {
    stop("df1 must contain 'tissue' and 'num' columns.")
  }

  # Calculate mean TPM within and maximum TPM outside tissue groups
  df2 <- df2 %>%
    rowwise() %>%
    mutate(
      tissues = unlist(strsplit(as.character(tissue), split = ",", fixed = TRUE)),
      combined_mean = list(df1 %>%
                             filter(tissue %in% tissues) %>%
                             summarise(across(starts_with("ENST"), ~sum(.x * num, na.rm = TRUE) / sum(num, na.rm = TRUE))) %>%
                             pivot_longer(everything(), names_to = "transcript", values_to = "mean_val") %>%
                             arrange(desc(mean_val)) %>%
                             slice_head(n = 1)),
      max_outside_mean = list(df1 %>%
                                filter(!(tissue %in% tissues)) %>%
                                summarise(across(starts_with("ENST"), max, na.rm = TRUE)) %>%
                                pivot_longer(everything(), names_to = "transcript", values_to = "max_val") %>%
                                arrange(desc(max_val)) %>%
                                slice_head(n = 1))
    ) %>%
    ungroup() %>%
    unnest(c(combined_mean, max_outside_mean))

  # Calculate log2 fold change
  df2$log2foldchange <- log2(df2$mean_val / df2$max_val)

  # Cleanup before return
  df2 <- df2 %>% select(-c(mean_val, max_val, tissues))

  return(df2)
}

##--------------------------------------------------------------------------------------------------------
