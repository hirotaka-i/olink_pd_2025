# Load necessary libraries
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(skimr))
library(optparse)

# Define the command line options with adjustment for multiple contrasts
option_list <- list(
  make_option(c("-d", "--data"), type = "character", default = NULL, help = "Path to the CSV data file", metavar = "FILE"),
  make_option(c("-p", "--protein_start"), type = "integer", default = NULL, help = "Starting column index for protein data"),
  make_option(c("-e", "--explanatory_vars"), type = "character", default = NULL, help = "Comma-separated list of explanatory variable names"),
  make_option(c("-c", "--contrast"), type = "character", default = NULL, help = "Contrasts for the analysis, separated by semicolons"),
  make_option(c("-l", "--log2_transform"), type = "logical", default = TRUE, help = "Apply log2 transformation (TRUE or FALSE)", metavar = "BOOLEAN")
)

# Parse the command line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Split the contrasts into a list if multiple contrasts are provided
contrast_list <- strsplit(opt$contrast, ";")[[1]]

# Function to perform the analysis
perform_DEA <- function(data_path, protein_start_col, explanatory_vars_string, contrast, log2_transform) {
  
    data <- read.csv(data_path, header = TRUE)
    # Split the explanatory_vars_string into a vector
    explanatory_vars <- unlist(strsplit(explanatory_vars_string, ","))

    # Check if explanatory variables are in the dataset
    missing_vars <- explanatory_vars[!explanatory_vars %in% colnames(data)]
    if (length(missing_vars) > 0) {
        stop("Error: The following explanatory variables are not in the dataset: ", paste(missing_vars, collapse = ", "))
    }

    # Extracting protein data and explanatory variables
    p1_prot <- as.matrix(data[, protein_start_col:ncol(data)])
    
    # Print number of proteins and the first and last 3 proteins
    print(paste("Number of proteins:", ncol(p1_prot), '[ FROM', head(colnames(p1_prot), 1), 'TO', tail(colnames(p1_prot), 1), ']'))

    # p_data <- data[,  c(names(data)[1], explanatory_vars)]
    p_data <- data[,  explanatory_vars]
    print(skim(p_data)) 

    # log2 transform (default==TRUE)
    if (log2_transform) {
        p1_prot <- log2(p1_prot + 1)  # Adding 1 to avoid log(0)
    }

    # remove the 3SD outliers
    p1_prot <- apply(p1_prot, 2, function(x) {
        if(all(is.na(x))) {
            return(rep(NA, length(x)))
        } else {
            return(ifelse(x > mean(x, na.rm = TRUE) + 3*sd(x, na.rm = TRUE) | x < mean(x, na.rm = TRUE) - 3*sd(x, na.rm = TRUE), NA, x))
        }
    })

    # Count the number of non-missing samples for each protein
    sample_counts_per_protein <- colSums(!is.na(p1_prot))
    names(sample_counts_per_protein) <- colnames(p1_prot)

    # Creating an ExpressionSet with data
    eset <- ExpressionSet(t(p1_prot), AnnotatedDataFrame(p_data))

    # Rest of the analysis
    design_formula <- as.formula(paste("~ 0 +", paste(colnames(p_data), collapse = " + ")))
    design <- model.matrix(design_formula, data = pData(eset))


    # Iterate over each contrast in the list
    for (contrast in contrast_list) {
      cm <- makeContrasts(contrasts=contrast, levels = design)
      fit <- lmFit(eset, design, maxit = 10000)
      fit2 <- contrasts.fit(fit, contrasts = cm)
      fit2 <- eBayes(fit2)
      stats <- topTable(fit2, number = nrow(fit2), sort.by = "P", adjust.method = "BH", confint = TRUE)
      
      # Add a column 'N' to stats for the number of non-missing samples per protein
      stats$N <- sample_counts_per_protein[row.names(stats)]
      result <- cbind(Probe = row.names(stats), stats)
    
      # Construct an output filename that includes the contrast to differentiate outputs
      sanitized_contrast <- gsub(" ", "_", gsub("-", "vs", contrast))  # Replace spaces and dashes for filename
      output_filename <- paste0("temp/", sanitized_contrast, ".csv")

      
      # Save the result with a unique filename for this contrast
      write.csv(result, file = output_filename, row.names = F)
      print(paste("Completed analysis for contrast:", contrast))
      print(paste("Number of significant proteins:", nrow(result[result$adj.P.Val < 0.05, ])))   

    }
}

perform_DEA(opt$data, opt$protein_start, opt$explanatory_vars, opt$contrast, opt$log2_transform)