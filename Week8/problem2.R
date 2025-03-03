library(tidyverse)

# Define the function to fit the second model
Myfunction2 <- function(df) {
  # Ensure treat and founder are factors
  df <- df %>%
    mutate(
      treat = as.factor(treat),
      founder = as.factor(founder)
    )

  # Check if founder has at least two levels
  if (n_distinct(df$founder) < 2) {
    return(NA)  # Return NA if insufficient levels
  }

  # Try fitting the ANOVA model, handle errors gracefully
  tryCatch({
    model <- anova(lm(asin(sqrt(freq)) ~ founder + treat %in% founder, data = df))
    pval <- model$"Pr(>F)"[2]  # Extract p-value for the nested term
    return(-log10(pval))  # Return -log10(p)
  }, error = function(e) {
    return(NA)  # Return NA if an error occurs
  })
}

# Read the full dataset
mal <- read_tsv("allhaps.malathion.200kb.txt.gz")

# Add treatment column
mal <- mal %>% mutate(treat = str_sub(pool, 2, 2))

# Run the scan on the full dataset
results <- mal %>%
  group_by(chr, pos) %>%
  nest() %>%
  mutate(log10p = map_dbl(data, Myfunction2)) %>%
  select(chr, pos, log10p)

# Write results to a file
write_tsv(results, "model2_results.tsv")

# Exit
quit()
