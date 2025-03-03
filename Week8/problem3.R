library(tidyverse)

# Read the results from both models
model1_results <- read_tsv("model1_results.tsv", col_types = cols(chr = col_character(), pos = col_double(), log10p = col_double()))
model2_results <- read_tsv("model2_results.tsv", col_types = cols(chr = col_character(), pos = col_double(), log10p = col_double()))

# Merge the results using full_join to retain all data
merged_results <- model1_results %>%
  rename(log10p_model1 = log10p) %>%
  full_join(model2_results %>% rename(log10p_model2 = log10p), 
            by = c("chr", "pos")) %>%
  # Replace NA values with 0 or another placeholder if needed
  mutate(across(starts_with("log10p_"), ~ replace_na(.x, 0)))  

# Write the merged results
write_tsv(merged_results, "merged_results.tsv")

# Exit
quit()
