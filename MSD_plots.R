library(dplyr)
library(tidyverse)
library(ggplot2)
library(magrittr)

ris_data <- read.csv("/media/varun/CRl_DATA/msd/plots/SPN Indian adult study/serotype_vs_resitance.csv")

# Function to calculate percentage
calculate_percentage <- function(x) {
  R_count <- sum(x == "R") / length(x) * 100
  I_count <- sum(x == "I") / length(x) * 100
  S_count <- sum(x == "S") / length(x) * 100
  return(c(R = R_count, I = I_count, S = S_count))
}

# Aggregate data by Serotype
aggregated_data <- aggregate(ris_data[, -1], by = list(Serotype = ris_data$Serotype), FUN = calculate_percentage)

# Reorder columns
aggregated_data <- aggregated_data[, c("Serotype", "Penicillin_RIS", "Cefotaxime_RIS", "Ceftriaxone_RIS", "Levofloxacin_RIS", "Moxifloxacin_RIS", "Erythromycin_RIS", "Clindamycin_RIS", "Vancomycin_RIS", "Tetracycline_RIS", "Chloramphenicol_RIS", "COT_RIS")]


#write the file
#write.csv(aggregated_data,"/media/varun/CRl_DATA/msd/per_ris_serotype.csv")
##there are spaces after aggregating data. That couldn't be resolved in R hence write it as a csv replace the spaces and clean it and rload the data

#read cleaned_data
aggr_data <- read.csv("/media/varun/CRl_DATA/msd/plots/SPN Indian adult study/per_ris_serotype_only_R.csv")

#convert data to long
aggregated_data_long <- pivot_longer(aggr_data, -Serotype, names_to = "Antibiotic", values_to = "Percentage")
#filtered_data <- aggregated_data_long[grep("_R", aggregated_data_long$Antibiotic), ]
#filtered_data$Antibiotic <- gsub("_R", "", filtered_data$Antibiotic)

aggregated_data_long$Serotype <- as.character(aggregated_data_long$Serotype)


# Define the desired order of serotypes
desired_order <- read.csv("/media/varun/CRl_DATA/msd/plots/SPN Indian adult study/serotype_order.csv", header = FALSE)
desired_order_serotype <- desired_order$V1
print(desired_order_serotype)

aggregated_data_long$Serotype <- factor(aggregated_data_long$Serotype, levels = desired_order_serotype)
aggregated_data_long <- aggregated_data_long %>% arrange(Serotype)



desired_order_antibiotic <- c("Penicillin", "Cefotaxime", "Ceftriaxone", "Levofloxacin", "Moxifloxacin",
                              "Erythromycin", "Clindamycin", "Vancomycin", "Tetracycline", "Chloramphenicol", "COT")


# Reorder Serotype and Antibiotic columns based on desired order
filtered_data <- aggregated_data_long %>%
  mutate(Serotype = factor(Serotype, levels = desired_order_serotype),
        Antibiotic = factor(Antibiotic, levels = desired_order_antibiotic)) 
#write.csv(filtered_data,"/media/varun/CRl_DATA/msd/R_vs_serotype_long.csv")


filtered_data_unique <- read.csv("/media/varun/CRl_DATA/msd/R_vs_serotype_long.csv")
#convert to wide format
wide_data <- pivot_wider(filtered_data_unique, names_from = Antibiotic, values_from = Percentage)

#write to csv
write.csv(wide_data,"/media/varun/CRl_DATA/msd/R_vs_serotype_ordered.csv")


# Plot the heatmaaggregated_data_long# Plot the heatmap
ggplot(filtered_data, aes(x = Antibiotic, y = Serotype, fill = Percentage)) +
  geom_tile() +
  scale_fill_gradient(low = "#ebf2fa", high = "#f65135") +
  geom_hline(yintercept = seq(0.5, length(desired_order_serotype) - 0.5), color = "grey") +
  geom_vline(xintercept = seq(0.5, length(desired_order_antibiotic) - 0.5), color = "grey") +
  labs(title = "Heatmap of Percentage Resistance",
       x = "Antibiotics",
       y = "Serotype",
       fill = "Value") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



vfdb_data <- read.csv("/media/varun/CRl_DATA/msd/plots/SPN Indian adult study/vfdb.csv")

# Calculate the percentage of each gene for "DISEASE"
percentage_genes_disease <- vfdb_data%>%
  filter(Invasiveness == "DISEASE") %>%
  summarise(across(starts_with("cba"):ends_with("zmpC"), ~ mean(. == "yes", na.rm = TRUE) * 100))

# Calculate the percentage of each gene for "CARRIAGE"
percentage_genes_carriage <- vfdb_data%>%
  filter(Invasiveness == "CARRIAGE") %>%
  summarise(across(starts_with("cba"):ends_with("zmpC"), ~ mean(. == "yes", na.rm = TRUE) * 100))

combined_percentage <- rbind(percentage_genes_disease, percentage_genes_carriage)

Invasiveness <- c("Disease", "carriage")
new_data <- cbind(Invasiveness, combined_percentage)
#write.csv(new_data,"/media/varun/CRl_DATA/msd/vfdb_per.csv")


data_long <- pivot_longer(new_data, -Invasiveness, names_to = "Gene")


# Plot
ggplot(data_long, aes(x = value, y = Gene, fill = Invasiveness)) +
  geom_bar(stat = "identity", position= "dodge", width = 1) +  # Adjust width
  geom_text(aes(label = paste0(round(value))), hjust = -0.5, size = 3, color = "black") +  # Add data labels
  scale_fill_manual(values = c("blue", "red")) +
  labs(title = "Invasiveness by Gene", x = "Percentage", y = "Gene") +
  theme_minimal() +
  theme(legend.title = element_blank())



##import ast data

ast <- read.csv("/media/varun/CRl_DATA/msd/plots/SPN Indian adult study/ast.csv")

# Columns ending with _MIC
mic_cols <- grep("_MIC$", names(ast), value = TRUE)

# Filter out rows with missing values in columns ending with _MIC
df <- ast[complete.cases(ast[mic_cols]), ]

# Loop through MIC columns and replace symbols
for (col in mic_cols) {
  df[[col]] <- as.numeric(gsub("[<>=]", "", df[[col]]))
}


# Convert MIC columns to numeric
df[mic_cols] <- lapply(df[mic_cols], as.numeric)
df[mic_cols] <- log(df[mic_cols])

# Plot violin plot
ggplot(df, aes(x = Invasiveness, y = Penicillin_MIC, fill = Invasiveness)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.2, fill = "#C5C9A4", color = "black") +
  labs(title = "Penicillin MIC Distribution for Disease and Carriage",
       x = "Invasiveness",
       y = "Penicillin MIC") +
  scale_fill_manual(values = c("DISEASE" = "red", "CARRIAGE" = "blue")) +
  theme_minimal()

# Plot violin plot
ggplot(df, aes(x = Invasiveness, y = Ceftriaxone_MIC, fill = Invasiveness)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "#C5C9A4", color = "black") +
  labs(title = "Ceftriaxone MIC Distribution for Disease and Carriage",
       x = "Invasiveness",
       y = "Ceftriaxone MIC") +
  scale_fill_manual(values = c("DISEASE" = "red", "CARRIAGE" = "blue")) +
  theme_minimal()

# Plot violin plot
ggplot(df, aes(x = Invasiveness, y = Cefotaxime_MIC, fill = Invasiveness)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "#C5C9A4", color = "black") +
  labs(title = "Cefotaxime MIC Distribution for Disease and Carriage",
       x = "Invasiveness",
       y = "Cefotaxime MIC") +
  scale_fill_manual(values = c("DISEASE" = "red", "CARRIAGE" = "blue")) +
  theme_minimal()


# Plot violin plot
ggplot(df, aes(x = Invasiveness, y = Erythromycin_MIC, fill = Invasiveness)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "#C5C9A4", color = "black") +
  labs(title = "Erythromycin MIC Distribution for Disease and Carriage",
       x = "Invasiveness",
       y = "Erythromycin MIC") +
  scale_fill_manual(values = c("DISEASE" = "red", "CARRIAGE" = "blue")) +
  theme_minimal()



# Plot violin plot
ggplot(df, aes(x = Invasiveness, y = COT_MIC, fill = Invasiveness)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "#C5C9A4", color = "black") +
  labs(title = "COT MIC Distribution for Disease and Carriage",
       x = "Invasiveness",
       y = "Erythromycin MIC") +
  scale_fill_manual(values = c("DISEASE" = "red", "CARRIAGE" = "blue")) +
  theme_minimal()




# Read the csv file
serotype_manifest <- read.csv("/media/varun/CRl_DATA/msd/plots/SPN Indian adult study/serotype_vs_disease.csv")
# Summarize the data
serotype_manifest_summary <- serotype_manifest %>%
  group_by(Serotype, manifest) %>%
  summarise(Count = n(), .groups = 'drop') 

serotype_manifest_summary$Serotype <- factor(serotype_manifest_summary$Serotype, levels = desired_order_serotype)
serotype_manifest_summary <- serotype_manifest_summary %>% arrange(Serotype)

# Plot the data for bar chart of each serotype coloured by invasive and non invasive isolates
ggplot(serotype_manifest_summary, aes(x = Serotype, y = Count, fill = manifest)) +
  geom_bar(stat = "identity", position = "stack") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate the x-axis labels to avoid overlap
  labs(title = "Counts of Disease and Carriage by Serotype", x = "Serotype", y = "Count") +
  scale_fill_brewer(palette = "Set1") # Use a color palette that is easily distinguishable

