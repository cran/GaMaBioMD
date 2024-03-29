% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compiled_code_gamabiomd.R
\name{tree_average_similarity}
\alias{tree_average_similarity}
\title{Generate Phylogenetic Tree and Color Palette for Average Similarity}
\usage{
tree_average_similarity(similarity_matrix)
}
\arguments{
\item{similarity_matrix}{A matrix containing pairwise similarities between samples.}
}
\value{
A list containing the phylogenetic tree and a color palette.
}
\description{
This function generates a Neighbor-Joining phylogenetic tree and a color palette based on the average similarity matrix.
}
\examples{
\donttest{
accession_ranges <- list(
  SRU1 = "AJ240966 to AJ240970",
  STU2 = "AB015240 to AB015245",
  WPU13 = "L11934 to L11939",
  INU20 = c("AF277467 to AF277470", "AF333080 to AF333085")
)

# Use the function to expand accession ranges
sam_acc <- expand_accession_ranges(accession_ranges)
print(sam_acc)

# 2 get_sequence_information
accessions_to_query <- sam_acc$accession
seq_info <- get_sequence_information(accessions_to_query, remove_dot_1 = TRUE)
print(seq_info)
result <- preprocess_for_alignment(sam_acc, seq_info)

# Access the resulting data frames
merged_data <- result$merged_data
main_data <- result$main_data
final_data <- result$final_data

# If you want to sample 10\% from each SampleID group:
sampled_data <- data_sampling(final_data, sample_proportion = 0.1)

alignment_results <- alignment_info(final_data, type = "global", verbose = 1)

# Access the resulting data frames
score_matrix <- alignment_results$score_matrix
normalized_score_matrix <- alignment_results$normalized_score_matrix

total_aligned_positions_matrix <- alignment_results$total_aligned_positions_matrix
number_of_matching_positions_matrix <- alignment_results$number_of_matching_positions_matrix

percent_similarity_matrix <- alignment_results$percent_similarity_matrix

alignment_results_list <- alignment_results$alignment_results_list

alignment_info_matrix <- alignment_results$alignment_info_matrix

output_directory <- tempdir()

# Save the list of alignment results to an RDS file
saveRDS(alignment_results_list, file.path(output_directory, "alignment_results_list.rds"))

# Save matrices to files
write.table(score_matrix, file.path(output_directory, "score_matrix.txt"), sep = "\t")
average_percent_similarity <- compute_average_similarity_matrix(percent_similarity_matrix)
print(average_percent_similarity)

output_directory <- tempdir()
width_inch <- 8
height_inch <- 6
dpi <- 300

clustering_result <- clustering_average_similarity(average_percent_similarity)

# Extract the dendrogram and clustered data
dend_colored <- clustering_result$dendrogram
clustered_data <- clustering_result$clustered_data
Cluster_SampleID_Percentage <- clustering_result$Cluster_SampleID_Percentage
Cluster_TotalPercentage <- clustering_result$Cluster_TotalPercentage

tiff_file <- file.path(output_directory, "6. hierarchical_clustering_dendrogram_colored.tiff")

# Save the dendrogram as a TIFF image
tiff(tiff_file, width = width_inch, height = height_inch, units = "in", res = dpi)
plot(dend_colored, main = "Colored Hierarchical Clustering Dendrogram")
dev.off()

# Save the clustered data frame to a CSV file
write.csv(clustered_data, file.path(output_directory, "7. clustered_data.csv"), row.names = FALSE)

# Example usage with similarity_matrix
result <- tree_average_similarity(average_percent_similarity)
tree <- result$tree
color_palette <- result$color_palette

output_directory <- tempdir()

tree_newick <- file.path(output_directory, "phylogenetic_tree_nj.nwk")

# Save the phylogenetic tree as a Newick file
# ape::write.tree(tree, file = tree_newick)

# Save the phylogenetic tree as a TIFF image
width_inch <- 8
height_inch <- 6
dpi <- 300

# Set the file name for the TIFF image
tiff_file <- file.path(output_directory, "phylogenetic_tree.tiff")

# Open the TIFF device
tiff(tiff_file, width = width_inch, height = height_inch, units = "in", res = dpi)

# Plot the phylogenetic tree vertically
plot(tree, main = "Neighbor-Joining Tree", cex = 1, direction = "downward")

# Close the TIFF device
dev.off()

# Set the file name for the TIFF image with rainbow-colored branches
tiff_file_rainbow <- file.path(output_directory, "phylogenetic_tree_rainbow.tiff")

# Open the TIFF device
tiff(tiff_file_rainbow, width = width_inch, height = height_inch, units = "in", res = dpi)

# Plot the phylogenetic tree vertically with rainbow-colored branches
plot(tree, main = "NJ Tree", cex = 1, direction = "downward", tip.color = color_palette)

# Close the TIFF device
dev.off()
}
}
