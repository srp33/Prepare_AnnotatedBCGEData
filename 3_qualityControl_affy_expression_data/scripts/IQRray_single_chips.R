## script to run IQRray for single chips

# define the datasets with Gene ST and Exon ST arrays
oligo_arrays <- c("GSE33692", "GSE86374", "GSE58644", "GSE118432", "GSE59772", "GSE81838")

IQRay_file <- NULL
for (gseID in huExon$gseID) {
final_score <- run_IQRray(gseID)
  if (is.null(IQRay_file)) {
    IQRay_file <- final_score
  } else {
    IQRay_file <- rbind(IQRay_file, final_score)
  }
}
write_tsv(IQRay_file, paste0(IQRray_file_path, "huExon.tsv"))
print("Saved to huExon.tsv")

IQRay_file <- NULL
for (gseID in huGene$gseID) {
  final_score <- run_IQRray(gseID)
  if (is.null(IQRay_file)) {
    IQRay_file <- final_score
  } else {
    IQRay_file <- rbind(IQRay_file, final_score)
  }
}
write_tsv(IQRay_file, paste0(IQRray_file_path, "huGene.tsv"))
print("Saved to huGene.tsv")

IQRay_file <- NULL
for (gseID in U95_2$gseID) {
final_score <- run_IQRray(gseID)
  if (is.null(IQRay_file)) {
    IQRay_file <- final_score
  } else {
    IQRay_file <- rbind(IQRay_file, final_score)
  }
}
write_tsv(IQRay_file, paste0(IQRray_file_path, "U95_2.tsv"))
print("Saved to U95_2.tsv")

IQRay_file <- NULL
for (gseID in U133_A_Early_Access$gseID) {
final_score <- run_IQRray(gseID)
  if (is.null(IQRay_file)) {
    IQRay_file <- final_score
  } else {
    IQRay_file <- rbind(IQRay_file, final_score)
  }
}
write_tsv(IQRay_file, paste0(IQRray_file_path, "U133A_Early_Access.tsv"))
print("Saved to U133A_Early_Access.tsv")

IQRay_file <- NULL
for (gseID in U133_A$gseID) {
final_score <- run_IQRray(gseID)
  if (is.null(IQRay_file)) {
    IQRay_file <- final_score
  } else {
    IQRay_file <- rbind(IQRay_file, final_score)
  }
}
write_tsv(IQRay_file, paste0(IQRray_file_path, "U133_A.tsv"))
print("Saved to U133_A.tsv")

IQRay_file <- NULL
for (gseID in U133_A2$gseID) {
final_score <- run_IQRray(gseID)
  if (is.null(IQRay_file)) {
    IQRay_file <- final_score
  } else {
    IQRay_file <- rbind(IQRay_file, final_score)
  }
}
write_tsv(IQRay_file, paste0(IQRray_file_path, "U133_A2.tsv"))
print("Saved to U133_A2.tsv")

IQRay_file <- NULL
for (gseID in U133_plus_2$gseID) {
final_score <- run_IQRray(gseID)
  if (is.null(IQRay_file)) {
    IQRay_file <- final_score
  } else {
    IQRay_file <- rbind(IQRay_file, final_score)
  }
}
write_tsv(IQRay_file, paste0(IQRray_file_path, "U133_plus_2.tsv"))
print("Saved to U133_plus_2.tsv")
