out_file_path = "/Data/IQRray_results/huExon.tsv"

if (!file.exists(out_file_path)) {
  IQRay_file <- NULL
  for (gseID in huExon$gseID) {
  final_score <- run_IQRray(gseID, TRUE)
    if (is.null(IQRay_file)) {
      IQRay_file <- final_score
    } else {
      IQRay_file <- rbind(IQRay_file, final_score)
    }
  }

  write_tsv(IQRay_file, out_file_path)
  print(paste0("Saved to ", out_file_path))
}

out_file_path = "/Data/IQRray_results/huGene.tsv"

if (!file.exists(out_file_path)) {
  IQRay_file <- NULL
  for (gseID in huGene$gseID) {
    final_score <- run_IQRray(gseID, TRUE)
    if (is.null(IQRay_file)) {
      IQRay_file <- final_score
    } else {
      IQRay_file <- rbind(IQRay_file, final_score)
    }
  }

  write_tsv(IQRay_file, out_file_path)
  print(paste0("Saved to ", out_file_path))
}

out_file_path = "/Data/IQRray_results/U95_2.tsv"

if (!file.exists(out_file_path)) {
  IQRay_file <- NULL
  for (gseID in U95_2$gseID) {
  final_score <- run_IQRray(gseID, FALSE)
    if (is.null(IQRay_file)) {
      IQRay_file <- final_score
    } else {
      IQRay_file <- rbind(IQRay_file, final_score)
    }
  }

  write_tsv(IQRay_file, out_file_path)
  print(paste0("Saved to ", out_file_path))
}

out_file_path = "/Data/IQRray_results/U133A_Early_Access.tsv"

if (!file.exists(out_file_path)) {
  IQRay_file <- NULL
  for (gseID in U133_A_Early_Access$gseID) {
  final_score <- run_IQRray(gseID, FALSE)
    if (is.null(IQRay_file)) {
      IQRay_file <- final_score
    } else {
      IQRay_file <- rbind(IQRay_file, final_score)
    }
  }

  write_tsv(IQRay_file, out_file_path)
  print(paste0("Saved to ", out_file_path))
}

out_file_path = "/Data/IQRray_results/U133_A.tsv"

if (!file.exists(out_file_path)) {
  IQRay_file <- NULL
  for (gseID in U133_A$gseID) {
  final_score <- run_IQRray(gseID, FALSE)
    if (is.null(IQRay_file)) {
      IQRay_file <- final_score
    } else {
      IQRay_file <- rbind(IQRay_file, final_score)
    }
  }

  write_tsv(IQRay_file, out_file_path)
  print(paste0("Saved to ", out_file_path))
}

out_file_path = "/Data/IQRray_results/U133_A2.tsv"

if (!file.exists(out_file_path)) {
  IQRay_file <- NULL
  for (gseID in U133_A2$gseID) {
  final_score <- run_IQRray(gseID, FALSE)
    if (is.null(IQRay_file)) {
      IQRay_file <- final_score
    } else {
      IQRay_file <- rbind(IQRay_file, final_score)
    }
  }

  write_tsv(IQRay_file, out_file_path)
  print(paste0("Saved to ", out_file_path))
}

out_file_path = "/Data/IQRray_results/U133_plus_2.tsv"

if (!file.exists(out_file_path)) {
  IQRay_file <- NULL
  for (gseID in U133_plus_2$gseID) {
  final_score <- run_IQRray(gseID, FALSE)
    if (is.null(IQRay_file)) {
      IQRay_file <- final_score
    } else {
      IQRay_file <- rbind(IQRay_file, final_score)
    }
  }

  write_tsv(IQRay_file, out_file_path)
  print(paste0("Saved to ", out_file_path))
}
