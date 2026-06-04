
bind_file <- function(IQRay_file, final_score) {
  if (is.null(IQRay_file)) {
    IQRay_file <- final_score
  } else {
    IQRay_file <- rbind(IQRay_file, final_score)
  }

  return(IQRay_file)
}