
# function computing arIQR quality score for oligo arrays

IQRray_oligo <- function(data) {
  # obtaining intensity values for perfect match (pm) probes
  pm_data <- oligo::pm(data)

  # ranking probe intensities for every array
  rank_data <- apply(pm_data, 2, rank)

  # obtaining names of probeset for every probe
  probeNames <- oligo::probeNames(data)

  # function computing IQR of mean probe ranks in probesets
  get_IQR <- function(rank_data, probeNames) {
      round(IQR(sapply(split(rank_data, probeNames), mean)), digits = 2)
  }

  # computing arIQR score
  IQRray_score <- apply(rank_data, 2, get_IQR, probeNames = probeNames)

  return(IQRray_score)
}

# function computing arIQR quality score for affy arrays

IQRray_affy <- function(data) {
  # obtaining intensity values for perfect match (pm) probes
  pm_data <- affy::pm(data)

  # ranking probe intensities for every array
  rank_data <- apply(pm_data, 2, rank)

  # obtaining names of probeset for every probe
  probeNames <- affy::probeNames(data)

  # function computing IQR of mean probe ranks in probesets
  get_IQR <- function(rank_data, probeNames) {
      round(IQR(sapply(split(rank_data, probeNames), mean)), digits = 2)
  }

  # computing arIQR score
  IQRray_score <- apply(rank_data, 2, get_IQR, probeNames = probeNames)

  return(IQRray_score)
}
