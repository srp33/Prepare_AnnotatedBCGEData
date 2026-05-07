
summariseVariables <- function(metadata) {
  add_to_char <- NULL

  #change variables to factors
  meta <- metadata %>%
    mutate_all(type.convert, as.is = TRUE) %>%
    mutate(Dataset_ID = factor(Dataset_ID)) %>%
    mutate(Sample_ID = factor(Sample_ID))

  numMatrix <-  meta %>%
    select_if(negate(is.character))

  charMatrix <- meta %>%
    select_if(negate(is.numeric))

  #summarize numeric variables
  if ((ncol(numMatrix) > 2)) {
    numMatrix <- numMatrix %>%
        pivot_longer(3:ncol(numMatrix),
                     names_to = "Variable",
                     values_to = "Value") %>%
        group_by(Dataset_ID, Variable) %>%
        summarise(Min = min(Value, na.rm = T),
                  Mean = mean(Value, na.rm = T),
                  Max = max(Value, na.rm = T),
                  Num_Unique = length(unique(Value)),
                  Unique_Summary = paste(sort(unique(Value), na.last = TRUE), collapse = ", "))

    #remove columns with 4 or less unique values
    for (i in 1:nrow(numMatrix)) {
      bCol <- numMatrix[i, ]
      if (bCol$Num_Unique <= 4) {
        add_to_char <- (rbind(add_to_char, bCol))
      }
    }
    #check for number of rows
    if (!(is.null(add_to_char))) {
      numMatrix <- numMatrix %>%
        anti_join(add_to_char)
    }
  } else (numMatrix <- NULL)

  #summarize character variables
  if ((ncol(charMatrix) > 2)) {
    charMatrix <- charMatrix %>%
    pivot_longer(3:ncol(charMatrix),
                 names_to = "Variable",
                 values_to = "Value") %>%
    group_by(Dataset_ID, Variable) %>%
    summarise(Num_Unique = length(unique(Value)),
              Unique_Summary = paste(sort(unique(Value), na.last = TRUE), collapse = ", "))
  } else (charMatrix <- NULL)

  #combine character variables and numeric variables with less than 4 unique values
  if (!is.null(charMatrix)) {
    charMatrix <- (rbind(charMatrix, add_to_char)) %>%
      dplyr::select(Dataset_ID, Variable, Unique_Summary, Num_Unique)
  }

  return(list(numSummary = numMatrix, charSummary = charMatrix))
}
