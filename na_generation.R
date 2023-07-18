na_generation <- function(x, na_rate = 0.1){
  if (is.matrix(x) == FALSE){
    stop("wrong data type for x")
  }
  total <- dim(x)[1]*dim(x)[2]
  number_na <- total*na_rate
  indexes <- sample(total, number_na)
  x[indexes] <- NA
  return(x)
}


na_generation_v2 <- function(matrix, na_matrix, na_rate){
  if (is.matrix(matrix) == FALSE | is.matrix(na_matrix) == FALSE){
    stop("wrong data type for x")
  }
  
  total <- dim(matrix)[1]*dim(matrix)[2]
  number_na <- total*na_rate
  na_pos <- which(na_matrix)
  val_pos <- which(na_matrix == FALSE)
  
  if(length(na_pos) >= number_na){
    indexes <- sample(na_pos, number_na)
  }
  
  else{
    val_num<- number_na - length(na_pos)
    new_na <- sample(val_pos, val_num)
    indexes <- c(na_pos, new_na)  
  }
  matrix[indexes] <- NA
  return(matrix)
}

