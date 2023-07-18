library(scp)
library(scpdata)
library(impute)
library(pcaMethods)
library(nipals)
library(missMDA)
library(microbenchmark)
library(MsCoreUtils)

start.time <- Sys.time()

################################################
                  #FUNCTIONS#
################################################


subset_generation <- function(dataset, row, col){
  #dataset should be a matrix
  original_row <-nrow(dataset)
  original_col <- ncol(dataset)
  features <- sample(original_row, row)
  samples <- sample(original_col, col)
  subset <- dataset[features, samples]
  return(subset)
}

subset_generation_wraper <- function(dim, dataset){
  col <- dim["cols"]
  row <- dim["rows"]
  return(subset_generation(dataset, row, col))
  }

pcaMethods_wraper <- function(matrix, method, npcs){
  matrix <- matrix[rowSums(is.na(matrix)) != ncol(matrix), ]
  matrix <- matrix[,colSums(is.na(matrix))<nrow(matrix)]
  pca <- pca(t(matrix), method = method, nPcs = npcs)
  return(pca)
}
svd_imputation_wraper <- function(matrix, npcs){
  matrix <- impute_knn(matrix, k = 3, colmax = 1)
  pca <- pcaMethods_wraper(matrix, method = "svd", npcs)
  return(pca)
}

nipals_wraper <- function(matrix, npcs){
  matrix <- matrix[rowSums(is.na(matrix)) != ncol(matrix), ]
  matrix <- matrix[,colSums(is.na(matrix))<nrow(matrix)]
  nipals <- nipals(t(matrix), maxiter = 5000, ncomp = npcs , startcol = 1, force.na = TRUE, scale =FALSE, center = FALSE)
  return(nipals)
}

# missMDA_wraper <- function(matrix){
#   comp <- estim_ncpPCA(matrix)
#   imputed_matrix <- imputePCA(matrix, ncp = ncomp$ncp)
#   miss_pca <- PCA(imputed_matrix$completeObs)
#   return(miss_pca)
# }

#extremely slow

benchmark_wraper <- function(npcs, matrix, imputed_matrix, times){
  bench <- microbenchmark(pcaMethods_wraper(matrix, "ppca", npcs),
                 pcaMethods_wraper(matrix, "nipals", npcs),
                 #pcaMethods_wraper(matrix, "bpca"),
                 pcaMethods_wraper(imputed_matrix, "svd", npcs),
                 svd_imputation_wraper(matrix, npcs),
                 nipals_wraper(matrix, npcs),
                 times = times)
  return(bench)
}

dataframe_maker <- function(bench, col, row){
  results <- as.data.frame(bench)
  dim <- data.frame("rows" = row, "cols"=col)
  merged <- cbind(results, dim)
  
  return(merged)
}

################################################
              # END OF FUNCTIONS#
################################################

leduc <- leduc2022()
leduc <- zeroIsNA(leduc, i = 1:138)
assay <- assay(leduc[["proteins_norm2"]])

leduc <- impute(leduc,
                  i = "proteins_norm2",
                  name = "prot_knn",
                  method = "knn",
                  k = 3, rowmax = 1, colmax= 1,
                  maxp = Inf
  )
imputed_assay <- assay(leduc[["prot_knn"]])

# the dimension of each subset size
#rows
rows <- c(300, 450, 600, 900, 1200, 1800, 2800)
cols <- c(200, 300, 450, 600, 900, 1200, 1500)

generated_subset <- mapply(function(rows, cols){
  subset_generation(dataset = assay, row = rows, col = cols)
  }, rows, cols)

generated_subset_imputed <- mapply(function(rows, cols){
  subset_generation(dataset = imputed_assay, row = rows, col = cols)
}, rows, cols)

# bench <- mapply(function(matrix, imputed_matrix){
#   benchmark_wraper(matrix = matrix, imputed_matrix = imputed_matrix, times = 10)
# }, generated_subset, generated_subset_imputed, SIMPLIFY = FALSE)

# Version to test the impact of Pcs number

bench_pc2 <- mapply(function(matrix, imputed_matrix){
    benchmark_wraper(npcs = 2, matrix = matrix, imputed_matrix = imputed_matrix, times = 10)
  }, generated_subset, generated_subset_imputed, SIMPLIFY = FALSE)

bench_pc15 <- mapply(function(matrix, imputed_matrix){
  benchmark_wraper(npcs = 10, matrix = matrix, imputed_matrix = imputed_matrix, times = 10)
}, generated_subset, generated_subset_imputed, SIMPLIFY = FALSE)

bench_pc50 <- mapply(function(matrix, imputed_matrix){
  benchmark_wraper(npcs = 30, matrix = matrix, imputed_matrix = imputed_matrix, times = 10)
}, generated_subset, generated_subset_imputed, SIMPLIFY = FALSE)

#bench <- lapply(generated_subset, benchmark_wraper, times = 2)
df_list_pc2 <- mapply(dataframe_maker, bench = bench_pc2, col = cols, row = rows)
df_list_pc15 <- mapply(dataframe_maker, bench = bench_pc15, col = cols, row = rows)
df_list_pc50 <- mapply(dataframe_maker, bench = bench_pc50, col = cols, row = rows)

df_pc2 <- as.data.frame(df_list_pc2[,1])
for (element in 2:ncol(df_list_pc2)){
  df_pc2 <- rbind(df_pc2, as.data.frame(df_list_pc2[,element]))
}
write.csv(df_pc2, file = "data_output/benchmark_run_pc2.csv")
df_pc15 <- as.data.frame(df_list_pc15[,1])
for (element in 2:ncol(df_list_pc15)){
  df_pc15 <- rbind(df_pc15, as.data.frame(df_list_pc15[,element]))
}
write.csv(df_pc15, file = "data_output/benchmark_run_pc10.csv")
df_pc50 <- as.data.frame(df_list_pc50[,1])
for (element in 2:ncol(df_list_pc50)){
  df_pc50 <- rbind(df_pc50, as.data.frame(df_list_pc50[,element]))
}
write.csv(df_pc50, file = "data_output/benchmark_run_pc30.csv")



end.time <- Sys.time()

time.taken <- end.time - start.time
time.taken