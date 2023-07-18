library(scp)
library(scpdata)
library(impute)
library(pcaMethods)
library(nipals)
source("na_generation.R")
library(purrr)
library(sva)
library(limma)
library(tidyverse)

start.time <- Sys.time()



################################################
                  #FUNCTIONS#
################################################

na_wraping <- function(matrix, rate){
  na_matrix <- is.na(matrix)
  matrix <- matrix
  matrix <- na_generation_v2(matrix, na_matrix, rate)
  return(matrix)
}

na_wraping_random <- function(matrix, rate){
  na_matrix <- is.na(matrix)
  matrix <- matrix
  matrix <- na_generation(matrix, rate)
  return(matrix)
}

na_replication <- function(matrix, rates, rep){
  list_rep <- list()
  for (rep in 1:rep){
    list_rep[[length(list_rep)+1]] <- lapply(rates, FUN = na_wraping, matrix = matrix)
  }
  
  return(list_rep)
}

na_replication_random <- function(matrix, rates, rep){
  list_rep <- list()
  for (rep in 1:rep){
    list_rep[[length(list_rep)+1]] <- lapply(rates, FUN = na_wraping_random, matrix = matrix)
  }
  
  return(list_rep)
}

pcaMethods_wraper <- function(matrix, method){
  matrix <- matrix[rowSums(is.na(matrix)) != ncol(matrix), ]
  matrix <- matrix[,colSums(is.na(matrix))<nrow(matrix)]
  tryCatch(
    {
      pca <- pca(t(matrix), method = method, scale = "uv", center = TRUE)
      return(pca)
    },
    error = function(e){
      print("An error occurs while computing the pca with pcaMethods, replacing with NA")
      return(NA)
    }
  )
}

nipals_wraper <- function(matrix){ #add the correction
  matrix <- matrix[rowSums(is.na(matrix)) != ncol(matrix),]
  matrix <- matrix[,colSums(is.na(matrix))<nrow(matrix)]
  tryCatch(
    {
      nipals <- nipals(t(matrix), maxiter = 5000, ncomp = 2 , startcol = 1, force.na = TRUE, scale = TRUE, center = TRUE)
      nipals$scores <- nipals$scores %*% diag(nipals$eig)
      nipals$eigenvectors <- nipals$loadings
      nipals$eigenvalues <- nipals$eig^2 / (nrow(matrix) - 1)
      nipals$loadings <- nipals$loadings %*% diag(sqrt(nipals$eigenvalues))
      return(nipals)
    },
    error = function(e){
      print("An error occurs while computing the pca with nipals, replacing with NA")
      return(NA)
    }
  )
  
}

pca_na_run <- function(matrix){
  print(".")
  unique_results <- list()
  unique_results[[length(unique_results)+1]] <- pcaMethods_wraper(matrix, method = "nipals")
  unique_results[[length(unique_results)+1]] <- pcaMethods_wraper(matrix, method = "ppca")
  #unique_results[[length(unique_results)+1]] <- pcaMethods_wraper(matrix, method = "bpca")
  #unique_results[[length(unique_results)+1]] <- nipals_wraper(matrix)
  return(unique_results)
}

kmeans_wraper <- function(pca){
  print("-")
  tryCatch(
    {
      if (class(pca) == "pcaRes"){
        kmn <- kmeans(scores(pca),2,100)
        return(kmn)
      }
      if (class(pca) == "list"){
        kmn <- kmeans(pca$score, 2, 100)
        return(kmn)
      }
    },
    error = function(e){
      print("There was an error while computing the kmn, replacing with NA")
      return(NA)
    }
  )
  }

################################################
              # END OF FUNCTIONS#
################################################


################################################
              # PCA computing #
################################################

leduc <- leduc2022()

leduc <- zeroIsNA(leduc, i = 1:138)

leduc <- impute(leduc,
                i = "proteins_norm2",
                name = "prot_knn",
                method = "knn",
                k = 3, rowmax = 1, colmax= 1,
                maxp = Inf
)

prot_knn <- getWithColData(leduc, i = "prot_knn")
batch <- colData(prot_knn)$Set
model <- model.matrix(~SampleType, data = colData(prot_knn))
assay(prot_knn) <- ComBat(dat = assay(prot_knn),
                                 batch = batch,
                                 mod = model)
leduc <- addAssay(leduc, y = prot_knn, name = "prot_batch")
leduc <- addAssayLinkOneToOne(leduc, from = "prot_knn", to = "prot_batch")

prot_knn <- getWithColData(leduc,"prot_batch")
prot_knn_matrix <- assay(prot_knn)

rates = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)

full_replication <- na_replication(prot_knn_matrix, rates, 8)
full_replication_random <- na_replication_random(prot_knn_matrix, rates, 8)
no_na_rep <- map(full_replication, 1)



na_pca <- lapply(full_replication, function(x){
  lapply(x, pca_na_run)
  })

na_pca_random <- lapply(full_replication_random, function(x){
  lapply(x, pca_na_run)
})

svd_pca <- lapply(no_na_rep, pcaMethods_wraper, method = "svd")

saveRDS(na_pca_random, file = "data_output/na_pca_random_rep_v3.rds")
saveRDS(na_pca, file = "data_output/na_pca_rep_v3.rds")
saveRDS(svd_pca, file = "data_output/svd_rep_v3.rds")

################################################
              # kmn computation #
################################################

na_kmn <- lapply(na_pca, function(x){
  lapply(x, function(y){
    lapply(y, kmeans_wraper)
  })
})

na_kmn_random <- lapply(na_pca_random, function(x){
  lapply(x, function(y){
    lapply(y, kmeans_wraper)
  })
})

saveRDS(na_kmn_random, file = "data_output/na_kmn_random_v3.rds")

saveRDS(na_kmn, file = "data_output/na_kmn_v3.rds")

svd_kmn <- lapply(svd_pca, kmeans_wraper)

saveRDS(svd_kmn, file = "data_output/svd_kmn_v3.rds")


###################################################
              # PCA for complete subset #
###################################################
prots <- getWithColData(leduc,
                        "proteins_norm2")

prots_batched <- removeBatchEffect(assay(prots), batch = colData(prots)$Set)
prots_batched <- as.data.frame(prots_batched)
prots_comp <- prots_batched %>% 
  select(where(~sum(is.na(.x))/length(.x) < 0.6))
prots_comp <- prots_comp[complete.cases(prots_comp),]
prots_comp <- as.matrix(prots_comp)
comp_rep <- na_replication_random(prots_comp, rates, 8)
no_na_comp_rep <- map(comp_rep, 1)

comp_pca <- lapply(comp_rep, function(x){
  lapply(x, pca_na_run)
})


svd_comp <- lapply(no_na_comp_rep, pcaMethods_wraper, method = "svd")

saveRDS(comp_pca, file = "data_output/comp_pca_v3.rds")
saveRDS(svd_comp, file = "data_output/svd_comp_v3.rds")

comp_pca_kmn <- lapply(comp_pca, function(x){
  lapply(x, function(y){
    lapply(y, kmeans_wraper)
  })
})

saveRDS(comp_pca_kmn, file = "data_output/comp_pca_kmn_v3.rds")

svd_comp_kmn <- lapply(svd_comp, kmeans_wraper)
saveRDS(svd_comp_kmn, file = "data_output/svd_comp_kmn_v3.rds")


end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken