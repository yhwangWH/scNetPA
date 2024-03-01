normalizeData <- function(data.raw, scale.factor = 10000, do.log = TRUE) {
  # Scale counts within a sample
  library.size <- Matrix::colSums(data.raw)
  #scale.factor <- median(library.size)
  expr <- Matrix::t(Matrix::t(data.raw) / library.size) * scale.factor
  if (do.log) {
    data.norm <-log1p(expr)
  }
  return(data.norm)
}

exp_mat_path = "GSM5687481_k562_rep1_matrix.mtx.gz"
data.raw = Matrix::readMM(exp_mat_path) 
data.norm = normalizeData(data.raw, scale.factor = 10000, do.log = TRUE) 


data.rowMeans = Matrix::rowMeans(data.norm)

write.table(data.rowMeans, file = "normalizedData.txt", row.names = F, col.names = F)

gene_name_path = "GSM5687481_k562_rep1_features.tsv.gz"
gene_name_list = read.table(gene_name_path)
gene_name_list2 = gene_name_list[,2]
write.table(gene_name_list2, file = "normalizedDataGeneNameList.txt"
            , quote = F
            , row.names = F, col.names = F)