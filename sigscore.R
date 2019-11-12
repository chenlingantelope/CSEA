library('Matrix')
library("VISION")

args = commandArgs(trailingOnly=TRUE)
signature = args[1]
exprs = args[2]
genenames = args[3]

# signature = 'sig/allsigs.txt'
# exprs = 'data/all_data.mtx'
# genenames = 'data/genenames.csv'

merged_exprs = readMM(exprs)
merged_exprs = t(merged_exprs)
exprs = as.matrix(merged_exprs)
cellsize = colSums(exprs)
ngenes = colSums(exprs>0)


genenames = read.csv(genenames,header=F,as.is=T)[,1]
rownames(exprs) = genenames

raw_exprs = exprs
n.umi = median(colSums(exprs))
exprs = apply(exprs, 2, function(x) (x * n.umi) / sum(x))

colnames(exprs) = as.character(c(1:dim(exprs)[2]))
colnames(raw_exprs) = as.character(c(1:dim(exprs)[2]))


vis <- Vision(data = exprs, nomodel = TRUE,
              unnormalizedData = raw_exprs,
              signatures = c(signature),
              cellsPerPartition = 1)


vis <- VISION:::calcSignatureScores(vis, sigData = vis@sigData,
	sig_norm_method = vis@sig_norm_method,
	sig_score_method = vis@sig_score_method)


write.csv(vis@sigScores,file=sprintf('%s.sigScores.csv', signature))
