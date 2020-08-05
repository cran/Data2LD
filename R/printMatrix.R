printMatrix <- function(instr, matval, rowind=1:m, colind=1:n, roundval=4) {
matval = as.matrix(matval)
matdim = dim(matval)
m = matdim[1]
n = matdim[2]
print(paste(instr,":"))
print(round(matval[rowind,colind],roundval))
}