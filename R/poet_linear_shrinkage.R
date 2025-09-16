#' poet_linear_shrinkage: Function to perform covariance matrix shrinkage using the POET method
#'
#' This function performs shrinkage of the covariance matrix using the POET (Principal Orthogonal complEment Thresholding) method, with an automatic selection of the number of latent factors and the optimal shrinkage parameter.
#'
#' @param S A symmetric matrix (p×p). It can be a covariance matrix or a
#'   correlation matrix.
#' @param n The sample size that produced `S`.
#' @param cutoff_method Character string specifying how to choose the number of
#'   factors `m`. One of `"D.ratio"` or `"ratio"` (default: `"D.ratio"`).
#' @param k_min,k_max Integers giving the search range for the number of
#'   factors `m`. Defaults are `k_min = 5` and
#'   `k_max = min(15, floor(nrow(S) / 2))`.
#' @param shrinkage_vec A vector of candidate shrinkage parameters. Default is `seq(0.025, 0.25, by = 0.025)`.
#' @param eigenmin The threshold for the minimum eigenvalues used in determining the optimal shrinkage parameter. Default is `1e-3`.
#'
#' @return A covariance matrix that has been shrunk using the POET method.
#'
#' @importFrom CppMatrix matrixInverse matrixMultiply matrixVectorMultiply matrixEigen matrixListProduct
#' @export
#'
poet_linear_shrinkage=function(S,n,cutoff_method="D.ratio",k_min=5,k_max=min(15,round(nrow(S)/2)),shrinkage_vec=seq(0.025,0.25,by=0.025),eigenmin=1e-3){
p=ncol(S)
S[is.na(S)]=0;
diag(S)=1
S[abs(S)<0.0001]=0
S=t(S)/2+S/2
eig=matrixEigen(S)
U=eig$vectors
d=eig$values
k_min <- max(k_min, 2L)
k_max <- min(k_max, p - 2L)

z <- rep(NA_real_, p)
if (cutoff_method == "D.ratio") {
j_lo <- max(k_min, 2L)
j_hi <- min(k_max, p - 1L)
if (j_lo <= j_hi) {
for (j in j_lo:j_hi) {
den <- (d[j] - d[j + 1L])
num <- (d[j - 1L] - d[j])
z[j - 1L] <- if (den != 0) num / den else NA_real_
}
}
} else if (cutoff_method == "ratio") {
j_lo <- max(k_min, 2L)
j_hi <- min(k_max, p)
if (j_lo <= j_hi) {
for (j in j_lo:j_hi) {
z[j - 1L] <- if (d[j] != 0) d[j - 1L] / d[j] else NA_real_
}
}
} else {
stop("Unknown cutoff_method.")
}

idx <- suppressWarnings(which.max(z))
if (length(idx) == 0L || is.na(idx)) {
pck <- max(k_min, 1L)
} else {
pck <- idx + 1L
}
pck <- min(max(pck, 1L), p - 1L)

Uk=U[,1:pck];dk=d[1:pck]
hatc=(sum(diag(S))-sum(dk))/(p-pck-pck*p/n)
dk <- pmax(dk - hatc*p/n, 0)
P=matrixMultiply(Uk,t(Uk)*dk)
E=S-P;e=diag(E);e[e<0]=max(hatc,0.01);diag(E)=e
eigenvec=shrinkage_vec
for(i in 1:length(shrinkage_vec)){
E1=E*(1-shrinkage_vec[i])+diag(diag(E))*shrinkage_vec[i]
eigenvec=min(matrixEigen(E1)$values)
if(eigenvec>eigenmin) break
}
hatS=P+E1
return(cov2cor(hatS))
}
