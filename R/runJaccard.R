#' @include utilities.R
globalVariables(names = 'i', package = 'SnapATAC', add = TRUE)
#' Calcualte Jaccard Index Matrix
#'
#' This function takes a snap object as input and calculates the jaccard index matrix
#' in which each entry Jij equals the intersection over the union between cell i and cell j.
#'
#' Calculating jaccard index becomes exponentially time-consuming and also memory
#' intense with the increase of cell number.
#'
#' The most memory and time consuming step of our procedure is the calculation of jaccard index matrix.
#' This step increases exponentially with the increase of cell number. Here, we calculate a partial jaccard
#' index matrix against a random subset of reference cells in lieu of the full spectrum of features.
#' Here we tested if we can perform the calculation of partial jaccard index matrix using a random
#' subset of max.var reference cells.
#'
#' @param obj A Snap obj
#' @param bin.downsample Percentage of bins to be downsampled to [1].
#' @param tmp.folder A non-empty character vector giving the directory name that saves the temp files
#' @param mat A character class that indicates what matrix slot is used to calculate jaccard index c("bmat", "pmat", "gmat")
#' @param max.var A numeric variable indicates the how many dimentions for jaccard index to be calcualted
#' @param seed.use A numeric variable indicates the random seed to use [10].
#'
#' @examples
#' data(demo.sp);
#' demo.sp = makeBinary(demo.sp, mat="bmat");
#' demo.sp = runJaccard(obj=demo.sp, mat="bmat", bin.downsample=1, tmp.folder=tempdir())
#'
#' @return Returns a Snap obj with the jaccard index matrix stored in obj@jmat
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom foreach foreach %dopar%
#' @importFrom stats lm approxfun density
#' @importFrom methods slot
#' @export
runJaccard <- function(obj, bin.downsample, tmp.folder, mat, max.var, seed.use) {
  UseMethod("runJaccard", obj);
}

#' @export
runJaccard.default <- function(
    obj,
	tmp.folder,
	mat = c("bmat", "pmat", "gmat"),
	max.var = 1000,
	ncell.chunk = 1000,
	do.par=FALSE,
	num.cores=1,
	seed.use=10
	){

	# check if mat is binary;
	mat = match.arg(mat);
	mat.use = slot(obj, mat);
	if((x=nrow(mat.use)) == 0L){
		stop("input matrix is empty");
	}

	if((x=max(mat.use)) > 1L){
		stop("input matrix is not a binary matrix, run 'makeBinary' first")
	}

	# check if empty rows exist in bmat;
	if(any(Matrix::rowSums(mat.use) == 0)){
		stop("input matrix contains empty rows, remove empty rows first")
	}

	max.var = min(max.var, nrow(mat.use));
	# randomly select a subset of cells as reference
	set.seed(seed.use);
	mat.ref = mat.use[sort(sample(seq(nrow(mat.use)), max.var)),];

	if(do.par==TRUE){
	    # input checking for parallel options
		if(num.cores > 1){
	        if (num.cores == 1) {
	          num.cores = 1
	        } else if (num.cores > detectCores()) {
	          num.cores <- detectCores() - 1
	          warning(paste0("num.cores set greater than number of available cores(", detectCores(), "). Setting num.cores to ", num.cores, "."))
	        }
	      } else if (num.cores != 1) {
	        num.cores <- 1
		}

		# step 2) slice the orginal obj into list
		id = seq(nrow(mat.use));
		id.ls = split(id, ceiling(seq(id)/ncell.chunk));

		if(length(id.ls) > 1){
			id.ls[[length(id.ls) - 1]] = c(id.ls[[length(id.ls) - 1]], id.ls[[length(id.ls)]]);
			# remove the last item of the list
			id.ls = id.ls[-length(id.ls)];
		}

		mat.list = splitBmat(mat.use, id.ls, num.cores, tmp.folder);

		cl <- makeCluster(num.cores);
		registerDoParallel(cl);

		jmat <- foreach(i=seq(mat.list), .verbose=FALSE,  .export=c("calJaccardSingle", "calJaccard"), .combine = "rbind") %dopar% {
			calJaccardSingle(mat.list[[i]], mat.ref)
		}
		stopCluster(cl);
		closeAllConnections();
		gc();

		rm(id.ls);
		lapply(mat.list, function(x){
			file.remove(x);
		})
		rm(mat.list);
		gc();

	}else{
		jmat = calJaccard(mat.use, mat.ref)
	}

	p1   <- Matrix::rowMeans(mat.use);
	p2   <- Matrix::rowMeans(mat.ref);
	rm(mat.ref);
	rm(mat.use);
	# remove large objects
	obj@jmat@jmat = jmat;
	obj@jmat@p1 = p1;
	obj@jmat@p2 = p2;
	obj@jmat@norm = FALSE;
	#obj@jmat@method = character();
	return(obj);
}


# Calculate jaccard index between two matrix
# @param X_i first matrix
# @param X_j second matrix
calJaccard <- function(X_i, X_j){
	A = Matrix::tcrossprod(X_i, X_j);
	bi = Matrix::rowSums(X_i);
	bj = Matrix::rowSums(X_j);
	jmat = as.matrix(A / (replicate(ncol(A), bi) + t(replicate(nrow(A), bj)) - A));
	return(jmat);
}
