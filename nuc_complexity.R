dyn.load( paste(dirname(sys.frame(1)$ofile), "src/nuc_complexity.so", sep="/") )

kmer.complexity <- function(seq, k, window.size, olap, comp.beg=FALSE){
    k <- as.integer(k)
    window.size <- as.integer(window.size)
    olap <- as.logical(olap)
    entropy <- .Call("kmer_complexity", seq[1], k, window.size, olap)
    names(entropy) <- c(paste("k", entropy[[length(entropy)]][,1], sep=""), "kws")
    colnames(entropy$kws) <- c("k", "ws")
    rownames(entropy$kws) <- names(entropy)[-length(entropy)]
    entropy$olap <- olap
    if(comp.beg){
        for(i in 1:nrow(entropy$kws)){
            ws <- entropy$kws[i,'ws']
            mod <- ws / (1:ws)
            entropy[[i]][ 1:ws ] <- entropy[[i]][ 1:ws ] * mod
        }
    }
    entropy
}

