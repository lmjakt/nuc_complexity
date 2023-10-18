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

plot.entropy <- function(ent, xlab="position", ylab="entropy", main="", cols=1:nrow(ent$kws), comp.w=TRUE){
    l <- length(ent)
    ent.range <- range( sapply(ent[3:l-2], range) )
    x <- 1:length(ent[[1]])
    if(!ent$olap)
        x <- x * ent$kws[1,'k']
    plot(1,1, xlim=range(x), ylim=ent.range, xlab=xlab, ylab=ylab, main=main, type='n')
    for(i in 1:nrow(ent$kws)){
        x <- 1:length(ent[[i]])
        if(!ent$olap)
            x <- x * ent$kws[i,'k']
        beg <- ifelse(comp.w, ent$kws[i,'ws'], 1)
        r <- beg:length(x)
        lines(x[r], ent[[i]][r], col=cols[i])
    }
}

## The following is in order to obtain an encoding that can be used
## with fft() to obtain the power spectrum for each nucleotide
fft.encode <- function(seq, center=TRUE){
    toN <- function(char){
        as.numeric( toupper(utf8ToInt(seq)) == utf8ToInt(char) )
    }
    m <- sapply(c(A="A", C="C", G="G", T="T"), toN)
    if(center)
        m <- t( t(m) - colMeans(m) )
    m
}
