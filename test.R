read.fa <- function(fn){
    lines <- readLines(fn)
    id.i <- grep("^>", lines)
    beg.i <- id.i + 1
    end.i <- c(id.i-1, length(lines))[-1]
    seq <- mapply(function(b, e){
        paste(lines[b:e], collapse="")},
        beg.i, end.i)
    names(seq) <- lines[id.i]
    seq
}

seq <- read.fa("test.fa") ## so slow compared to command line grep..

nchar(seq)
## >scaffold_1 >scaffold_2 >scaffold_3 >scaffold_4 
##    48952293    44112721    42283995    42113633 

## dyn.load( "src/nuc_complexity.so")

## use the wrapper function:
source("nuc_complexity.R")

entropy.ol <- kmer.complexity(seq[1], c(5, 7), 1024, TRUE)
entropy.no <- kmer.complexity(seq[1], c(5, 7), 1024, FALSE)

## entropy.ol <- .Call("dna_entropy", seq[1], c(5L, 7L), 1024L, TRUE)
## entropy.no <- .Call("dna_entropy", seq[1], c(5L, 7L), 1024L, FALSE)

par(mfrow=c(2,1))
plot(entropy.ol[[1]], type='l')
lines(1:length(entropy.no[[1]]) * 5, entropy.no[[1]], col='red')
##
plot(entropy.ol[[2]], type='l')
lines( 1:length(entropy.no[[2]]) * 7, entropy.no[[2]], col='red')

## 
entropy.ol <- kmer.complexity(seq[1], c(5, 7), 1024, olap=FALSE, comp.beg=FALSE)
entropy.ol.bc <- kmer.complexity(seq[1], c(5, 7), 1024, olap=FALSE, comp.beg=TRUE)
## 
par(mfrow=c(2,1))
x <- 1:100000
plot(x*entropy.ol$kws[1,'k'], entropy.ol[[1]][x], type='l')
lines(x*entropy.ol$kws[1,'k'], entropy.ol.bc[[1]][x], type='l', col='red')
abline(v=entropy.ol$kws[1,'ws'] * entropy.ol$kws[1,'k'], lty=2)
##
plot(x*entropy.ol$kws[2,'k'], entropy.ol[[2]][x], type='l')
lines(x*entropy.ol$kws[2,'k'], entropy.ol.bc[[2]][x], type='l', col='red')
abline(v=entropy.ol$kws[2,'ws'] * entropy.ol$kws[2,'k'], lty=2)

par(mfrow=c(4,1))
for(i in 1:length(seq)){
    entropy.ol.bc <- kmer.complexity(seq[i], c(5, 7), 1024, olap=FALSE, comp.beg=TRUE)
    x <- with(entropy.ol.bc, 1:length(k5) * 5)
    plot(x, entropy.ol.bc$k5, type='l')
}

i <- 4
system.time(
    entropy.ol <- kmer.complexity(seq[i], c(5), 1024, olap=TRUE, comp.beg=TRUE)
)
##  user  system elapsed 
## 4.075   0.185   4.260 

i <- 4
system.time(
    entropy.no <- kmer.complexity(seq[i], c(5), 1024, olap=FALSE, comp.beg=TRUE)
)
##  user  system elapsed 
## 1.695   0.050   1.744 

pdf("entropy_example.pdf", width=12, height=5)
with(entropy.ol, plot( 1:length(k5), k5, type='l', lwd=0.1, xlab='position', ylab='Entropy'))
with(entropy.no, lines( 1:length(k5) * kws['k5','k'], k5, type='l', lwd=0.1, col='red'))
dev.off()


entropy.ol <- .Call("dna_entropy", substring(seq[1], 1, 100000), c(3L), 0L, TRUE)
entropy.no <- .Call("dna_entropy", substring(seq[1], 1, 100000), c(3L), 0L, FALSE)

x1 <- 1:nrow(entropy.ol)
x2 <- 3 * 1:length(entropy.no[[1]])

par(mfrow=c(1,1))
plot(x1, -entropy.ol[,1], type='l')
lines(x2, -entropy.no[[1]], type='l', col='red')

par(mfrow=c(2,1))
plot(x1, -entropy.ol[,1], type='l')
plot(x2, -entropy.no[[1]], type='l')


entropy.ol <- .Call("dna_entropy", substring(seq[1], 1, 100000), c(3L), 64L * 3, TRUE)
entropy.no <- .Call("dna_entropy", substring(seq[1], 1, 100000), c(3L), 64L, FALSE)

par(mfrow=c(1,1))
plot(x1, -entropy.ol[,1], type='l')
lines(x2, -entropy.no[[1]], type='l', col='red')


system.time(
    entropy <- .Call("dna_entropy", seq[4], c(1L, 2L, 3L, 5L), 0L, TRUE)
)
##   user  system elapsed 
##  15.985   0.957  16.940 


system.time(
    entropy <- .Call("dna_entropy", seq[4], c(1L, 2L, 3L, 5L), 0L, FALSE)
)
##   user  system elapsed 
## 11.253   0.338  11.589
##
## surprisingly, not that much faster; but maybe because of the very short
## kmer lengths used:


system.time(
    entropy <- .Call("dna_entropy", seq[4], c(3L, 5L), 0L, TRUE)
)
##   user  system elapsed 
##  8.113   0.448   8.561 


system.time(
    entropy <- .Call("dna_entropy", seq[4], c(3L, 5L), 0L, FALSE)
)
##   user  system elapsed 
##  3.843   0.083   3.926 

system.time(
    entropy <- .Call("dna_entropy", seq[4], c(7L), 0L, TRUE)
)
##  user  system elapsed 
## 4.226   0.227   4.453 

system.time(
    entropy <- .Call("dna_entropy", seq[1], c(5L, 7L), 1000L, FALSE)
)
##  user  system elapsed 
##  1.490   0.029   1.519 

x1 <- 1:1e6
x2 <- seq(1, 1e6, 5)

plot(x2, -entropy[[1]][1:length(x2)], type='l')

plot(-entropy[[1]], type='l')
plot(-entropy[[2]], type='l')
