# nuc_complexity: windowed kmer entropy estimates

This code contains a single function, `kmer.complexity` that calculates
entropy ($E$) based on the kmer content across specified windows:

$$
E = -\sum_{i=1..n}log(f_i)f_i
$$

Where $f_i$ is the frequency of kmer $i$ and $n$ is the number of k-mers
for a specified window.

``` R
kmer.complexity(seq, k, window.size, olap, comp.beg=FALSE)
```

## Arguments

The function takes 6 arguments:

1. seq  
   A character vector. Entropy will be determined for only the first element of this vector.
2. k  
   A numeric vector specifying the size of the k-mers to be used for the analysis. Any
   number of values of k can be specified, but they should be within 1 and 16.
3. window.size  
   The size of windows across which entropy should be calculated. If set to 0, the window
   size will be set to $4^k$ for each `k` specified. If a vector of window sizes is specified
   this will be recycled for each value of `k`.
4. olap  
   Should the complexity be calculated for all overlapping k-mers or only for distinct ones?
   That is, should the kmers be moved by a step of `1` or `k` during the determination.
5. comp.beg  
   Compensate the entropy calculated for the incomplete windows. If true, the frequency $f_i$
   will be determined for the number of k-mers counted until a complete window is obtained. 
   This means that the initial entropy will always be the maximal possible, and this will then
   drop until a complete window has been obtained.
   
The function uses $log_2$ in order that the unit of entropy is in bits. This is convenient as a
completely random nucleid acid sequence will have an entropy of $2k$. Since it is possible to
specify windows that are smaller than the k-mer, the max entropy $E_{max}$ is:

$$
E_{max} = \min
	\begin{cases}
	2k \\
	log_2{s}
	\end{cases}
$$

where $s$ is the window size used.

## Value

`kmer.complexity` returns a named list containing the following components:

1. `k<k>`: A vector of entropy for each value of `k` specified. These will be named "k<k>"
   with `<k>` representing the size of `k`. That is, if `k` is 5, the entry will
   be named `k5`.
2. `kws`: a matrix with columns `k` and `ws` giving the kmer and window sizes for which
   entropy has been estimated.
3. `olap`: A logical value indicating whether entropy was estimated
   for overlapping kmers.

## Example

Here I have run the function scaffold 4 of an assembly of the genome of
*L. piscatorius*. This scaffold appears to include telomere sequences at both
ends and contains about 42 million base pairs.

``` R
## seq contains genome assembly scaffolds
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
```

It takes about 4 and 1.7 seconds respectively for overlapping and non-overlapping
kmers with a k of 5. There is no discernible differences in the estimates for
overlapping and non-overlapping kmers at this resolution.

![Complexity along a scaffold from the *L. piscatorius* genome. Black and red lines
indicate estimates based on overlapping and non-overlapping kmers respectively.
Low entropy at the beginning and end of the scaffold reflect 
telomere sequences](entropy_example.png)


## Warning

### No checks for ambiguity symbols

The function uses a 2 bit representation of kmers by modifying a single 32 bit
unsigned integer using the following bitwise operations:

``` c
kmer = (kmer << 2) | ( (seq[i] >> 1) & 3 );
```

This means that any stretches containing ambiguity symbols (especially
`N`) will give incorrect entropy estimates without warning. The most
likely problem is that all windows that contain `N`s will have a
reduced entropy. I do not consider this to be a major issue since the
estimates themselves do not consider sequence content and this is
something that should be extracted as a secondary step.

### Window size

The window size is specified by the number of k-mers counted; this means that
the physical window size (the size of the region of the chromosome) is different
for overlapping and non-overlapping k-mers.


