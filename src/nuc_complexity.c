#include <R.h>
#include <Rinternals.h>

/*
  Calculate sequence complexity over specified window sizes
  The functions given here will take a sequence (expected to
  be long) and calculate the entropy of the sequence based
  on the kmer frequency spectrum within sliding windows.
  
  The kmer counts will be kept in a single array with 
  the kmer index simply being the two-bit representation
  of A, C, G, T obtained from concatenation of bits 2 and 3
  of the individual bytes. If c is the nucleotide and b is
  the two bit representation:
  
  b <- (c >> 1) & 3
  This maps as:
  A: 0, C: 1, G: 3, T: 2.
  
  Hence to update the leading kmer represented as an unsigned
  integer (or possibly long) we can simply do
  
  kmer_right <- (kmer_right << 2) | ((c >> 1) & 3)
  
  and to update the counts:
  
  counts[kmer_right]++
  
  Note that we will also need to keep track of a lagging kmer; that
  is the kmer at the beginning of the window for which we calculate
  the complexity.
  This will be updated exactly as above but taking the character c
  from the last character of the first kmer of the window. 
  
  Note that here we need to decrement
  the counts before we update the kmer:
  
  counts[kmer_left]--
  kmer_left <- (kmer_right << 2) | ((c >> 1) & 3)
  
  
  The complexity of each window can be calculated as the entropy (E) of
  a single letter from it's alphabet, where a letter here refers to
  a kmer:
  
  E = -sum( log2(f_i) * f_i )
  
  where f_i is the frequency of word i and the sum is obtained from
  the frequencies of all possible words. 

  This means that to update E as we move the window we need to subtract
  the previous frequency term and add the new one for both the leading
  and lagging kmers (kmer_right and kmer_left)
  
  E <- E - (log2(f_io) * f_io) + (log2(f_no) * f_no)
  
  Where f_io and f_in are the old and new frequencies for the word
  whose count was updated. The old frequencies need to be calculated from
  the counts table prior to the updating of the counts.

  Note that a single word can be used to represent different kmers by simply
  masking appropriately.

*/

/* A struct that holds the count arrays for a set of k-mer sizes
   In this incarnation it only supports k up to 16. But we probably
   do not want to use such long k-mers for this purpose anyway.
   k_n       number of distinct kmers used
   k         the kmer lengths (k)
   count_l   the lengths of the count arrays
   counts    the counts. Note that if k=1, this might need to be 64 bit
   masks     masks the upper bits of a word.
   entropy   the negative entropy for the current state.
 */ 
#define MAX_K 16
struct k_mer_counts {
  unsigned int k_n;
  unsigned int *k;
  unsigned int *counts_l;
  unsigned int *window_s;
  unsigned int **counts;
  unsigned int *masks;
  double *entropy;
};

  
// Returns a zeroed counts for all k <= MAX_K;
// if all k are too large, k_n will be 0.
// This should be checked by the caller.
struct k_mer_counts init_kmer_counts(unsigned int k_n, unsigned int *k){
  struct k_mer_counts kc;
  memset(&kc, 0, sizeof(kc));
  kc.k = calloc( k_n, sizeof(unsigned int));
  kc.counts_l = calloc( k_n, sizeof(unsigned int));
  kc.window_s = calloc( k_n, sizeof(unsigned int));
  kc.counts = calloc( k_n, sizeof(unsigned int*));
  kc.masks = calloc( k_n, sizeof(unsigned int));
  kc.entropy = calloc(k_n, sizeof(double));
  for(unsigned int i=0; i < k_n; ++i){
    if(k[i] <= MAX_K){
      kc.k[ kc.k_n ] = k[i];
      kc.counts_l[ kc.k_n ] = (1 << (2 * k[i]));
      kc.window_s[ kc.k_n ] = kc.counts_l[ kc.k_n ];
      kc.counts[ kc.k_n ] = calloc( kc.counts_l[ kc.k_n ], sizeof(unsigned int) );
      kc.masks[ kc.k_n ] = kc.counts_l[ kc.k_n ] - 1;
      kc.k_n++;
    }
  }
  return(kc);
}

void clear_init_kmer_counts(struct k_mer_counts *kc){
  for(size_t i=0; i < kc->k_n; ++i)
    free( kc->counts[i] );
  free( kc->counts );
  free( kc->k );
  free( kc->counts_l );
  free( kc->window_s );
  free( kc->masks );
  free( kc->entropy );
  kc->k_n = 0;
}

// add a word to the indexed counts
// This will also calculate the updated entropy
// for the i-th k-mer size;
void add_word( struct k_mer_counts *kc, unsigned int word, size_t i ){
  double of = (double)kc->counts[i][ word & kc->masks[i] ] / (double)kc->window_s[i];
  double nf = (double)(kc->counts[i][ word & kc->masks[i] ] + 1) / (double)kc->window_s[i];
  kc->counts[i][ word & kc->masks[i] ]++;
  kc->entropy[i] -= ( log2(nf) * nf - (of > 0 ? (log2(of) * of) : 0) );
}

void rm_word( struct k_mer_counts *kc, unsigned int word, size_t i ){
  if(!kc->counts[i][ word & kc->masks[i] ]){
    Rprintf("Error at %d with word: %x\n", i, word & kc->masks[i]);
    return;
  }
  double of = (double)kc->counts[i][ word & kc->masks[i] ] / (double)kc->window_s[i];
  double nf = (double)(kc->counts[i][ word & kc->masks[i] ] - 1) / (double)kc->window_s[i];
  kc->counts[i][ word & kc->masks[i] ]--;
  kc->entropy[i] -= ( (nf > 0 ? log2(nf) * nf : 0) - (log2(of) * of) );
}

// entropy should be a preallocated block of doubles; where each column represents a kmer
// size. There should be one row for each residue of seq_l.
// This can be allocated with allocMatrix(REALSXP, seq_l, kc->k_n )
// That may be a massive matrix, but.. 
void scan_sequence(const char *seq, size_t seq_l, struct k_mer_counts *kc, double **entropy){
  unsigned int kmer_right;
  unsigned int *kmer_left = calloc( kc->k_n, sizeof(unsigned int) );
  for(size_t i=0; i < seq_l; ++i){
    kmer_right = (kmer_right << 2) | ( (seq[i] >> 1) & 3 );
    // then do different things depending on whether i >= the various ks in kc
    for(int j=0; j < kc->k_n; ++j){
      if( i+1 < kc->k[j] )
	continue;
      add_word( kc, kmer_right, j);
      if( i >= kc->window_s[j] )
	kmer_left[j] = (kmer_left[j] << 2) | ((seq[i - kc->window_s[j]] >> 1) & 3);
      if( i+1 >= kc->window_s[j] + kc->k[j] )
	rm_word( kc, kmer_left[j], j);
      entropy[j][i] = kc->entropy[j];
      //      entropy[ i + j * seq_l ] = kc->entropy[j];
    }
  }
  free(kmer_left);
}

// no: no overlap
// this means that the number of entropy values calculated is seq_l / k and that
// we need to have an array of arrays of variable size.
// it would be better to have this as an option to a single
// function. But first check what needs to be different.
void scan_sequence_no(const char *seq, size_t seq_l, struct k_mer_counts *kc, double **entropy){
  unsigned int kmer_right;
  unsigned int *kmer_left = calloc( kc->k_n, sizeof(unsigned int) );
  for(size_t i=0; i < seq_l; ++i){
    kmer_right = (kmer_right << 2) | ( (seq[i] >> 1) & 3 );
    // then do different things depending on whether i >= the various ks in kc
    for(int j=0; j < kc->k_n; ++j){
      if( i+1 < kc->k[j] )
	continue;
      if( i >= kc->window_s[j] * kc->k[j] )
	kmer_left[j] = (kmer_left[j] << 2) | ((seq[i - kc->window_s[j] * kc->k[j]] >> 1) & 3);
      if( (1 + i) % kc->k[j] != 0)
	continue;
      add_word( kc, kmer_right, j);
      if( i+1 >= (1 + kc->window_s[j]) * kc->k[j] )
	rm_word( kc, kmer_left[j], j);
      entropy[j][ i / kc->k[j] ] = kc->entropy[j];
    }
  }
  free(kmer_left);
}


// Although seq_r could contain multiple sequences, this function will only
// consider the first one
// seq_r: one or more sequences, the first of which will be analysed
// k_r;   one or more k_mer sizes
// ws_r:  window sizes; if the first value is 0, then default sizes will be used
// olap_r: count overlapping kmers; if false count only distinct kmers
SEXP kmer_complexity(SEXP seq_r, SEXP k_r, SEXP ws_r, SEXP olap_r){
  if(TYPEOF( seq_r ) != STRSXP || length(seq_r) < 1 )
    error("seq_r must be a character vector of length at least one");
  if(TYPEOF( k_r ) != INTSXP || length(k_r) < 1 )
    error("k_r must be an integer vector of length at least one");
  if(TYPEOF( ws_r ) != INTSXP || length(ws_r) < 1 )
    error("ws_r must be an integer vector of length at least one");
  if(TYPEOF(olap_r) != LGLSXP || length(olap_r) < 1)
    error("olap_r should be a logical vector of length 1");
  SEXP seq = STRING_ELT(seq_r, 0);
  int seq_l = length(seq);
  int k_n = length(k_r);
  int *k = INTEGER( k_r );
  int ws_n = length(ws_r);
  int *ws = INTEGER( ws_r );
  int olap = asLogical(olap_r);
  // go on without error checking.. 
  struct k_mer_counts kc = init_kmer_counts((unsigned int)k_n, (unsigned int*)k);
  if(ws[0]){
    for(size_t i=0; i < kc.k_n; ++i)
      kc.window_s[i] = ws[ i % ws_n ];
  }
  SEXP entropy_r = PROTECT(allocVector(VECSXP, kc.k_n+1));
  // the last of these should contain a matrix giving the respective
  // kmer and window sizes:
  SET_VECTOR_ELT(entropy_r, kc.k_n, allocMatrix(INTSXP, kc.k_n, 2));
  int *k_w_size = INTEGER(VECTOR_ELT(entropy_r, kc.k_n));
  for(size_t i=0; i < kc.k_n; ++i){
    k_w_size[i] = kc.k[i];
    k_w_size[i + kc.k_n] = kc.window_s[i];
  }
  double **entropy = malloc(sizeof(double*) * kc.k_n);
  for(int i=0; i < kc.k_n; ++i){
    size_t entropy_l = olap ? seq_l : seq_l / kc.k[i];
    SET_VECTOR_ELT( entropy_r, i, allocVector(REALSXP, entropy_l));
    entropy[i] = REAL(VECTOR_ELT(entropy_r, i));
    memset(entropy[i], 0, sizeof(int) * entropy_l); // unnecessary.
  }
  if(olap)
    scan_sequence(CHAR(seq), seq_l, &kc, entropy);
  else
    scan_sequence_no(CHAR(seq), seq_l, &kc, entropy);
  free(entropy);
  clear_init_kmer_counts(&kc);
  UNPROTECT(1);
  return(entropy_r);
}

static const R_CallMethodDef callMethods[] = {
  {"kmer_complexity", (DL_FUNC)&kmer_complexity, 4},
  {NULL, NULL, 0}
};

void R_init_nuc_complexity(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
