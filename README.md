# LDPC
Nonbinary LDPC codes encoding and decoding implementation

What is included:

1. Regular Nonbinary LDPC code builder. 
It gets code length *N*, redundancy *r*, Galois field exponent *m*, column degree *Dc* and row degree *Dr* as input. 
As a result, a check matrix of *(n, &le; (r-k))* regular *(Dr, Dc)* LDPC code over GF(2^m) it built as the block matrix consisting of random permutation-and-multiplication matrices with fixed dimensions and uniformly distributed nonzero elements.
This matrix is further stored in the specification file in following format:
  1. First row contains space-separated length, redundancy and field exponent
  2. Second row contains space-separated letter 'R' meaning 'Regular', row and column degrees
  3. Next r rows contain space-separated indices of nonzero elements along with their values in the corresponding rows of code check matrix
  
2. LDPC encoder. It extracts check matrix from specification file and gets a generator matrix from it. This matrix is further used for code encoding.

3. Log-domain sum-product decoding algorithm (see H. Wymeersch, H. Steendam, and M. Moeneclaey "Log-domain decoding
of LDPC codes over GF(q)").
