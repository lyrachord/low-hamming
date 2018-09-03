# Low Hamming key recovery attack

C++ and GMP bignum library

In C/, "make" and then ./o.o 
if list does not exist, generate it with "make list && ./list.o" and save output to precomp/<block_size>-<max_weight_per_block>


## Benchmarks for log2(n)= 2048

### Brute-force security = 115 bits
log2(p) = 1024
Hamming(p) = 16
Poisson = 1/64
Key ocurrence: High (1/10 ?)

Block_size = 100, scope = 3, max_weight = 3.
Runtime: 2^16, 98 seconds (26 minutes GPU)

### Brute-force security = 138 bits
log2(p) = 1024
Hamming(p) = 20
Poisson = 1/52
Key ocurrence: 1/10 (?)

Block_size = 70, scope = 4, max_weight = 3.
Runtime: 2^29, 45 seconds (68 days GPU)


## Benchmarks for log2(n)= 4096

### Brute-force security = 88 bits
log2(p) = 2048
Hamming(p) = 10
Poisson = 1/205
Key ocurrence: Medium (25% ?)

Block_size = 205, scope = 6, max_weight = 2.
Runtime: 2^24, 36 seconds (1 day GPU)

### Brute-force security = 124 bits
log2(p) = 2048
Hamming(p) = 15
Poisson: 1/136
Key occurrence: Rare (1/150?)

Block_size = 200, scope = 6, max_weight = 2
Runtime: 2^26, 40 seconds (7 days GPU)


