# Low Hamming key recovery attack

C++ and GMP bignum library


## Benchmarks for log2(n)= 2048

### Brute-force security = 138 bits
log2(p) = 1024
Hamming(p) = 20

Runtime: 2^37, 3.5 seconds (1 GPU year)
Block_size = 70, scope = 6, max_weight = 3.


### Brute-force security = 115 bits
log2(p) = 1024
Hamming(p) = 16

Runtime: 2^25, 4 seconds (9 GPU hours)
Block_size = 85, scope = 4, max_weight = 2.

## Benchmarks for log2(n)= 4096

### Brute-force security = 88 bits
log2(p) = 2048
Hamming(p) = 10

Runtime: 2^24, 36 seconds (1 GPU day)
Block_size = 205, scope = 6, max_weight = 2.

### Brute-force security = 124 bits
log2(p) = 2048
Hamming(p) = 15

Runtime: 2^26, 46 seconds (8 GPU days)
Block_size = 200, scope = 6, max_weight = 2


