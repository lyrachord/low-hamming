# Low Hamming key recovery attack

C++ and GMP bignum library

In C/, "make" and then ./o.o. If list does not exist, generate it with "make list && ./list.o" and save output to precomp/<block_size>-<max_weight_per_block>


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

# Idea: Use Gröebner basis.

Algebraic Cryptanalysis of Hidden Field Equation (HFE) Cryptosystems Using Grobner Bases (Joux)
Using Algebraic Geometry (Cox)

# References


## Factoring Unbalanced Moduli with Known Bits 

### ICISC 2009, Eric Brier, David Naccache, Mehdi Tibouchi

  LLL method.
  Unbalanced n = pq > q^3. Factorize when (i) 2log(q) contiguous bits of p are known, (ii) two chunks of known bits of p totaling 1.5 log(q) bits.


## Implicit Factoring: On Polynomial Time Factoring Given Only an Implicit Hint
### PKC 2009, Alexander May, Maike Ritzenhofen

  Good introduction.
  Not necessarily balanced N=pq, factorize with access to an oracle that outputs N'=p'q' such that p,p' share the t LSBs.
  We have t >= k/(k-1) alpha where k queries have been made, and the outputs are k moduli N_1, ... , N_k, each with log(q_i) = alpha.
  TODO: Lattice reduction of section 3.

## Further Results on Implicit Factoring in Polynomial Time
### 2009 Advances in Mathematics of Communications

  Extension of May et al PKC 2009. Same problem, improved results.

## Factoring RSA moduli with weak prime factors
### C2SI-Berger2015 

  A prime p is called weak if there exist parameters a, u_0, u_1, ... , u_k, M_1, ... , M_k such that ap = u_0 + u_1M_1 + ... + u_kM_k.
  This may yield our results: a = 1, u_0 = 1, u_i = 1, M_i = 2^i for some index i\in I.
  TODO: Check bounds on parameters.
  Does not cover the case "HD(p,q) small". 

## Factoring RSA Modulus Using Prime Reconstruction from Random Known Bits
### AFRICACRYPT 2010, Miatra, Sarkar, Sen Gupta


  N = pq balanced, when some random subset of the bits are known in the LSB halves. And give directions for MSB.

## On Factoring Arbitrary Integers with Known Bits
### 2007 Workshop "Kryptologie in Theorie und Praxis"

  TODO:
  "Rivest and Shamir [RS86] showed in 1985, that for an RSA-modulus N = pq an amount of 1/3 log N of the bits of p is sufficient to factor N. This result was improved by Coppersmith [Cop96] in 1996 to 3/10 log N bits, and in 1997 again by Coppersmith [Cop97] to 1/4 log N bits. In 1999, Boneh, Durfee and Howgrave-Graham [BDHG99] generalized the Coppersmith result to moduli of the form N=p^k q. They showed that k/(k+1)^2 log N bits are sufficient to find the factorization of N in polynomial time. One should notice that this result coincides with the one of Coppersmith for the RSA case, where  k = 1."

## Factoring RSA Moduli. Part I.
### Blog. 
  
  https://windowsontheory.org/2012/05/15/979/
  TODO: Download this before it dissapears. Part II too.


## Twenty Years of Attacks on the RSA Cryptosystem
### 1999 Dan Boneh

  Nice reading. 
  Theorem 10: (Coppersmith)
  Given the log(n)/4 LBSs or the log(n)/4 MSBs of p, one can factorize RSA moduli (note that this is RSA, not factoring, i.e. it uses ed=1 mod phi(n)).



## A Coding-Theoretic Approach to Recovering Noisy RSA Keys
### Asiacrypt 2012

  Cold boot attacks, when a key is recovered but allowing error in some bits. I.e. recover N=pq from N=(p+e)(q+e') where e,e' are of small Hamming weight.
  In the introduction they describe a very similar algorithm to ours (Correcting Errors in RSA Private Keys, see below)


## Correcting Errors in RSA Private Keys
### Crypto 2010

  TODO. The low Hamming hypothesis reduces to this (where the case that the recovered "noisy" key is p=0, and we want to "correct" h bits)...


## Using LLL-Reduction for Solving RSA and Factorization Problems: A Survey
### Alexander May, chapter of the book The LLL Algorithm: Survey and Applications

  Nice read. TODO. 






