#Factoring Unbalanced Moduli with Known Bits 
## ICISC 2009, Eric Brier, David Naccache, Mehdi Tibouchi

  Poorly written. LLL method.
  Unbalanced n = pq > q^3. Factorize when (i) 2log(q) contiguous bits of p are known, (ii) two chunks of known bits of p totaling 1.5 log(q) bits.


#Implicit Factoring: On Polynomial Time Factoring Given Only an Implicit Hint
## PKC 2009, Alexander May, Maike Ritzenhofen

  Good introduction.
  Not necessarily balanced N=pq, factorize with access to an oracle that outputs N'=p'q' such that p,p' share the t LSBs.
  We have t >= k/(k-1) alpha where k queries have been made, and the outputs are k moduli N_1, ... , N_k, each with log(q_i) = alpha.
  TODO: Lattice reduction of section 3.

# Factoring RSA Modulus Using Prime Reconstruction from Random Known Bits
## AFRICACRYPT 2010, Miatra, Sarkar, Sen Gupta


  N = pq balanced, when some random subset of the bits are known in the LSB halves. And give directions for MSB.

#On Factoring Arbitrary Integers with Known Bits
##2007 Workshop "Kryptologie in Theorie und Praxis"

  TODO:
  "Rivest and Shamir [RS86] showed in 1985, that for an RSA-modulus N = pq an amount of 1/3 log N of the bits of p is sufficient to factor N. This result was improved by Coppersmith [Cop96] in 1996 to 3/10 log N bits, and in 1997 again by Coppersmith [Cop97] to 1/4 log N bits. In 1999, Boneh, Durfee and Howgrave-Graham [BDHG99] generalized the Coppersmith result to moduli of the form N=p^k q. They showed that k/(k+1)^2 log N bits are sufficient to find the factorization of N in polynomial time. One should notice that this result coincides with the one of Coppersmith for the RSA case, where  k = 1."

# Factoring RSA Moduli. Part I.
## Blog. 
  
  https://windowsontheory.org/2012/05/15/979/
  TODO: Download this before it dissapears. Part II too.


# Twenty Years of Attacks on the RSA Cryptosystem
## 1999 Dan Boneh

  Nice reading. Does not say much about this.
  TODO: Coppersmith


# A Coding-Theoretic Approach to Recovering Noisy RSA Keys
## Asiacrypt 2012

  Cold boot attacks, when a key is recovered but allowing error in some bits. I.e. recover N=pq from N=(p+e)(q+e') where e,e' are of small Hamming weight.
  In the introduction they describe a very similar algorithm to ours (Correcting Errors in RSA Private Keys, see below)


# Correcting Errors in RSA Private Keys
## Crypto 2010

  TODO. The low Hamming hypothesis reduces to this (where the case that the recovered "noisy" key is p=0, and we want to "correct" h bits)...

  


