# Readme file for Symbolic Sequence Processing repository

### Daniel J. Greenhoe
### https://github.com/dgreenhoe/symbolic-sequence-processing.git


## Project Desription

This project is C++ code to support symbolic sequence processing, 
as described in the paper Greenhoe(2016).
Below is an abstract for that paper:


## Abstract for Greenhoe(2016)

A real-valued random variable X is a measurable function that maps from a 
probability space (Omega,E,P) to (R,<=,d,+,.,B) 
where B is the usual Borel sigma-algebra on the "real line" (X,<=,d,+,.), 
and where R is the set of real numbers, <= is the standard linear order relation on R, 
d(x,y)=|x-y| is the usual metric on R, and (R, +, .) is the standard field on R. 
Greenhoe(2015) has demonstrated that this definition of random variable is 
often a poor choice for computing statistics when the stochastic process that X maps 
from has structure that is dissimilar to that of the real line. 
Greenhoe(2015) has further proposed an alternative statistical system, 
that rather than mapping a stochastic process to the real line, 
instead maps to a weighted graph that has order and metric geometry structures 
similar to that of the underlying stochastic process. 
In particular, ideally the structure X maps from and the structure X maps to are, 
with respect to each other, both isomorphic and isometric. 
Such a mapping is useful for analysis of a single random variable---for example 
the expectation EX of X can be defined simply as the center of its weighted graph. 
However, the mapping has limitations with regards to a sequence of random variables in 
performing sequence analysis (using for example Fourier analysis or wavelet analysis) 
or in performing sequence processing (using for example FIR filtering or IIR filtering). 
Rather than mapping to a weighted graph, this present paper proposes instead mapping to an 
ordered distance linear space Y=(R^n,<=,d,+,.,R,+,x), where (R,+,x) is a field, + is the 
vector addition operator on R^n x R^n, and . is the scalar-vector multiplication operator on R x R^n. 
The linear space component of Y provides a much more convenient (as compared to the weighted graph) 
framework for sequence analysis and processing. 
The ordered set and distance space components of Y allow one to preserve the order structure and distance geometry 
inherent in the underlying stochastic process, which in turn likely provides a less distorted (as compared to the real line) 
framework for analysis and sequence processing. 

## References

  * Greenhoe(2015): Greenhoe, Daniel J.,
    "Order and metric geometry compatible stochastic processing", 
    2015 February 19, https://doi.org/10.7287/peerj.preprints.844v1

  * Greenhoe(2016): Greenhoe, Daniel J.,
    "Order and Metric Compatible Symbolic Sequence Processing",
    2016 May 18, https://doi.org/10.7287/peerj.preprints.2052v1
