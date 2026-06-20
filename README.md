# An Improved Upper Bound on Romanoff's Constant

## Summary

Romanoff's Constant is the asymptotic density of the set of natural numbers that can be expressed as the sum of a prime and a power of two.
See https://teorth.github.io/optimizationproblems/constants/45a.html for an overview.

This repository improves the best known upper bound from the one provided by [Griego](https://github.com/sebastian-griego/c45-romanoff-certificate):

$$C_{45} < 0.490249407811155$$

to

$$C_{45} < 0.490248930138967$$

which is an improvement of roughly $4.8 \cdot 10^{-7}$

## Method

We use optimization techniques to improve Griego's prime set. The new 37-prime set is:

$$\{3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 47, 61, 73, 89, 167, 223, 233, 263, 359, 383, 431, 439, 479, 503, 1103, 1399, 1433, 1913, 2089, 2351, 2687, 4513, 6361, 8191, 9719\}$$


## Verification

### Requirements

- C++ compiler (g++ 13.3.0)
- Python (Python 3.12.3)

### Instructions

```sh
./verify.sh
```

Computes the upper bound using exact (rational) arithmetic.

## Notes

Verification scripts were adopted from Griego.
Prepared with assistance from DeepSeek V4.
