# Abstract/Introduction
The Fourier Transform is used in signal processing to decompose signals into their component frequencies - called the 'fundamentals' by computing a   linear combination of sinusoidal functions. The Fast Fourier Transform is amethod of computing the Discrete Fourier Transform that reduces the complexity from $O(n^{2})$ to $O(n\log{n})$. For large N, this algorithm is significantly more efficient than DFT. Here, I implement a novel (but naive) implementation of DFT parallelized with OpenMPI and compare its scaling capabilities to Gnu Scientific Libraries FFT.

# Scaling
The Tukey-Cooley algorithm has performance of $O(n\log{n})$ while DFT has performance of $O(N^{2})$. This difference is extremely significant, for large N the it can be the difference between seconds and days.
