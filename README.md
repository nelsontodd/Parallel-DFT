# Abstract/Introduction
The Fourier Transform is used in signal processing to decompose signals into their component frequencies - called the 'fundamentals' by computing a   linear combination of sinusoidal functions. The Fast Fourier Transform is amethod of computing the Discrete Fourier Transform that reduces the complexity from $O(n^{2})$ to $O(n\log{n})$. For large N, this algorithm is significantly more efficient than DFT. Here, I implement a novel (but naive) implementation of DFT parallelized with OpenMPI and compare its scaling capabilities to Gnu Scientific Libraries FFT.

![alt text](https://github.com/nelsontodd/Parallel-DFT/blob/master/Paper/titlepage.png)
![alt text](https://github.com/nelsontodd/Parallel-DFT/blob/master/Paper/abstract.png)
![alt text](https://github.com/nelsontodd/Parallel-DFT/blob/master/Paper/code_and_discussion.png)
![alt text](https://github.com/nelsontodd/Parallel-DFT/blob/master/Paper/final_results.png)
![alt text](https://github.com/nelsontodd/Parallel-DFT/blob/master/Paper/scaling_infor.png)
