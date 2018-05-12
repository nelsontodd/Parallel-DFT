/************************************************************************
Comparing the performance increases we get from a dumb application of 
MPI to double radix FFT's, for both real and complex data. All implementations of FFT are taken from gsl (https://www.gnu.org/software/gsl/doc/html/fft.html)
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <omp.h>
#include <mpi.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

//#define REAL(z,i) (*(z + 2*(i)))
//#define IMAG(z,i) (*(z + 2*(i) + 1))
#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])
#define N 2048


static void bitreverse (double *data, const size_t stride, const size_t n,size_t logn)
{
  /* This is the Goldrader bit-reversal algorithm */

  size_t i;
  size_t j = 0;

  logn = 0 ; /* not needed for this algorithm */

  for (i = 0; i < n - 1; i++)
    {
      size_t k = n / 2 ;

      if (i < j)
        {
          const double tmp_real = REAL(data,i);
          const double tmp_imag = IMAG(data,i);
          REAL(data,i) = REAL(data,j);
          IMAG(data,i) = IMAG(data,j);
          REAL(data,j) = tmp_real;
          IMAG(data,j) = tmp_imag;
        }

      while (k <= j) 
        {
          j = j - k ;
          k = k / 2 ;
        }

      j += k ;
    }
}


static void create_complex_polynomial_signal(double *data, int length) {
  
  int i;
  //Random Coeffecients
  double a =(double) random()/RAND_MAX;
  double b =(double) random()/RAND_MAX;
  double c =(double) random()/RAND_MAX;
  for (i = 0; i < length; i++)
    {
       REAL(data,i) = a*(double)i*i + b*(double)i + c;
       IMAG(data,i) = a*(double)i*i + b*(double)i + c;
    }
}
static void complex_box_signal(double *data, int length) {
  
  int i;
  for (i = 0; i < length; i++)
    {
       REAL(data,i) = 0.0; IMAG(data,i) = 0.0;
    }

  REAL(data,0) = 1.0;

  for (i = 1; i <= length-5; i++)
    {
       REAL(data,i) = REAL(data,length-i) = 1.0;
    }
}

void main(int argc, char *argv[]) { 

  int i;
  FILE *orig, *complex_fft_output, *parallel_complex_fft_output;
  orig               = fopen( "original_data.txt", "w" ); // Open file for writing
  complex_fft_output = fopen( "complex_fft.txt", "w" ); // Open file for writing
  parallel_complex_fft_output = fopen( "parallel_complex_fft.txt", "w" ); // Open file for writing
  //Init MPI
  int world_size;
  int world_rank;
  int irec = 5;
  MPI_Init(&argc,&argv);
  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  MPI_Comm_rank(MPI_COMM_WORLD,&world_rank); 
  int elements_per_proc = (int) (2.0*(double)N)/(double)world_size;
  int sharded_arr_size          = (int) ((double)elements_per_proc/2.0);
  {
    fprintf(stderr, "");
  }


  double *complex_polynomial;
  double *test_element;
  double *original_data;
  double *shard = (double *)malloc(elements_per_proc * sizeof(double));
  if (world_rank == 0) {
    complex_polynomial = (double *)malloc((2*N) * sizeof(double));
    original_data      = (double *)malloc((2*N) * sizeof(double));
    test_element       = (double *)malloc((2*N) * sizeof(double));
    create_complex_polynomial_signal(complex_polynomial, N);
    memcpy(original_data, complex_polynomial, 2*N * sizeof(double));
    bitreverse(complex_polynomial, 1, N, 0);
  }
  bitreverse(shard, 1, sharded_arr_size, 0);
  MPI_Scatter(complex_polynomial, elements_per_proc, MPI_DOUBLE, shard,
               elements_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  printf ("\n\n");
  gsl_fft_complex_radix2_forward (shard, 1, sharded_arr_size);

  //Collect back into root to print it out
  MPI_Gather(shard, elements_per_proc, MPI_DOUBLE, test_element, elements_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  if (world_rank == 0) {
    
    //gsl_fft_complex_radix2_forward (complex_polynomial, 1, N);
    //fprintf(orig,"Output of original data: \n");
    for (i = 0; i < N; i++)
    {
      fprintf(orig,"%d %e %e\n", i,
      REAL(original_data,i)/sqrt(N),
      IMAG(original_data,i)/sqrt(N));
    }
    //fprintf(complex_fft_output,"Collected shards: \n");
    gsl_fft_complex_radix2_inverse (test_element, 1, N);
    for (i = 0; i < N; i++)
    {
      fprintf(parallel_complex_fft_output, "%d %e %e\n", i,
      REAL(test_element,i)/sqrt(N),
      IMAG(test_element,i)/sqrt(N));
    }

    gsl_fft_complex_radix2_forward (complex_polynomial, 1, N);
    gsl_fft_complex_radix2_inverse (complex_polynomial, 1, N);
    for (i = 0; i < N; i++)
    {
      fprintf(complex_fft_output, "%d %e %e\n", i,
      REAL(complex_polynomial,i)/sqrt(N),
      IMAG(complex_polynomial,i)/sqrt(N));
    }
  }
  //for (i = 0; i < sharded_arr_size; i++)
  //{
      //printf ("%d %e %e\n", i,
              //REAL(shard,i), IMAG(shard,i));
  //}

  MPI_Finalize();
}
