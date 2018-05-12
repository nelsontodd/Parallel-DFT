/************************************************************************
Comparing the performance increases we get from a dumb application of 
MPI to double radix FFT's, for both real and complex data. All implementations of FFT are taken from gsl (https://www.gnu.org/software/gsl/doc/html/fft.html)
*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <mpi.h>
#include <omp.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])
//In order for GSL's FFT to work this must be 2^m for some (large) m
#define N 32768 //2^15

static void dft_kernel(double *input, double* sumreal, double* sumimag, int k, int min, int size) {
  double loc_real_sum = *sumreal;
  double loc_imag_sum = *sumimag;
  for (int t = 0; t < size; t++) {  // For each input element
      double angle = 2 * M_PI * (t+min) * k / N;
      loc_real_sum +=  REAL(input, t) * cos(angle) + IMAG(input, t)* sin(angle);
      loc_imag_sum += -REAL(input, t) * sin(angle) + IMAG(input, t)* cos(angle);
  }
  *sumreal = loc_real_sum;
  *sumimag = loc_imag_sum;
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

void main(int argc, char *argv[]) { 

  //Open files for writing
  FILE *orig_data, *serial_fft_out, *parallel_dft_out;
  orig_data    = fopen( "output/generated_signal.txt", "w" ); // Open file for writing
  serial_fft_out   = fopen( "output/serial_fft.txt", "w" ); // Open file for writing
  parallel_dft_out = fopen( "output/parallel_dft.txt", "w" ); // Open file for writing

  //Init MPI
  int world_size, world_rank, min;
  int irec = 5;
  double global_imag_sum, global_real_sum, local_real_sum, local_imag_sum;
  double times[6], start, end;
  MPI_Init(&argc,&argv);
  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  MPI_Comm_rank(MPI_COMM_WORLD,&world_rank); 
  int elements_per_proc = (int) (2.0*(double)N)/(double)world_size;
  int sharded_arr_size  = (int) ((double)elements_per_proc/2.0);
  {
    //Throwaway kernel
    fprintf(stderr, "");
  }


  double *complex_polynomial;
  double *parallel_dft;
  double *serial_dft;
  double *original_data;
  double *shard = (double *)malloc(elements_per_proc * sizeof(double));
  double *dft_shard = (double *)malloc(elements_per_proc * sizeof(double));
  if (world_rank == 0) {
    complex_polynomial = (double *)malloc((2*N) * sizeof(double));
    serial_dft         = (double *)malloc((2*N) * sizeof(double));
    original_data      = (double *)malloc((2*N) * sizeof(double));
    parallel_dft       = (double *)malloc((2*N) * sizeof(double));
    create_complex_polynomial_signal(complex_polynomial, N);
    memcpy(original_data, complex_polynomial, 2*N * sizeof(double));
  }
  MPI_Scatter(original_data, elements_per_proc, MPI_DOUBLE, shard,
               elements_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  start = omp_get_wtime();
  min =  sharded_arr_size*world_rank;
  for (int k = 0; k < N; k++) {  // For each output element
      local_real_sum = 0;
      local_imag_sum = 0;
      dft_kernel(shard, &local_real_sum, &local_imag_sum, k, min, sharded_arr_size);
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Reduce(&local_real_sum, &global_real_sum, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
      MPI_Reduce(&local_imag_sum, &global_imag_sum, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
      if (world_rank == 0) {
        REAL(parallel_dft, k) = global_real_sum;
        IMAG(parallel_dft, k) = global_imag_sum;
      }
  }
  end = omp_get_wtime();
  times[0]  = end-start;

  if (world_rank == 0) {

    int i;
    for (i = 0; i < 10; i++)
    {
      fprintf(orig_data, "%d %e %e\n", i,
      REAL(complex_polynomial,i)/sqrt(N),
      IMAG(complex_polynomial,i)/sqrt(N));
    }
    for (i = 0; i < 10; i++)
    {
      fprintf(parallel_dft_out, "%d %e %e\n", i,
      REAL(parallel_dft,i)/sqrt(N),
      IMAG(parallel_dft,i)/sqrt(N));
    }

    start = omp_get_wtime();
    gsl_fft_complex_radix2_forward (complex_polynomial, 1, N);
    end   = omp_get_wtime();
    times[1]  = end-start;

    for (i = 0; i < 10; i++)
    {
      fprintf(serial_fft_out, "%d %e %e\n", i,
      REAL(complex_polynomial,i)/sqrt(N),
      IMAG(complex_polynomial,i)/sqrt(N));
    }
    //Timing info
    fprintf(stdout, "Timing.\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "Serial FFT  : %f\n", times[1]);
    fprintf(stdout, "Parallel DFT: %f\n", times[0]);
  }
  MPI_Finalize();

}
