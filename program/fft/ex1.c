#include <fftw3.h>
#include <stdio.h>
#include <math.h>
int main(){
  const int N=500;
  const double PI = acos(-1.0);
  double f_0 = 10.0;
  
  int i;
  double *in;
  fftw_complex *out;
  fftw_plan p;
  in = (double*) fftw_malloc(sizeof(double) * N);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N); // just need N/2+1
  p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_MEASURE);


  // Initialization & write the input
  FILE *infile = fopen("inData.dat", "w");
  for(i=0;i<N;i++){
    if(i<N/2) in[i] = 0.0;
	 else
       in[i] = 1.0;

	out[i][0] = 0.0;
	out[i][1] = 0.0;
    fprintf(infile, "%f\t%f\n", in[i], 0.0);
  }
  fclose(infile);


  fftw_execute(p); /* repeat as needed */


  // write the output
  FILE *outfile = fopen("outData.dat", "w");
  for(i=0;i<N;i++){
    fprintf(outfile, "%f\t%f\n", out[i][0], out[i][1]);
  }
  fclose(outfile);



  fftw_destroy_plan(p);
  fftw_free(in); fftw_free(out);
}
