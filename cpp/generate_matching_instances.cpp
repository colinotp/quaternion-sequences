#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <fftw3.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <array>
#define MAX_N 70
//#define PRINT_PSDS

double* fft_signal;
fftw_complex* fft_result;
fftw_plan plan;

int paf(const int n, const int* A, const int s)
{	int res = 0;
	for(int i=0; i<n; i++)
		res += A[i]*A[(i+s)%n];
	return res;
}

void fprintseqn(FILE* f, const int n, int* A)
{	
	for(int i=0; i<n; i++)
		fprintf(f, "%d ", A[i]);
	fprintf(f, "\n");
}

void fprintpafs(FILE* f, const int n, int* A, const int lineno)
{	
	for(int i=0; i<=n/2; i++)
		fprintf(f, "%d ", paf(n, A, i));
	fprintf(f, ": %d \n", lineno);
}

#ifdef PRINT_PSDS
void fprintpsds(FILE* f, const int n, fftw_complex* P)
{	
	for(int i=0; i<=n/2; i++)
		fprintf(f, "%.5f ", fft_result[i][0]*fft_result[i][0]+fft_result[i][1]*fft_result[i][1]);
	fprintf(f, "\n");
}
#endif

int main(int argc, char** argv)
{
	if(argc<=1)
		fprintf(stderr, "Need order of matchings to compute\n"), exit(0);

	const int n = atoi(argv[1]);

	const char seqnsfilename[] = "matchings/%d.%d.seqns.txt";
	const char psdsfilename[] = "matchings/%d.%d.psds.txt";
	const char pafsfilename[] = "matchings/%d.%d.pafs.txt";

	mkdir("matchings", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	mkdir("timings", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

	printf("ORDER %d: Generate sequences\n", n);

	fft_signal = (double*)malloc(sizeof(double)*n);
	fft_result = (fftw_complex*)malloc(sizeof(fftw_complex)*n);
	plan = fftw_plan_dft_r2c_1d(n, fft_signal, fft_result, FFTW_ESTIMATE);

	FILE *seqnsfile[2*n];
	#ifdef PRINT_PSDS
	FILE *psdsfile[2*n];
	#endif
	FILE *pafsfile[2*n];

	int counts[2*n];
	
	char filename[100];

	for(int i=0; i*i<=4*n; i++)
	{
		if(i%2 != n%2)
			continue;

		sprintf(filename, seqnsfilename, n, i);
		seqnsfile[i] = fopen(filename, "w");
		#ifdef PRINT_PSDS
		sprintf(filename, psdsfilename, n, i);
		psdsfile[i] = fopen(filename, "w");
		#endif
		sprintf(filename, pafsfilename, n, i);
		pafsfile[i] = fopen(filename, "w");
		counts[i] = 0;
	}

	clock_t start = clock();
	int count = 0, psdfiltcount = 0;

	std::array<int, MAX_N> A = {};
	for(int i=0; i<n; i++)
	{	A[i] = -1;
		fft_signal[i] = A[i];
	}

	int rowsum = -n;

	while(1)
	{	
		bool filtered = false;
		
		if(rowsum < 0)
			filtered = true;
		else
		{
			fftw_execute(plan);

			for(int i=0; i<n; i++)
			{	const double psd_i = fft_result[i][0]*fft_result[i][0] + fft_result[i][1]*fft_result[i][1];
				if(psd_i > (4*n + 0.001))
				{	filtered = true;
					psdfiltcount++;
					break;
				}
			}
		}

		if(!filtered)
		{
			fprintseqn(seqnsfile[rowsum], n, A.data());
			#ifdef PRINT_PSDS
			fprintpsds(psdsfile[rowsum], n, fft_result);
			#endif
			fprintpafs(pafsfile[rowsum], n, A.data(), count);
			count++;
			counts[rowsum]++;
		}

		int i;
		for(i=0; i<n; i++)
		{	A[i] += 2;
			rowsum += 2;
			fft_signal[i] = A[i];
			if(A[i]==3)
			{	A[i] = -1;
				rowsum -= 4;
				fft_signal[i] = A[i];
			}
			else
				break;
		}
		if(i==n)
			break;
	}

	sprintf(filename, "timings/%d.gentime", n);
	FILE* f = fopen(filename, "w");
	fprintf(f, "%.2f\n", (clock() - start)/(float)CLOCKS_PER_SEC);
	fclose(f);

	printf("  %d sequences of length %d generated in %.2f seconds, %d sequences filtered using PSD test\n", count, n, (clock() - start)/(float)CLOCKS_PER_SEC, psdfiltcount);

	for(int i=0; i*i<=4*n; i++)
	{
		if(i%2 != n%2)
			continue;

		fclose(seqnsfile[i]);
		#ifdef PRINT_PSDS
		fclose(psdsfile[i]);
		#endif
		fclose(pafsfile[i]);
		printf("  rowsum %d count: %d\n", i, counts[i]);
	}

	fftw_destroy_plan(plan);
	free(fft_signal);
	free(fft_result);

}
