#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <fftw3.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "decomps.h"

const int MAX_N = 70;
const int HALF_MAX_N = 1+MAX_N/2;
#include <array>
#include <vector>

void fprintpair(FILE* f, int n, int* A, int iA, int iB)
{	for(int i=0; i<n; i++)
		fprintf(f, "%d ", A[i]);
	fprintf(f, ": %d %d\n", iA, iB);
}

int main(int argc, char** argv)
{
	if(argc<=1)
		fprintf(stderr, "Need order of matchings to compute and optionally the compression factor\n"), exit(0);

	const int n = atoi(argv[1]);

	int d;
	if(argc>2)
	{	d = atoi(argv[2]);
	}
	else
	{	d = 1;
	}
	const int l = n/d;
	const int pafslen = l/2+1;

	const char seqnsfilename[] = "matchings/%d.%d.%d.comp.txt";
	const char pafsinfilename[] = "matchings/%d.%d.%d.pafs.txt";
	const char pafsoutfilename[] = "matchings/%d.%d.%d.%s.pafs.txt";

	mkdir("matchings", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	mkdir("timings", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

	printf("ORDER %d: Generate compression sequence pairs\n", n);

	double* fft_signal = (double*)malloc(sizeof(double)*l);
	fftw_complex* fft_result = (fftw_complex*)malloc(sizeof(fftw_complex)*l);
	fftw_plan plan = fftw_plan_dft_r2c_1d(l, fft_signal, fft_result, FFTW_ESTIMATE);

	char filename[100];
	FILE* seqnsfile;
	FILE* pafsfile;

	std::vector<std::array<int, MAX_N>> seqns[2*n];
	std::vector<std::array<int, HALF_MAX_N>> pafs[2*n];
	std::vector<std::array<double, HALF_MAX_N>> psdslist[2*n];

	std::array<int, MAX_N> X = {};
	std::array<int, HALF_MAX_N> P = {};
	std::array<double, HALF_MAX_N> psds = {};

	int in;

	for(int i=0; i*i<=4*n; i++)
	{
		if(i%2 != n%2)
			continue;

		sprintf(filename, seqnsfilename, n, i, l);
		if((seqnsfile = fopen(filename, "r")) == NULL)
			fprintf(stderr, "Could not read file %s\n", filename), exit(0);

		clock_t start = clock();

		int count = 0;
		int k = 0;
		while(fscanf(seqnsfile, "%d ", &in)>0)
		{	X[k] = in;
			k++;
			if(k%l==0)
			{	k = 0;

				for(int j=0; j<l; j++)
					fft_signal[j] = X[j];
				fftw_execute(plan);
				for(int j=1; j<=l/2; j++)
					psds[j] = fft_result[j][0]*fft_result[j][0]+fft_result[j][1]*fft_result[j][1];

				seqns[i].push_back(X);
				psdslist[i].push_back(psds);
				count++;
			}
		}

		printf("  Computed %d rowsum %d PSDs in %.2f seconds\n", count, i, (clock() - start)/(float)CLOCKS_PER_SEC);

		fclose(seqnsfile);

		sprintf(filename, pafsinfilename, n, i, l);
		if((pafsfile = fopen(filename, "r")) == NULL)
			fprintf(stderr, "Could not read file %s\n", filename), exit(0);

		k = 0;
		while(fscanf(pafsfile, "%d ", &in)>0)
		{	P[k] = in;
			k++;
			if(k%pafslen==0)
			{	k = 0;
				pafs[i].push_back(P);
			}
		}

		fclose(pafsfile);
	}

	
	for(int c = 0; c < decomps_len[n]; c++)
	{
		const std::vector<std::array<int, MAX_N>> Aseqns = seqns[decomps[n][c][0]];
		const std::vector<std::array<int, MAX_N>> Bseqns = seqns[decomps[n][c][1]];
		const std::vector<std::array<int, MAX_N>> Cseqns = seqns[decomps[n][c][2]];
		const std::vector<std::array<int, MAX_N>> Dseqns = seqns[decomps[n][c][3]];
		const std::vector<std::array<double, HALF_MAX_N>> Apsdslist = psdslist[decomps[n][c][0]];
		const std::vector<std::array<double, HALF_MAX_N>> Bpsdslist = psdslist[decomps[n][c][1]];
		const std::vector<std::array<double, HALF_MAX_N>> Cpsdslist = psdslist[decomps[n][c][2]];
		const std::vector<std::array<double, HALF_MAX_N>> Dpsdslist = psdslist[decomps[n][c][3]];
		const std::vector<std::array<int, HALF_MAX_N>> Apafs = pafs[decomps[n][c][0]];
		const std::vector<std::array<int, HALF_MAX_N>> Bpafs = pafs[decomps[n][c][1]];
		const std::vector<std::array<int, HALF_MAX_N>> Cpafs = pafs[decomps[n][c][2]];
		const std::vector<std::array<int, HALF_MAX_N>> Dpafs = pafs[decomps[n][c][3]];

		clock_t start = clock();

		long ABcount = 0;
		long CDcount = 0;

		std::array<int, HALF_MAX_N> ABpafs;
		std::array<double, HALF_MAX_N> ABpsds = {};
		std::array<int, HALF_MAX_N> CDpafs;
		std::array<double, HALF_MAX_N> CDpsds = {};
		
		FILE* pairfile;

		sprintf(filename, pafsoutfilename, n, c, l, "AB");
		if((pairfile = fopen(filename, "w")) == NULL)
			fprintf(stderr, "Could not create file %s\n", filename), exit(0);

		for(int iA = 0; iA < Apsdslist.size(); iA++)
		{	
			const std::array<double, HALF_MAX_N> Apsds = Apsdslist[iA];

			for(int iB = 0; iB < Bpsdslist.size(); iB++)
			{	
				const std::array<double, HALF_MAX_N> Bpsds = Bpsdslist[iB];

				bool tobreak = false;
				for(int j=1; j<=l/2; j++)
				{	ABpsds[j] = Apsds[j]+Bpsds[j];
					if(ABpsds[j] > 4*n+0.01)
					{	tobreak = true;
						break;
					}
				}
				if(tobreak)
					continue;

				for(int i = 0; i < pafslen; i++)
					ABpafs[i] = Apafs[iA][i] + Bpafs[iB][i];

				fprintpair(pairfile, pafslen, ABpafs.data(), iA, iB);
				ABcount++;
			}
		}

		fclose(pairfile);

		sprintf(filename, pafsoutfilename, n, c, l, "CD");
		if((pairfile = fopen(filename, "w")) == NULL)
			fprintf(stderr, "Could not create file %s\n", filename), exit(0);

		for(int iC = 0; iC < Cpsdslist.size(); iC++)
		{	
			const std::array<double, HALF_MAX_N> Cpsds = Cpsdslist[iC];

			for(int iD = 0; iD < Dpsdslist.size(); iD++)
			{	
				const std::array<double, HALF_MAX_N> Dpsds = Dpsdslist[iD];

				bool tobreak = false;
				for(int j=1; j<=l/2; j++)
				{	CDpsds[j] = Cpsds[j]+Dpsds[j];
					if(CDpsds[j] > 4*n+0.01)
					{	tobreak = true;
						break;
					}
				}
				if(tobreak)
					continue;

				for(int i = 0; i < pafslen; i++)
				{	CDpafs[i] = - (Cpafs[iC][i] + Dpafs[iD][i]);
					if(i==0)
						CDpafs[i] += 4*n;
				}

				fprintpair(pairfile, pafslen, CDpafs.data(), iC, iD);
				CDcount++;
			}
		}

		fclose(pairfile);

		sprintf(filename, "timings/%d.%d.%d.genpairtime", n, c, l);
		FILE* f = fopen(filename, "w");
		fprintf(f, "%.2f\n", (clock() - start)/(float)CLOCKS_PER_SEC);
		fclose(f);

		printf("  Case %d: %ld AB and %ld CD paired PAF sequences of length %d generated in %.2f seconds\n", c, ABcount, CDcount, pafslen, (clock() - start)/(float)CLOCKS_PER_SEC);

	}

	fftw_destroy_plan(plan);
	free(fft_signal);
	free(fft_result);

}
