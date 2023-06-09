#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <tuple>
#include <array>
#include <set>
#include <sys/types.h>
#include <sys/stat.h>
#include "coprimelist.h"
#include "decomps.h"

#define MAX_N 70

void swap(std::array<int, MAX_N> &A, std::array<int, MAX_N> &B)
{	std::array<int, MAX_N> tmp;
	tmp = A;
	A = B;
	B = tmp;
}

void negateA(int n, std::array<int, MAX_N> &A)
{	for(int i=0; i<n; i++)
		A[i] = -A[i];
}

int rowsum(int n, std::array<int, MAX_N> A)
{	int result;
	for(int i=0; i<n; i++)
		result += A[i];
	return result;
}

int firstnonzero(int n, std::array<int, MAX_N> A)
{	int i;
	for(i=0; i<n; i++)
		if(A[i]!=0)
			return A[i];
	return A[n-1];
}

std::array<int, MAX_N> permuteA(int n, int k, std::array<int, MAX_N> A)
{	std::array<int, MAX_N> result = {};
	for(int i=0; i<n; i++)
	{	result[i] = A[(i*k)%n];
	}
	return result;
}

std::tuple<std::array<int, MAX_N>, std::array<int, MAX_N>, std::array<int, MAX_N>, std::array<int, MAX_N>> minrep(int n, int l, std::array<int, MAX_N> A, std::array<int, MAX_N> B, std::array<int, MAX_N> C, std::array<int, MAX_N> D)
{	
	std::set<std::tuple<std::array<int, MAX_N>, std::array<int, MAX_N>, std::array<int, MAX_N>, std::array<int, MAX_N>>> equivseqns;

	for(int j=0; j<coprimelist_len[n]; j++)
	{	int k = coprimelist[n][j];
		std::array<int, MAX_N> permutedA = permuteA(l, k, A);
		std::array<int, MAX_N> permutedB = permuteA(l, k, B);
		std::array<int, MAX_N> permutedC = permuteA(l, k, C);
		std::array<int, MAX_N> permutedD = permuteA(l, k, D);

		if(permutedA[0]<0)
		{	negateA(l, permutedA);
			negateA(l, permutedB);
		}

		if(permutedB[0]<0)
		{	negateA(l, permutedB);
			negateA(l, permutedD);
		}

		if(permutedC[0]<0)
		{	negateA(l, permutedC);
			negateA(l, permutedD);
		}

		/*if(permutedA>permutedB)
		{	swap(permutedA, permutedB);
			swap(permutedC, permutedD);
		}
		if(permutedA>permutedC)
		{	swap(permutedA, permutedC);
			swap(permutedB, permutedD);
		}*/

		equivseqns.insert(make_tuple(permutedA, permutedB, permutedC, permutedD));
	}

	return *(equivseqns.begin());
}

void fprintseqn(FILE* f, int n, std::tuple<std::array<int, MAX_N>, std::array<int, MAX_N>, std::array<int, MAX_N>, std::array<int, MAX_N>> seqn)
{	
	std::array<int, MAX_N> A;
	std::array<int, MAX_N> B;
	std::array<int, MAX_N> C;
	std::array<int, MAX_N> D;
	std::tie(A, B, C, D) = seqn;

	for(int i=0; i<n; i++)
		fprintf(f, "%d ", A[i]);
	for(int i=0; i<n; i++)
		fprintf(f, "%d ", B[i]);
	for(int i=0; i<n; i++)
		fprintf(f, "%d ", C[i]);
	for(int i=0; i<n; i++)
		fprintf(f, "%d ", D[i]);
	fprintf(f, "\n");
}

/*void fprintseqn(FILE* f, int n, int* A)
{	
	for(int i=0; i<n; i++)
		fprintf(f, "%d ", A[i]);
	fprintf(f, "\n");
}*/

int pcf(const int n, const std::array<int, MAX_N> &A, const std::array<int, MAX_N> &B, const int s)
{	int res = 0;
	for(int i=0; i<n; i++)
		res += A[i]*B[(i+s)%n];
	return res;
}

int main(int argc, char** argv)
{
	int n1, n2;

	if(argc==1)
		fprintf(stderr, "Need order of matched files to remove equivalences from and optionally the case # to solve and the compression factor\nPassing -1 as the case # will solve all cases\n"), exit(0);

	const int n = atoi(argv[1]);

	int casetosolve = -1;
	if(argc>2)
	{	casetosolve = atoi(argv[2]);
	}

	int d;
	if(argc>3)
	{	d = atoi(argv[3]);
	}
	else
	{	d = 1;
	}
	const int l = n/d;

	int result;
	char filename[100];
	const char seqnsfilename[] = "matchedpairs/%d.%d.%d";
	const char seqnsoutfilename[] = "matchedseqns/%d.%d.%d.inequiv";

	mkdir("timings", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	mkdir("matchedseqns", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

	printf("ORDER %d: Remove equivalent matched compressions\n", n);
	
	FILE* seqnsfile, * seqnsoutfile;

	for(int c = 0; c < decomps_len[n]; c++)
	{
		if(casetosolve != -1 && casetosolve != c)
			continue;

		sprintf(filename, seqnsoutfilename, n, c, l);
		seqnsoutfile = fopen(filename, "w");

		clock_t start = clock();

		int in, totalcount = 0, inequivcount = 0;
		std::set<std::tuple<std::array<int, MAX_N>, std::array<int, MAX_N>, std::array<int, MAX_N>, std::array<int, MAX_N>>> inequivseqns;
		std::array<int, MAX_N> A = {};
		std::array<int, MAX_N> B = {};
		std::array<int, MAX_N> C = {};
		std::array<int, MAX_N> D = {};

		sprintf(filename, seqnsfilename, n, c, l);
		seqnsfile = fopen(filename, "r");

		while(fscanf(seqnsfile, "%d ", &in)>0)
		{	
			A[0] = in;
			int res;
			for(int i=1; i<l; i++)
			{	res = fscanf(seqnsfile, "%d ", &in);
				A[i] = in;
			}
			for(int i=0; i<l; i++)
			{	res = fscanf(seqnsfile, "%d ", &in);
				B[i] = in;
			}
			for(int i=0; i<l; i++)
			{	res = fscanf(seqnsfile, "%d ", &in);
				C[i] = in;
			}
			for(int i=0; i<l; i++)
			{	res = fscanf(seqnsfile, "%d ", &in);
				D[i] = in;
			}

			// Loop twice, the second time negating A
			for(int neg=0; neg<2; neg++)
			{	bool to_break = false;

				for(int k=0; k<n; k++)
				{	const int pAB = pcf(l, A, B, k);
					const int pBA = pcf(l, B, A, k);
					const int pCD = pcf(l, C, D, k);
					const int pDC = pcf(l, D, C, k);
					if(pAB - pBA + pCD - pDC != 0)
					{	to_break = true;
						break;
					}
				}

				if(!to_break) {
					std::tuple<std::array<int, MAX_N>, std::array<int, MAX_N>, std::array<int, MAX_N>, std::array<int, MAX_N>> repseqn = minrep(n, l, A, B, C, D);
					if(inequivseqns.count(repseqn)==0)
					{	inequivseqns.insert(repseqn);
						fprintseqn(seqnsoutfile, l, make_tuple(A, B, C, D));
						inequivcount++;
					}

					totalcount++;
				}

				// Need to also negate one of the sequences for a complete search, as negating a single sequence is not a equivalence of Quasi-Williamson sequences
				negateA(l, D);
			}
		}

		fclose(seqnsfile);

		#ifdef NOCOMPRESSION
		sprintf(filename, "timings/%d.%d.%d.equivpairstimenocompress", n, c, l);
		#else
		sprintf(filename, "timings/%d.%d.equivpairstime", n, c);
		#endif
		FILE* f = fopen(filename, "w");
		fprintf(f, "%.2f\n", (clock() - start)/(float)CLOCKS_PER_SEC);
		fclose(f);

		printf("  Case %d: %d/%d inequivalent matched sequences of length %d output in %.2f seconds\n", c, inequivcount, totalcount, l, (clock() - start)/(float)CLOCKS_PER_SEC);

		fclose(seqnsoutfile);
	}

}
