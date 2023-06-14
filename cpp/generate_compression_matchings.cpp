#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <array>
#include <set>
#include <map>

#define MAX_N 70

int paf(int n, const int* A, int s)
{	int res = 0;
	for(int i=0; i<n; i++)
		res += A[i]*A[(i+s)%n];
	return res;
}

void printA(int n, const int* A)
{	for(int i=0; i<n; i++)
		printf("%d ", A[i]);
	printf(": ");
	for(int i=0; i<n; i++)
		printf("%d ", paf(n, A, i));
	printf("\n");
}

void fprintseqn(FILE* f, int n, const int* A)
{	
	for(int i=0; i<n; i++)
		fprintf(f, "%d ", A[i]);
	fprintf(f, "\n");
}

void fprintpafs(FILE* f, int n, const int* A)
{	
	for(int i=0; i<=n/2; i++)
		fprintf(f, "%d ", paf(n, A, i));
	fprintf(f, "\n");
}

std::array<int, MAX_N> compress(int n, int l, std::array<int, MAX_N> A)
{	
	std::array<int, MAX_N> result = {};
	for(int i=0; i<n; i++)
	{	result[i%l] += A[i];
	}
	return result;
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

	const char seqnsfilename[] = "matchings/%d.%d.seqns.txt";
	const char pafsfilename[] = "matchings/%d.%d.%d.pafs.txt";
	const char compfilename[] = "matchings/%d.%d.%d.comp.txt";
	const char mapfilename[] = "matchings/%d.%d.%d.map.txt";

	mkdir("matchings", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	mkdir("timings", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

	printf("ORDER %d: Generate compressions\n", n);

	char filename[100];

	FILE *seqnsfile[2*n];
	FILE *pafsfile[2*n];
	FILE *compfile[2*n];
	FILE *mapfile[2*n];

	for(int i=0; i*i<=4*n; i++)
	{
		if(i%2 != n%2)
			continue;

		int totalcount = 0;
		int compcount = 0;
		int in;

		sprintf(filename, seqnsfilename, n, i);
		seqnsfile[i] = fopen(filename, "r");

		sprintf(filename, pafsfilename, n, i, l);
		pafsfile[i] = fopen(filename, "w");
		sprintf(filename, compfilename, n, i, l);
		compfile[i] = fopen(filename, "w");
		sprintf(filename, mapfilename, n, i, l);
		mapfile[i] = fopen(filename, "w");

		clock_t start = clock();

		std::set<std::array<int, MAX_N>> myset;
		std::multimap<std::array<int, MAX_N>, int> compmap;
		std::array<int, MAX_N> A = {};

		int k = 0;
		while(fscanf(seqnsfile[i], "%d ", &in)>0)
		{	
			A[k] = in;
			k++;
			if(k%n==0)
			{	k = 0;
				std::array<int, MAX_N> compressA = compress(n, l, A);

				compmap.insert(std::pair<std::array<int, MAX_N>, int>(compressA, totalcount));

				if(myset.count(compressA)==0)
				{	
					myset.insert(compressA);
					compcount++;
				}
				totalcount++;
			}
		}

		for(auto it=myset.cbegin(); it!=myset.cend(); it++)
		{	fprintseqn(compfile[i], l, it->data());
			fprintpafs(pafsfile[i], l, it->data());

			std::pair<std::multimap<std::array<int, MAX_N>, int>::iterator, std::multimap<std::array<int, MAX_N>, int>::iterator> ret;
			ret = compmap.equal_range(*it);
			for(auto it=ret.first; it!=ret.second; it++)
			{	fprintf(mapfile[i], "%d ", it->second);
			}
			fprintf(mapfile[i], "\n");
		}

		printf("  rowsum %d: %d sequences read and %d compressions of length %d output in %.2f seconds\n", i, totalcount, compcount, l, (clock() - start)/(float)CLOCKS_PER_SEC);

		sprintf(filename, "timings/%d.%d.%d.comptime", n, i, l);
		FILE* f = fopen(filename, "w");
		fprintf(f, "%.2f\n", (clock() - start)/(float)CLOCKS_PER_SEC);
		fclose(f);

		fclose(seqnsfile[i]);
		fclose(pafsfile[i]);
		fclose(compfile[i]);
		fclose(mapfile[i]);
	}

}
