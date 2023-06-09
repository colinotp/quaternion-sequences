#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <iostream>
#include <tuple>
#include <array>
#include <set>
#include <sys/types.h>
#include <sys/stat.h>
#include "coprimelist.h"
#include "decomps.h"

#define MAX_N 70
#define VERBOSE 0
#define PRINTPAFS 0
#define PRINTTABLE 0

#define NDEBUG

int paf(int n, const std::array<int, MAX_N> &A, int s)
{	int res = 0;
	for(int i=0; i<n; i++)
		res += A[i]*A[(i+s)%n];
	return res;
}

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

void altnegateA(int n, std::array<int, MAX_N> &A)
{	for(int i=1; i<n; i+=2)
		A[i] = -A[i];
}

std::array<int, MAX_N> permuteA(int n, int k, int s, const std::array<int, MAX_N> &A)
{	std::array<int, MAX_N> result = {};
	for(int i=0; i<n; i++)
	{	result[i] = A[(i*k+s)%n];
	}
	return result;
}

std::tuple<std::array<int, MAX_N>, std::array<int, MAX_N>, std::array<int, MAX_N>, std::array<int, MAX_N>> minrep(int n, int l, const std::array<int, MAX_N> &A, const std::array<int, MAX_N> &B, const std::array<int, MAX_N> &C, const std::array<int, MAX_N> &D)
{	
	std::set<std::tuple<std::array<int, MAX_N>, std::array<int, MAX_N>, std::array<int, MAX_N>, std::array<int, MAX_N>>> equivseqns;

	for(int j=0; j<coprimelist_len[n]; j++)
	{	int k = coprimelist[n][j];
		std::array<int, MAX_N> permutedA = permuteA(l, k, 0, A);
		std::array<int, MAX_N> permutedB = permuteA(l, k, 0, B);
		std::array<int, MAX_N> permutedC = permuteA(l, k, 0, C);
		std::array<int, MAX_N> permutedD = permuteA(l, k, 0, D);

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

std::tuple<std::array<int, MAX_N>, std::array<int, MAX_N>, std::array<int, MAX_N>, std::array<int, MAX_N>> minrep_full(int n, int l, std::array<int, MAX_N> &A, std::array<int, MAX_N> &B, std::array<int, MAX_N> &C, std::array<int, MAX_N> &D)
{	
	std::tuple<std::array<int, MAX_N>, std::array<int, MAX_N>, std::array<int, MAX_N>, std::array<int, MAX_N>> rep1 = minrep(n, n, A, B, C, D);
	
	if(n % 2 == 1)
		return rep1;
	else
	{
		altnegateA(n, A);
		altnegateA(n, B);
		altnegateA(n, C);
		altnegateA(n, D);
		
		std::tuple<std::array<int, MAX_N>, std::array<int, MAX_N>, std::array<int, MAX_N>, std::array<int, MAX_N>> rep2 = minrep(n, n, A, B, C, D);

		altnegateA(n, A);
		altnegateA(n, B);
		altnegateA(n, C);
		altnegateA(n, D);
		
		if(rep1 < rep2)
			return rep1;
		else
			return rep2;
	}
}

void fprintseqn(FILE* f, int n, const std::tuple<std::array<int, MAX_N>, std::array<int, MAX_N>, std::array<int, MAX_N>, std::array<int, MAX_N>> &seqn)
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

#if PRINTPAFS==1
int paf(int n, int s, std::array<int, MAX_N> A)
{	int res = 0;
	for(int i=0; i<n; i++)
		res += A[i]*A[(i+s)%n];
	return res;
}

int rowsum(const int n, const std::array<int, MAX_N> &A)
{	int result = 0;
	for(int i=0; i<n; i++)
		result += A[i];
	return result;
}
#endif

int pcf(const int n, const std::array<int, MAX_N> &A, const std::array<int, MAX_N> &B, const int s)
{	int res = 0;
	for(int i=0; i<n; i++)
		res += A[i]*B[(i+s)%n];
	return res;
}

void fprettyprintseqn(FILE* f, int n, const std::tuple<std::array<int, MAX_N>, std::array<int, MAX_N>, std::array<int, MAX_N>, std::array<int, MAX_N>> &seqn)
{	
	std::array<int, MAX_N> A;
	std::array<int, MAX_N> B;
	std::array<int, MAX_N> C;
	std::array<int, MAX_N> D;
	std::tie(A, B, C, D) = seqn;

	for(int i=0; i<n; i++)
		fprintf(f, "%c", A[i] == 1 ? '+' : '-');
	fprintf(f, " ");
	for(int i=0; i<n; i++)
		fprintf(f, "%c", B[i] == 1 ? '+' : '-');
	fprintf(f, " ");
	for(int i=0; i<n; i++)
		fprintf(f, "%c", C[i] == 1 ? '+' : '-');
	fprintf(f, " ");
	for(int i=0; i<n; i++)
		fprintf(f, "%c", D[i] == 1 ? '+' : '-');

	#if PRINTPAFS==1
	fprintf(f, " %d %d %d %d |", rowsum(n, A), rowsum(n, B), rowsum(n, C), rowsum(n, D));

	//altnegateA(n, A);
	//altnegateA(n, B);
	//altnegateA(n, C);
	//altnegateA(n, D);

	//fprintf(f, " %d %d %d %d", rowsum(n, A), rowsum(n, B), rowsum(n, C), rowsum(n, D));
	
	/*fprintf(f, " (");
	for(int i=0; i<n; i++)
		fprintf(f, "%d%s", paf(n, A, i)+paf(n, B, i)+paf(n, C, i)+paf(n, D, i), i<n-1 ? " " : "");
	fprintf(f, ") ");

	fprintf(f, "(");
	for(int i=0; i<n; i++)
		fprintf(f, "%d%s", paf(n, A, i), i<n-1 ? " " : "");
	fprintf(f, ") ");

	fprintf(f, "(");
	for(int i=0; i<n; i++)
		fprintf(f, "%d%s", paf(n, B, i), i<n-1 ? " " : "");
	fprintf(f, ") ");

	fprintf(f, "(");
	for(int i=0; i<n; i++)
		fprintf(f, "%d%s", paf(n, C, i), i<n-1 ? " " : "");
	fprintf(f, ") ");

	fprintf(f, "(");
	for(int i=0; i<n; i++)
		fprintf(f, "%d%s", paf(n, D, i), i<n-1 ? " " : "");
	fprintf(f, ") ");*/

	fprintf(f, " (");
	for(int i=0; i<n; i++)
		fprintf(f, "%d%s", pcf(n, A, B, i)-pcf(n, B, A, i)+pcf(n, C, D, i)-pcf(n, D, C, i), i<n-1 ? " " : "");
	fprintf(f, ") ");

	fprintf(f, "(");
	for(int i=0; i<n; i++)
		fprintf(f, "%d%s", pcf(n, A, B, i), i<n-1 ? " " : "");
	fprintf(f, ") ");

	fprintf(f, "(");
	for(int i=0; i<n; i++)
		fprintf(f, "%d%s", pcf(n, B, A, i), i<n-1 ? " " : "");
	fprintf(f, ") ");

	fprintf(f, "(");
	for(int i=0; i<n; i++)
		fprintf(f, "%d%s", pcf(n, C, D, i), i<n-1 ? " " : "");
	fprintf(f, ") ");

	fprintf(f, "(");
	for(int i=0; i<n; i++)
		fprintf(f, "%d%s", pcf(n, D, C, i), i<n-1 ? " " : "");
	fprintf(f, ") ");
	#endif

	fprintf(f, "\n");
}

int main(int argc, char** argv)
{
	int n1, n2;

	if(argc==1)
		fprintf(stderr, "Need order of exhaust file to remove equivalences from\n"), exit(0);
	else if(argc==2)
	{	n1 = atoi(argv[1]);
		n2 = n1;
	}
	else
	{	n1 = atoi(argv[1]);
		n2 = atoi(argv[2]);
	}

	int result;
	char filename[100];
	#if defined(NOCOMPRESSION)
	const char seqnsfilename[] = "matchedpairs/%d.%d.%d";
	const char seqnsoutfilename[] = "matchedpairs/%d.inequiv";
	#elif defined(ADDER)
	const char seqnsfilename[] = "exhaustadder/%d.%d";
	const char seqnsoutfilename[] = "exhaustadder/%d.inequiv";
	#elif defined(USECOS)
	const char seqnsfilename[] = "exhaustusecos/%d.%d";
	const char seqnsoutfilename[] = "exhaustusecos/%d.inequiv";
	#elif defined(SORT)
	const char seqnsfilename[] = "exhaustsort/%d.%d";
	const char seqnsoutfilename[] = "exhaustsort/%d.inequiv";
	#elif defined(UNCOMP)
	const char seqnsfilename[] = "uncomp/%d.%d";
	const char seqnsoutfilename[] = "uncomp/%d.inequiv";
	#else
	const char seqnsfilename[] = "exhaust/%d.%d";
	const char seqnsoutfilename[] = "exhaust/%d.inequiv";
	#endif
	const char seqnsprettyoutfilename[] = "exhaust/quasiwilliamson-%d.txt";

	mkdir("timings", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

	#if PRINTTABLE==1
	std::cout << "<tr>\n";
	std::cout << "<td>order <var>n</var></td>\n";
	std::cout << "<td># sequences</td>\n";
	std::cout << "<td>filesize</td>\n";
	std::cout << "</tr>\n";
	#endif

	for(int n=n1; n<=n2; n++)
	{
		#if ONLY_EVEN
		if(n%2 == 1 && n1 != n2)
			continue;
		#endif
		//if(n%2 != 0 && n%3 != 0 && n%5 != 0)
		//	continue;

		FILE* seqnsfile, * seqnsoutfile, * seqnsprettyoutfile;

		int invalidcount = 0;
		
		sprintf(filename, seqnsprettyoutfilename, n);
		seqnsprettyoutfile = fopen(filename, "w");

		sprintf(filename, seqnsoutfilename, n);
		seqnsoutfile = fopen(filename, "w");

		int in, totalcount = 0, inequivcount = 0, notquasicount = 0;
		std::set<std::tuple<std::array<int, MAX_N>, std::array<int, MAX_N>, std::array<int, MAX_N>, std::array<int, MAX_N>>> inequivseqns;
		std::array<int, MAX_N> A = {};
		std::array<int, MAX_N> B = {};
		std::array<int, MAX_N> C = {};
		std::array<int, MAX_N> D = {};
	
		clock_t start = clock();
	
		for(int c = 0; c < decomps_len[n]; c++)
		{
			#ifdef NOCOMPRESSION
			sprintf(filename, seqnsfilename, n, c, n);
			#else
			sprintf(filename, seqnsfilename, n, c);
			#endif
			seqnsfile = fopen(filename, "r");
			if(seqnsfile==NULL)
				continue;

			while(fscanf(seqnsfile, "%d ", &in)>0)
			{	
				A[0] = in;
				int res;
				for(int i=1; i<n; i++)
				{	res = fscanf(seqnsfile, "%d ", &in);
					A[i] = in;
				}
				for(int i=0; i<n; i++)
				{	res = fscanf(seqnsfile, "%d ", &in);
					B[i] = in;
				}
				for(int i=0; i<n; i++)
				{	res = fscanf(seqnsfile, "%d ", &in);
					C[i] = in;
				}
				for(int i=0; i<n; i++)
				{	res = fscanf(seqnsfile, "%d ", &in);
					D[i] = in;
				}

				bool tobreak = false;

				for(int s=1; s<=n/2; s++)
				{	if(paf(n, A, s) + paf(n, B, s) + paf(n, C, s) + paf(n, D, s) != 0)
					{	invalidcount++;
						tobreak = true;
						break;
					}
					if(pcf(n, A, B, s) + pcf(n, C, D, s) != pcf(n, B, A, s) + pcf(n, D, C, s))
					{	notquasicount++;
						tobreak = true;
						break;
					}
				}

				totalcount++;

				if(tobreak)
					continue;

				std::tuple<std::array<int, MAX_N>, std::array<int, MAX_N>, std::array<int, MAX_N>, std::array<int, MAX_N>> repseqn = minrep_full(n, n, A, B, C, D);
				if(inequivseqns.count(repseqn)==0)
				{	inequivseqns.insert(repseqn);

					std::array<int, MAX_N> rA;
					std::array<int, MAX_N> rB;
					std::array<int, MAX_N> rC;
					std::array<int, MAX_N> rD;
					std::tie(rA, rB, rC, rD) = repseqn;

					for(int s=1; s<=n/2; s++)
					{	if(paf(n, rA, s) + paf(n, rB, s) + paf(n, rC, s) + paf(n, rD, s) != 0 || pcf(n, rA, rB, s) + pcf(n, rC, rD, s) != pcf(n, rB, rA, s) + pcf(n, rD, rC, s))
						{	fprettyprintseqn(stderr, n, make_tuple(A, B, C, D));
							fprettyprintseqn(stderr, n, make_tuple(rA, rB, rC, rD));
							exit(0);
						}
					}

					fprintseqn(seqnsoutfile, n, make_tuple(A, B, C, D));
					fprettyprintseqn(seqnsprettyoutfile, n, make_tuple(A, B, C, D));
					#if VERBOSE==1
					fprettyprintseqn(stdout, n, make_tuple(A, B, C, D));
					#endif
					inequivcount++;
				}
			}

			fclose(seqnsfile);

			if(invalidcount > 0)
				printf("  WARNING: Read %d sequences which were not Williamson\n", invalidcount);
		}
		fclose(seqnsprettyoutfile);
		fclose(seqnsoutfile);

		#ifdef NOCOMPRESSION
		sprintf(filename, "timings/%d.equivexhausttimenocompress", n);
		#elif defined(UNCOMP)
		sprintf(filename, "timings/%d.equivexhausttimeuncomp", n);
		#else
		sprintf(filename, "timings/%d.equivexhausttime", n);
		#endif
		FILE* f = fopen(filename, "w");
		fprintf(f, "%.2f\n", (clock() - start)/(float)CLOCKS_PER_SEC);
		fclose(f);

		#if PRINTTABLE==0
		printf("Order %d: %d/%d inequivalent sequences (%d not quasi-Williamson) of length %d output in %.2f seconds\n", n, inequivcount, totalcount, notquasicount, n, (clock() - start)/(float)CLOCKS_PER_SEC);
		#else
		FILE *fin;
		char command[512];
		char buff[512];
		sprintf(command, "ls -hl exhaust/williamson-%d.txt | awk '{ print $5; }'", n);
		fin = popen(command, "r");
		while(fgets(buff, sizeof(buff), fin)!=NULL)
		{	buff[strlen(buff)-1] = '\0';
			std::cout << "<tr>\n";
			std::cout << "<td>" << n << "</td>\n";
			//std::cout << "<td><a href=\"williamson-" << n << ".txt\">order " << n << "</a></td>\n";
			std::cout << "<td><a href=\"williamson-" << n << ".txt\">" << inequivcount << "</a></td>\n";
			std::cout << "<td>" << buff << "</td>\n";
			std::cout << "</tr>\n";
		}
		pclose(fin);
		#endif
	}

}
