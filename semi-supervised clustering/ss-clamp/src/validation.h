#ifndef VALIDATION_H_
#define VALIDATION_H_

#include "minimumcostmaximumflow.h"

namespace Validation {
	pair< double, vector<int> > accuracy(vector< vector<unsigned int> >& tab) {
		unsigned int N = tab.size() - 1, M = tab[0].size() - 1, source = 0, sink = N + M + 1;
		vector<int> match(N, -1);
		MinimumCostMaximumFlow mcmf;
		for(unsigned int i = 0; i < N; ++i) {
			mcmf.add(source,i+1,0,1);
		}
		int cflow = (N == M ? 1 : N);
		for(unsigned int j = 0; j < M; ++j) {
			mcmf.add(j+N+1,sink,0,cflow);
		}
		for(unsigned int i = 0; i < N; ++i) {
			for(unsigned int j = 0; j < M; ++j) {
				int w = tab[i][j];
				mcmf.add(i+1,j+N+1,-w,1);
			}
		}
		pair<int,int> ret = mcmf.process(source,sink);
		for(unsigned int i = 0; i < N; ++i) {
			for(unsigned int j = 0; j < M; ++j) {
				if(mcmf.cap[(i*M+j)<<1] != cflow) {
					match[i] = j;
				}
			}
		}
		ASSERT(ret.second == (int)(N), "invalid flow");
		return make_pair((double)(-ret.first) / tab[N][M], match);
	}

	unsigned int comb(unsigned int x) {
		return x * (x - 1) / 2;
	}

	double adjusted_rand_index(vector< vector<unsigned int> >& tab) {
		unsigned int N = tab.size() - 1, M = tab[0].size() - 1;
		double term[4] = {0}, temp[2] = {0};
		double pot = 1.0 / comb(tab[N][M]);
		for(unsigned int i = 0; i < N; ++i) {
			for(unsigned int j = 0; j < M; ++j) {
				term[1] += comb(tab[i][j]);
			}
		}
		for(unsigned int i = 0; i < N; ++i) {
			temp[0] += comb(tab[i][M]);
		}
		for(unsigned int j = 0; j < M; ++j) {
			temp[1] += comb(tab[N][j]);
		}
		term[2] = pot * (temp[0] * temp[1]);
		term[3] = 0.5 * (temp[0] + temp[1]);
		return (term[1] - term[2]) / (term[3] - term[2]);
	}

	double f_measure(vector< vector<unsigned int> >& tab) {
		double F = 0;
		unsigned int N = tab.size() - 1, M = tab[0].size() - 1;
		for(unsigned int j = 0; j < M; ++j) {
			double vmax = 0;
			for(unsigned int i = 0; i < N; ++i) {
				if(tab[i][j] != 0) {
					double rappel = (double)(tab[i][j]) / tab[N][j];
					double precision = (double)(tab[i][j]) / tab[i][M];
					vmax = max(vmax, 2 * rappel * precision / (rappel + precision));
				}
			}
			F += vmax * tab[N][j];
		}
		return F / tab[N][M];
	}
};

#endif
