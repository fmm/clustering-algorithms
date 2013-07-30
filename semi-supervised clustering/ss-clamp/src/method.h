#ifndef METHOD_H_
#define METHOD_H_

#include <sys/time.h>
#include "database.h"
#include "mersennetwister.h"
#include "parameter.h"

struct Method {

	Parameter &params;

	struct Timer {
		struct timeval start, now;
		void start_timer() {
			gettimeofday(&start, NULL);
		}
		unsigned long long elapsed() {
			gettimeofday(&now, NULL);
			return (now.tv_sec-start.tv_sec)*1000000ULL+(now.tv_usec-start.tv_usec);
		}
	} timer;

	MersenneTwister random;

	struct Answer {
		// global
		unsigned int initialization;
		unsigned int iteration;
		double criterion;
		double restriction;
		vector<Cluster> cluster;
		Matrix U;
		Matrix Relevance;
		// auxiliar
		vector<Prototype> prototype;
		double beta;
		double eps;
		// compare function
		bool operator<(const Answer& answer) const {
			return Util::cmp(criterion, answer.criterion, eps) < 0;
		}
	};

	Method(Parameter &params) : params(params) {
		timer.start_timer();
		random = MersenneTwister(params.seed);
	}

	virtual void initialize(Answer &answer, unsigned int init, unsigned int iter) {
		answer.initialization = init;
		answer.iteration = iter;
		answer.criterion = INF;
		answer.restriction = INF;
		answer.cluster = vector<Cluster>(params.C);
		answer.U = Matrix(params.N,Row(params.C));
		answer.Relevance = Matrix(params.C,Row(params.T));
		answer.eps = params.eps_for_criterion;
	}

	virtual double compute_criterion(Answer &answer) = 0;

	virtual void srand(Answer &answer) = 0;

	virtual bool optimize(Answer &answer) {
		update_clusters(answer);
		double old_criterion = answer.criterion;
		answer.criterion = compute_criterion(answer);
		double new_criterion = answer.criterion;
		// TODO: change to throw exception		
		ASSERT(
				new_criterion <= old_criterion or log(fabs(new_criterion - old_criterion)) <= answer.eps,
				"failed to minimize criterion"
				);
		answer.criterion = new_criterion;
		if(fabs(new_criterion - old_criterion) > answer.eps) {
			// TODO: store results
			return true;
		}
		// local optimum was reached
		return false;
	}

	const void update_clusters(Answer &answer) {
		for(unsigned int k = 0; k < params.C; ++k) {
			answer.cluster[k].clear();
		}
		for(unsigned int i = 0; i < params.N; ++i) {
			vector< pair<double,int> > v(params.C);
			for(unsigned int k = 0; k < params.C; ++k) {
				v[k] = make_pair(answer.U[i][k],k);
			}
			answer.cluster[max_element(v.begin(),v.end())->second].insert(i);
		}
	}

	const vector< vector<unsigned int> > compute_confusion_matrix(Answer &answer) {
		unsigned int k = params.C;
		unsigned int p = params.priori_cluster.size();
		vector< vector<unsigned int> > table(k+1,vector<unsigned int>(p+1,0));
		for(unsigned int i = 0; i < k; ++i) {
			for(unsigned int j = 0; j < p; ++j) {
				for(Cluster::const_iterator iter = answer.cluster[i].begin(); iter != answer.cluster[i].end(); iter++) {
					table[i][j] += params.priori_cluster[j].count(*iter);
				}
			}
		}
		for(unsigned int i = 0; i < k; ++i) {
			table[i][p] = answer.cluster[i].size();
		}
		for(unsigned int j = 0; j < p; ++j) {
			table[k][j] = params.priori_cluster[j].size();
		}
		table[k][p] = params.N;
		return table;
	}

	// TODO: considering that every values is already defined
	Answer process() {
		Answer best;
		initialize(best,0,0);
		for(unsigned int init = 1; init <= params.initialization; ++init) {
			if(timer.elapsed() > params.time_limit) {
				WARNING("time limit was reached at initialization #" + Util::cast<string>(init));
				break;
			}
			Answer now;
			initialize(now, init, 0);
			srand(now);
			for(unsigned int iter = 1; iter <= params.maximum_iteration; ++iter) {
				if(optimize(now)) {
					// TODO: print line
				} else {
					break;
				}
			}
			// optimize the best result
			best = min(best, now);
		}
		// TODO: print result in a pdf file
		return best;
	}

};

#endif
