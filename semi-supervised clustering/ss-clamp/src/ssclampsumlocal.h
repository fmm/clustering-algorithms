#ifndef SS_CLAMP_SUM_LOCAL_H_
#define SS_CLAMP_SUM_LOCAL_H_

#include "ssclamp.h"

struct SSClampSumLocal : public SSClamp {

  SSClampSumLocal(Parameter &params) : SSClamp(params) {
  }

  virtual void set_default_relevance(Answer &answer) {
    for(unsigned int k = 0; k < params.C; ++k) {
      for(unsigned int t = 0; t < params.T; ++t) {
        answer.Relevance[k][t] = 1.0 / params.T;
      }
    }
  }

  virtual inline double dissimilarity(Answer &answer, unsigned int i, unsigned int j, unsigned int k, unsigned int t) {
    return params.table[t][i][j] * pow(answer.Relevance[k][t], params.relevance_v);
  }

  virtual void update_relevance(Answer& answer) {
    vector< vector<double> > dist(params.C,vector<double>(params.T,0.0));
    for(unsigned int t = 0; t < params.T; ++t) {
      for(unsigned int k = 0; k < params.C; ++k) {
        for(unsigned int i = 0; i < params.N; ++i) {
          for(unsigned int p = 0; p < params.P; ++p) {
            dist[k][t] += params.table[t][i][answer.prototype[k][p]] * pow(answer.U[i][k], 2.0);
          }
        }
      }
    }
    for(unsigned int k = 0; k < params.C; ++k) {
      for(unsigned int t = 0; t < params.T; ++t) {
        VALIDATE_DENOMINATOR(dist[k][t]);
      }
    }
    for(unsigned int k = 0; k < params.C; ++k) {
      for(unsigned int t = 0; t < params.T; ++t) {
        answer.Relevance[k][t] = 0.0;
        for(unsigned int h = 0; h < params.T; ++h) {
          answer.Relevance[k][t] += pow(dist[k][t] / dist[k][h], 1.0 / (params.relevance_v - 1));
        }
        VALIDATE_DENOMINATOR(answer.Relevance[k][t]);
        answer.Relevance[k][t] = 1.0 / answer.Relevance[k][t];
      }
    }		
  }

};

#endif
