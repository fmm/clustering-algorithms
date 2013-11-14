#ifndef SS_CLAMP_PRODUCT_GLOBAL_H_
#define SS_CLAMP_PRODUCT_GLOBAL_H_

#include "ssclamp.h"

struct SSClampProductGlobal : public SSClamp {

  SSClampProductGlobal(Parameter &params) : SSClamp(params) {
  }

  virtual void update_relevance(Answer& answer) {
    ASSERT(Util::cmp(answer.beta) == 0, "unimplemented");
    double num = 1.0;
    vector<double> denominator(params.T,0.0);
    for(unsigned int t = 0; t < params.T; ++t) {
      for(unsigned int k = 0; k < params.C; ++k) {
        for(unsigned int i = 0; i < params.N; ++i) {
          if(params.mask[i]) {
            for(unsigned int p = 0; p < params.P; ++p) {
              denominator[t] += params.table[t][i][answer.prototype[k][p]] * pow(answer.U[i][k], 2.0);
            }
          }
        }
      }
      VALIDATE_DENOMINATOR(denominator[t]);
      num *= pow(denominator[t], 1.0 / params.T);
    }
    for(unsigned int k = 0; k < params.C; ++k) {
      for(unsigned int t = 0; t < params.T; ++t) {
        answer.Relevance[k][t] = num / denominator[t];
      }
    }		
  }

};

#endif
