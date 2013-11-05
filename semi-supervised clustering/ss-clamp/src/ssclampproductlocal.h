#ifndef SS_CLAMP_PRODUCT_LOCAL_H_
#define SS_CLAMP_PRODUCT_LOCAL_H_

#include "ssclamp.h"

// this is a very reasonable value for the maximum relevance
#define MAX_RELEVANCE (1e2)

struct SSClampProductLocal : public SSClamp {

  SSClampProductLocal(Parameter &params) : SSClamp(params) {
  }

  virtual void update_relevance(Answer& answer) {
    for(unsigned int k = 0; k < params.C; ++k) {
      // initially the relevance for each table can be updated
      vector<bool> can_update(params.T, true);
      unsigned int count = params.T;
      // we keep this loop while there's a bad relevance value
      for(bool bad = true; bad;) {
        double num = 1.0;
        vector<double> denominator(params.T,0.0);
        for(unsigned int t = 0; t < params.T; ++t) {
          if(can_update[t]) {
            for(unsigned int i = 0; i < params.N; ++i) {
              if(params.mask[i]) {
                for(unsigned int p = 0; p < params.P; ++p) {
                  denominator[t] += params.table[t][i][answer.prototype[k][p]] * pow(answer.U[i][k], 2.0);
                }
              }
            }
            num *= pow(denominator[t], 1.0 / count);
          } else {
            VALIDATE_DENOMINATOR(answer.Relevance[k][t]);
            num /= pow(answer.Relevance[k][t], 1.0 / count);
          }
        }
        bad = false;
        // check if the new values are good
        for(unsigned int t = 0; t < params.T; ++t) {
          if(can_update[t]) {
            if(Util::cmp(denominator[t]) <= 0 or Util::cmp(num / denominator[t], MAX_RELEVANCE) >= 0) {
              bad = true;
              can_update[t] = 0;
              count--;
            }
          }
        }
        if(!bad) {
          // update relevance
          for(unsigned int t = 0; t < params.T; ++t) {
            if(can_update[t]) {
              VALIDATE_DENOMINATOR(denominator[t]);
              answer.Relevance[k][t] = num / denominator[t];
            }
          }
        }
      }
    }
  }

};

#undef MAX_RELEVANCE

#endif
