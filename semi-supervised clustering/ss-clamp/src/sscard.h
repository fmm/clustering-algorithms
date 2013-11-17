#ifndef SS_CARD_H_
#define SS_CARD_H_

#include "method.h"

struct SSCARD : public Method {

  const static int Q = 2;

  // constructor
  SSCARD(Parameter &params) : Method(params) {
    // objetive function doesn't minimize all the time
    check_before_update = false;
  }

  virtual inline double dissimilarity(Answer &answer, unsigned int i, unsigned int j, unsigned int k, unsigned int t) {
    return params.table[t][i][j] * pow(answer.Relevance[k][t],SSCARD::Q);
  }

  virtual double compute_criterion(Answer &answer) {
    double criterion = 0.0, restriction = 0.0, relevance_bias = 0.0;
    for(unsigned int k = 0; k < params.C; ++k) {
      double numerator = 0;
      for(unsigned int i = 0; i < params.N; ++i) {
        for(unsigned int j = 0; j < params.N; ++j) {
          double value = 0;
          for(unsigned int t = 0; t < params.T; ++t) {
            value += dissimilarity(answer,i,j,k,t);
          }
          numerator += pow(answer.U[i][k] * answer.U[j][k], 2.0) * value;
        }
      }
      double denominator = 0;
      for(unsigned int i = 0; i < params.N; ++i) {
        denominator += pow(answer.U[i][k], 2.0);
      }
      VALIDATE_DENOMINATOR(denominator);
      criterion += numerator / 2.0 / denominator;
    }
    for(auto x : params.must_link) {
      unsigned int l = x.first, m = x.second;
      for(unsigned int r = 0; r < params.C; ++r) {
        for(unsigned int s = 0; s < params.C; ++s) {
          if(r != s) {
            restriction += answer.U[l][r] * answer.U[m][s];
          }
        }
      }
    }
    for(auto x : params.cannot_link) {
      unsigned int l = x.first, m = x.second;
      for(unsigned int r = 0; r < params.C; ++r) {
        restriction += answer.U[l][r] * answer.U[m][r];
      }
    }
    answer.restriction = restriction;
    for(unsigned int k = 0; k < params.C; ++k) {
      for(unsigned int t = 0; t < params.T; ++t) {
        relevance_bias += pow(answer.Relevance[k][t],2.0);
      }
    }
    answer.relevance_bias = relevance_bias;
    return criterion + answer.alpha * restriction + answer.beta * relevance_bias;
  }

  virtual void initialize(Answer &answer, unsigned int init, unsigned int iter) {
    Method::initialize(answer, init, iter);
  }

  virtual void set_default_relevance(Answer &answer) {
    for(unsigned int k = 0; k < params.C; ++k) {
      for(unsigned int t = 0; t < params.T; ++t) {
        answer.Relevance[k][t] = 1.0 / params.T;
      }
    }
  }

  virtual void srand(Answer &answer) {
    // initial fuzzy partition
    for(unsigned int i = 0; i < params.N; ++i) {
      double sum = 0;
      for(unsigned int k = 0; k < params.C; ++k) {
        answer.U[i][k] = random.rand_unsigned();
        sum += answer.U[i][k];
      }
      VALIDATE_DENOMINATOR(sum);
      for(unsigned int k = 0; k < params.C; ++k) {
        answer.U[i][k] /= sum;
      }
    }
    // initial relevance weight
    set_default_relevance(answer);
    // update partition
    update_clusters(answer);
    answer.criterion = compute_criterion(answer);
  }

  virtual bool optimize(Answer &answer) {
    update_membership(answer);
    update_relevance(answer);
    return Method::optimize(answer);
  }

  virtual void update_relevance(Answer &answer) {
    // equation [22]
    // preprocess
    Matrix D(params.C,Row(params.T,0));
    for(unsigned int k = 0; k < params.C; ++k) {
      for(unsigned int t = 0; t < params.T; ++t) {
        for(unsigned int i = 0; i < params.N; ++i) {
          for(unsigned int j = 0; j < params.N; ++j) {
            D[k][t] += pow(answer.U[i][k] * answer.U[j][k], 2.0) * params.table[t][i][j];
          }
        }
      }
    }
    // assignment
    for(unsigned int k = 0; k < params.C; ++k) {
      bool good = true;
      for(unsigned int t = 0; t < params.T; ++t) {
        if(Util::cmp(D[k][t]) <= 0) {
          good = false;
        }
      }
      if(good) {
        for(unsigned int t = 0; t < params.T; ++t) {
          double value = 0;
          for(unsigned int p = 0; p < params.T; ++p) {
            VALIDATE_DENOMINATOR(D[k][p]);
            value += pow(D[k][t] / D[k][p], 1.0 / (SSCARD::Q-1));
          }
          VALIDATE_DENOMINATOR(value);
          answer.Relevance[k][t] = 1.0 / value;
        }
      }
    }
  }

  void update_membership(Answer& answer) {
    Matrix Upow2(params.N,Row(params.C,0));
    for(unsigned int i = 0; i < params.N; ++i) {
      for(unsigned int k = 0; k < params.C; ++k) {
        Upow2[i][k] = pow(answer.U[i][k], 2.0);
      }
    }
    // equation [5]
    Matrix a(params.N,Row(params.C,0));
    // first term
    for(unsigned int i = 0; i < params.N; ++i) {
      for(unsigned int k = 0; k < params.C; ++k) {
        double numerator = 0, denominator = 0;
        for(unsigned int j = 0; j < params.N; ++j) {
          double value = 0;
          for(unsigned int t = 0; t < params.T; ++t) {
            value += dissimilarity(answer,i,j,k,t);
          }
          numerator += Upow2[j][k] * value;
        }
        for(unsigned int j = 0; j < params.N; ++j) {
          denominator += Upow2[j][k];
        }
        VALIDATE_DENOMINATOR(denominator);
        a[i][k] += 2.0 * numerator / denominator;
      }
    }
    // second term
    for(unsigned int k = 0; k < params.C; ++k) {
      double numerator = 0, denominator = 0;
      for(unsigned int i = 0; i < params.N; ++i) {
        for(unsigned int j = 0; j < params.N; ++j) {
          double value = 0;
          for(unsigned int t = 0; t < params.T; ++t) {
            value += dissimilarity(answer,i,j,k,t);
          }
          numerator += Upow2[i][k] * Upow2[j][k] * value;
        }
      }
      for(unsigned int i = 0; i < params.N; ++i) {
        denominator += Upow2[i][k];
      }
      VALIDATE_DENOMINATOR(denominator);
      for(unsigned int i = 0; i < params.N; ++i) {
        a[i][k] -= numerator / denominator / denominator;
      }
    }
    // beta-spread transform
    vector<Row> v(params.N,Row(params.C,0));
    for(unsigned int k = 0; k < params.C; ++k) {
      double denominator = 0;
      for(unsigned int i = 0; i < params.N; ++i) {
        v[i][k] = Upow2[i][k];
        denominator += v[i][k];
      }
      VALIDATE_DENOMINATOR(denominator);
      for(unsigned int i = 0; i < params.N; ++i) {
        v[i][k] /= denominator;
      }
    }
    // fix
    double delta = 0;
    for(unsigned int i = 0; i < params.N; ++i) {
      for(unsigned int k = 0; k < params.C; ++k) {
        if(Util::cmp(a[i][k]) < 0) {
          double normalize = 0;
          for(unsigned int j = 0; j < params.N; ++j) {
            normalize += pow(v[j][k] - (i == j), 2.0);
          }
          VALIDATE_DENOMINATOR(normalize);
          normalize = -a[i][k] / normalize;
          delta = max(delta, normalize);
        }
      }
    }
    if(Util::cmp(delta) > 0) {
      for(unsigned int i = 0; i < params.N; ++i) {
        for(unsigned int k = 0; k < params.C; ++k) {
          double normalize = 0;
          for(unsigned int j = 0; j < params.N; ++j) {
            normalize += pow(v[j][k] - (i == j), 2.0);
          }
          a[i][k] += 2 * normalize * delta;
          ASSERT(Util::cmp(a[i][k]) >= 0, "dist matrix values should be non-negative");
        }
      }
    }
    // u-rfcm
    Matrix u_rfcm(params.N,Row(params.C,0));
    for(unsigned int i = 0; i < params.N; ++i) {
      for(unsigned int k = 0; k < params.C; ++k) {
        if(Util::cmp(a[i][k]) > 0) {
          double denominator = 0;
          for(unsigned int h = 0; h < params.C; ++h) {
            if(Util::cmp(a[i][h]) > 0) {
              denominator += a[i][k] / a[i][h];
            }
          }
          VALIDATE_DENOMINATOR(denominator);
          u_rfcm[i][k] = 1.0 / denominator;
        } else {
          u_rfcm[i][k] = 0.0;
        }
      }
    }
    // equation[5]
    Matrix c(params.N,Row(params.C,0));
    for(unsigned int i = 0; i < params.N; ++i) {
      for(unsigned int j = 0; j < params.N; ++j) {
        // must-link
        if(params.must_link.count(Pair(i,j))) {
          for(unsigned int r = 0; r < params.C; ++r) {
            for(unsigned int s = 0; s < params.C; ++s) {
              if(r != s) {
                c[i][r] += answer.U[j][s];
              }
            }
          }
        }
        if(params.must_link.count(Pair(j,i))) {
          for(unsigned int r = 0; r < params.C; ++r) {
            for(unsigned int s = 0; s < params.C; ++s) {
              if(r != s) {
                c[i][r] += answer.U[j][s];
              }
            }
          }
        }
        // cannot link
        if(params.cannot_link.count(Pair(i,j))) {
          for(unsigned int r = 0; r < params.C; ++r) {
            c[i][r] += answer.U[j][r];
          }
        }
        if(params.cannot_link.count(Pair(j,i))) {
          for(unsigned int r = 0; r < params.C; ++r) {
            c[i][r] += answer.U[j][r];
          }
        }
      }
    }
    // next term
    vector<double> ct(params.N,0.0);
    for(unsigned int i = 0; i < params.N; ++i) {
      double numerator = 0, denominator = 0;
      for(unsigned int k = 0; k < params.C; ++k) {
        if(Util::cmp(a[i][k]) > 0) {
          numerator += c[i][k] / a[i][k];
          denominator += 1.0 / a[i][k];
        }
      }
      VALIDATE_DENOMINATOR(denominator);
      ct[i] = numerator / denominator;
    }
    // u-const
    Matrix u_const(params.N,Row(params.C,0));
    for(unsigned int i = 0; i < params.N; ++i) {
      for(unsigned int k = 0; k < params.C; ++k) {
        if(Util::cmp(a[i][k]) > 0) {
          u_const[i][k] = answer.alpha * (ct[i] - c[i][k]) / a[i][k];
        }
      }
    }
    // assignment
    for(unsigned int i = 0; i < params.N; ++i) {
      for(unsigned int k = 0; k < params.C; ++k) {
        answer.U[i][k] = u_rfcm[i][k] + u_const[i][k];
      }
    }
    // fix and normalize
    if(true) {
      // fix values
      for(unsigned int i = 0; i < params.N; ++i) {
        for(unsigned int k = 0; k < params.C; ++k) {
          if(Util::cmp(answer.U[i][k],0.0) <= 0) {
            answer.U[i][k] = 0.0;
          } else if(Util::cmp(answer.U[i][k],1.0) >= 0) {
            answer.U[i][k] = 1.0;
          }
        }
        // normalize
        double sum = accumulate(answer.U[i].begin(),answer.U[i].end(),(double)(0.0));
        VALIDATE_DENOMINATOR(sum);
        for(unsigned int k = 0; k < params.C; ++k) {
          answer.U[i][k] /= sum;
        }
      }
    }
  }

};

#endif
