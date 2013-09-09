#ifndef SS_CLAMP_H_
#define SS_CLAMP_H_

#include "method.h"

struct SSClamp : public Method {

  // constructor
  SSClamp(Parameter &params) : Method(params) {
  }
  
  virtual inline vector<double> classify(Answer &answer, Matrix &table) {
    ASSERT(table.size() == params.T, "unexpected number of tables");
    ASSERT(table[0].size() == params.N, "unexpected number of individuals");
    vector<double> dist(params.C,0.0), u(params.C,0.0);
    for(unsigned int t = 0; t < params.T; ++t) {
      for(unsigned int k = 0; k < params.C; ++k) {
        for(unsigned int p = 0; p < params.P; ++p) {
          dist[k] += table[t][answer.prototype[k][p]] * relevance(answer, k, t);
        }
      }
    }
    vector<int> is_zero;
    for(unsigned int k = 0; k < params.C; ++k) {
      if(Util::cmp(dist[k]) <= 0) {
        is_zero.push_back(k);
      }
    }
    if(is_zero.size()) {
      for(unsigned int k = 0; k < is_zero.size(); ++k) {
        u[is_zero[k]] = 1.0 / is_zero.size();
      }
    } else {
      for(unsigned int k = 0; k < params.C; ++k) {
        for(unsigned int h = 0; h < params.C; ++h) {
          VALIDATE_DENOMINATOR(dist[h]);
          u[k] += dist[k] / dist[h];
        }
        VALIDATE_DENOMINATOR(u[k]);
        u[k] = 1.0 / u[k];
      }
    }
    return u;
  }

  virtual inline double relevance(Answer &answer, unsigned int k, unsigned int t) {
    return answer.Relevance[k][t];
  }

  inline double dissimilarity(Answer &answer, unsigned int i, unsigned int j, unsigned int k, unsigned int t) {
    ASSERT(params.mask[i] and params.mask[j], "invalid individuals");
    return params.table[t][i][j] * relevance(answer, k, t);
  }

  virtual double compute_criterion(Answer &answer) {
    double criterion = 0.0, restriction = 0.0;
    Matrix d(params.N,Row(params.C,0));
    for(unsigned int t = 0; t < params.T; ++t) {
      for(unsigned int i = 0; i < params.N; ++i) {
        if(params.mask[i]) {
          for(unsigned int k = 0; k < params.C; ++k) {
            for(unsigned int p = 0; p < params.P; ++p) {
              d[i][k] += dissimilarity(answer, i, answer.prototype[k][p], k, t);
            }
          }
        }
      }
    }
    for(unsigned int i = 0; i < params.N; ++i) {
      if(params.mask[i]) {
        for(unsigned int k = 0; k < params.C; ++k) {
          criterion += pow(answer.U[i][k], 2.0) * d[i][k];
        }
      }
    }
    for(auto x : params.must_link) {
      unsigned int l = x.first, m = x.second;
      if(params.mask[l] and params.mask[m]) {
        for(unsigned int r = 0; r < params.C; ++r) {
          for(unsigned int s = 0; s < params.C; ++s) {
            if(r != s) {
              restriction += answer.U[l][r] * answer.U[m][s];
            }
          }
        }
      }
    }
    for(auto x : params.cannot_link) {
      unsigned int l = x.first, m = x.second;
      if(params.mask[l] and params.mask[m]) {
        for(unsigned int r = 0; r < params.C; ++r) {
          restriction += answer.U[l][r] * answer.U[m][r];
        }
      }
    }
    answer.restriction = restriction;
    return criterion + params.alpha * restriction;
  }

  virtual void initialize(Answer &answer, unsigned int init, unsigned int iter) {
    Method::initialize(answer, init, iter);
    answer.prototype = vector<Prototype>(params.C);
  }
  
  virtual void set_default_relevance(Answer &answer) {
    for(unsigned int k = 0; k < params.C; ++k) {
      for(unsigned int t = 0; t < params.T; ++t) {
        answer.Relevance[k][t] = 1.0;
      }
    }
  }

  virtual void srand(Answer &answer) {
    // random prototype
    vector<unsigned int> new_prototype;
    for(unsigned int i = 0; i < params.N; ++i) {
      if(params.mask[i]) {
        new_prototype.push_back(i);
      }
    }
    ASSERT(
        params.C * params.P <= new_prototype.size(),
        "invalid number of prototypes"
        );
    // randomize vector
    for (unsigned int i = 1; i < new_prototype.size(); ++i) {
      unsigned int j = random.rand_unsigned() % (i + 1);
      swap(new_prototype[i],new_prototype[j]);
    }
    for(unsigned int k = 0; k < params.C; ++k) {
      answer.prototype[k].clear();
      for(unsigned int p = 0; p < params.P; ++p) {
        answer.prototype[k].push_back(new_prototype[k * params.P + p]);
      }
    }
    // initial relevance weight
    set_default_relevance(answer);
    // initial fuzzy partition
    for(unsigned int i = 0; i < params.N; ++i) {
      if(params.mask[i]) {
        vector<double> dist(params.C,0.0);
        for(unsigned int t = 0; t < params.T; ++t) {
          for(unsigned int k = 0; k < params.C; ++k) {
            for(unsigned int p = 0; p < params.P; ++p) {
              dist[k] += dissimilarity(answer, i, answer.prototype[k][p], k, t);
            }
          }
        }
        vector<int> is_zero;
        for(unsigned int k = 0; k < params.C; ++k) {
          if(Util::cmp(dist[k]) <= 0) {
            is_zero.push_back(k);
          }
        }
        if(is_zero.size()) {
          for(unsigned int k = 0; k < params.C; ++k) {
            answer.U[i][k] = 0.0;
          }
          for(unsigned int k = 0; k < is_zero.size(); ++k) {
            answer.U[i][is_zero[k]] = 1.0 / is_zero.size();
          }
        } else {
          for(unsigned int k = 0; k < params.C; ++k) {
            answer.U[i][k] = 0;
            for(unsigned int h = 0; h < params.C; ++h) {
              VALIDATE_DENOMINATOR(dist[h]);
              answer.U[i][k] += dist[k] / dist[h];
            }
            VALIDATE_DENOMINATOR(answer.U[i][k]);
            answer.U[i][k] = 1.0 / answer.U[i][k];
          }
        }
      }
    }
    update_clusters(answer);
    answer.criterion = compute_criterion(answer);
  }

  virtual bool optimize(Answer &answer) {
    update_prototypes(answer);
    update_relevance(answer);
    update_membership(answer);
    return Method::optimize(answer);
  }

  void update_prototypes(Answer &answer) {
    for(unsigned int k = 0; k < params.C; ++k) {
      vector<double> dist(params.N,0.0);
      for(unsigned int t = 0; t < params.T; ++t) {
        for(unsigned int i = 0; i < params.N; ++i) {
          if(params.mask[i]) {
            for(unsigned int j = 0; j < params.N; ++j) {
              if(params.mask[j]) {
                dist[i] += dissimilarity(answer, i, j, k, t) * pow(answer.U[j][k],2.0);
              }
            }
          }
        }
      }
      vector< pair<double,int> > v;
      for(unsigned int i = 0; i < params.N; ++i) {
        if(params.mask[i]) {
          v.push_back(make_pair(dist[i], i));
        }
      }
      nth_element(v.begin(),v.begin()+params.P,v.end());
      for(unsigned int p = 0; p < params.P; ++p) {
        answer.prototype[k][p] = v[p].second;
      }
    }
  }

  virtual void update_relevance(Answer &answer) = 0;

  void update_membership(Answer& answer) {
    for(unsigned int i = 0; i < params.N; ++i) {
      if(params.mask[i]) {
        vector<double> a(params.C,0), b(params.C,0);
        // compute a[]
        for(unsigned int k = 0; k < params.C; ++k) {
          for(unsigned int t = 0; t < params.T; ++t) {
            for(unsigned int p = 0; p < params.P; ++p) {
              a[k] += dissimilarity(answer, i, answer.prototype[k][p], k, t);
            }
          }
        }
        // compute b[]
        for(unsigned int m = 0; m < params.N; ++m) {
          if(params.mask[m]) {
            // must link
            if(params.must_link.count(Pair(i,m))) {
              for(unsigned int r = 0; r < params.C; ++r) {
                for(unsigned int s = 0; s < params.C; ++s) {
                  if(r != s) {
                    b[r] += answer.U[m][s];
                  }
                }
              }
            }
            // cannot link
            if(params.cannot_link.count(Pair(i,m))) {
              for(unsigned int r = 0; r < params.C; ++r) {
                b[r] += answer.U[m][r];
              }
            }
          }
        }
        for(unsigned int l = 0; l < params.N; ++l) {
          if(params.mask[l]) {
            // must link
            if(params.must_link.count(Pair(l,i))) {
              for(unsigned int r = 0; r < params.C; ++r) {
                for(unsigned int s = 0; s < params.C; ++s) {
                  if(r != s) {
                    b[s] += answer.U[l][r];
                  }
                }
              }
            }
            // cannot link
            if(params.cannot_link.count(Pair(l,i))) {
              for(unsigned int s = 0; s < params.C; ++s) {
                b[s] += answer.U[l][s];
              }
            }
          }
        }
        for(unsigned int k = 0; k < params.C; ++k) {
          a[k] *= 2;
          b[k] *= params.alpha;
        }
        vector<unsigned int> V;
        for(unsigned int k = 0; k < params.C; ++k) {
          if(Util::cmp(a[k]) <= 0) {
            V.push_back(k);
          }
        }
        if(V.size()) {
          // maximize membership for a[k]=0 in order to do this gamma must be
          // the minimum b[k] for all k with a[k]=0 
          double gamma = INF;
          for(unsigned int k = 0; k < V.size(); ++k) {
            gamma = min(gamma, b[V[k]]);
          }
          double membership = 1.0;
          V.clear();
          for(unsigned int k = 0; k < params.C; ++k) {
            if(Util::cmp(gamma,b[k]) >= 0) {
              if(Util::cmp(a[k]) > 0) {
                answer.U[i][k] = (gamma - b[k]) / a[k];
                membership -= answer.U[i][k];
              } else {
                V.push_back(k);
              }
            } else {
              answer.U[i][k] = 0.0;
            }
          }
          if(Util::cmp(membership) >= 0) {
            // the sum should be equal to `membership'
            // here this value is being equally distributed
            for(unsigned int k = 0; k < V.size(); ++k) {
              answer.U[i][V[k]] = membership / V.size();
            }
            V.clear();
          } else {
            // the b[k] greater equal than gamma should have membership
            // value equal to zero according to the Kuhn-Tucker conditions
            V = vector<unsigned int>(params.C,true);
            for(unsigned int k = 0; k < params.C; ++k) {
              if(Util::cmp(b[k],gamma) >= 0) {
                answer.U[i][k] = 0.0;
                V[k] = false;
              }
            }
          }
        } else {
          V = vector<unsigned int>(params.C,true);
        }
        if(V.size()) {
          while(true) {
            bool updated = false;
            double gamma = ({
              double num = 0, den = 0;
              for(unsigned int k = 0; k < params.C; ++k) {
                if(V[k]) {
                  VALIDATE_DENOMINATOR(a[k]);
                  num += (b[k] / a[k]);
                  den += (1.0 / a[k]);
                }
              }
              VALIDATE_DENOMINATOR(den);
              (1.0 + num) / den;
            });
            for(unsigned int k = 0; k < params.C; ++k) {
              if(V[k]) {
                VALIDATE_DENOMINATOR(a[k]);
                answer.U[i][k] = (gamma - b[k]) / a[k];
                if(Util::cmp(answer.U[i][k]) <= 0) {
                  V[k] = false;
                  updated = true;
                }
              } else {
                answer.U[i][k] = 0;
              }
            }
            if(!updated) {
              break;
            }
          }
        }
      }
      // just to have a better precision
      if(true) {
        // check
        double sum = accumulate(answer.U[i].begin(),answer.U[i].end(),(double)(0.0));
        ASSERT(log(fabs(sum-1.0)) <= answer.eps, "membership values were invalid");
        // extra-normalize
        for(unsigned int k = 0; k < params.C; ++k) {
          if(Util::cmp(answer.U[i][k],0.0) <= 0) {
            answer.U[i][k] = 0.0;
          } else if(Util::cmp(answer.U[i][k],1.0) >= 0) {
            answer.U[i][k] = 1.0;
          }
        }
      }
    }
  }

};

#endif
