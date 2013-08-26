#ifndef METHOD_H_
#define METHOD_H_

#include "database.h"
#include "mersennetwister.h"
#include "parameter.h"
#include "validation.h"

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
    double eps;
    // compare function
    bool operator<(const Answer& answer) const {
      return Util::cmp(criterion, answer.criterion, eps) < 0;
    }
  };

  bool check_before_update;

  Method(Parameter &params) : params(params) {
    timer.start_timer();
    random = MersenneTwister(params.seed);
    check_before_update = true;
  }

  void save_iteration(Answer &answer, string key) {
    vector< vector<unsigned int> > confusing_matrix = compute_confusion_matrix(answer);
    Matrix priori_matrix = compute_priori_matrix();
    string sql = "INSERT INTO answer("
      "algorithm_id,"
      "initialization,"
      "iteration,"
      "criterion,"
      "restriction,"
      "accuracy,"
      "adjusted_rand_index,"
      "f_measure,"
      "fuzzy_rand_index_campello,"
      "fuzzy_rand_index_hullermeier)"
      "VALUES(" +
      string("\"" + key + "\"") + "," +
      Util::cast<string>(answer.initialization) + "," +
      Util::cast<string>(answer.iteration) + "," +
      Util::cast<string>(answer.criterion) + "," +
      Util::cast<string>(answer.restriction) + "," +
      Util::cast<string>(Validation::accuracy(confusing_matrix).first) + "," +
      Util::cast<string>(Validation::adjusted_rand_index(confusing_matrix)) + "," +
      Util::cast<string>(Validation::f_measure(confusing_matrix)) + "," +
      Util::cast<string>(Validation::fuzzy_rand_index_campello(answer.U,priori_matrix)) + "," +
      Util::cast<string>(Validation::fuzzy_rand_index_hullermeier(answer.U,priori_matrix)) + 
      ");";
    params.database.execute(sql);              
  }

  void save_best(Answer &answer, string key) {
    string sql = "";
    // save partition
    for(unsigned int i = 0; i < params.N; ++i) {
      for(unsigned int k = 0; k < params.C; ++k) {
        sql = "INSERT INTO partition("
          "algorithm_id,"
          "individual,"
          "cluster,"
          "value)"
          "VALUES(" +
          string("\"" + key + "\"") + "," +
          Util::cast<string>(i) + "," +
          Util::cast<string>(k) + "," +
          Util::cast<string>(answer.U[i][k]) +
          ");";
        params.database.execute(sql);
      }
    }
    // save relevance vector
    for(unsigned int k = 0; k < params.C; ++k) {
      for(unsigned int t = 0; t < params.T; ++t) {
        sql = "INSERT INTO relevance("
          "algorithm_id,"
          "cluster,"
          "matrix,"
          "value)"
          "VALUES(" +
          string("\"" + key + "\"") + "," +
          Util::cast<string>(k) + "," +
          Util::cast<string>(t) + "," +
          Util::cast<string>(answer.Relevance[k][t]) +
          ");";
        params.database.execute(sql);
      }
    }
    // update parameters used during the algorithm
    sql = "UPDATE algorithm SET "
      "used_pwc_file=\"" + params.pwc_file + "\"," +
      "used_clusters=" + Util::cast<string>(params.C) + "," +
      "used_prototypes=" + Util::cast<string>(params.P) + "," + 
      "used_alpha=" + Util::cast<string>(params.alpha) + "," +
      "best_initialization=" + Util::cast<string>(answer.initialization) + " " +
      "WHERE sha1=\"" + key + "\";";
    params.database.execute(sql);
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
    if(check_before_update) {
      ASSERT(
          new_criterion <= old_criterion or log(fabs(new_criterion - old_criterion)) <= answer.eps,
          "failed to minimize criterion"
          "(old=" + Util::cast<string>(old_criterion) + 
          ",new=" + Util::cast<string>(new_criterion) + ")"
          );
    }
    answer.criterion = new_criterion;
    if(fabs(new_criterion - old_criterion) > answer.eps) {
      ++answer.iteration;
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

  const Matrix compute_priori_matrix() {
    Matrix Priori(params.N,Row(params.C,0));
    for(unsigned int i = 0; i < params.N; ++i) {
      for(unsigned int k = 0; k < params.C; ++k) {
        Priori[i][k] = params.priori_cluster.size() > k and params.priori_cluster[k].count(i);
      }
    }
    return Priori;
  }

  string status(double percentage) {
    stringstream message;	  
    time_t ending = timer.elapsed() * (1 / percentage - 1);
    message << std::setw(2) << std::setfill('0')
      << Util::cast<int>(percentage*100) << "%, time left = ";
    // days
    if(unsigned long long d = ending / (1000000ULL * 60 * 60 * 24)) {
      message << d << "d:";
      ending %= (1000000ULL * 60 * 60 * 24);
    }
    // hours
    unsigned long long h = ending / (1000000ULL * 60 * 60);
    message << std::setw(2) << std::setfill('0')
      << h << "hs:";
    ending %= (1000000ULL * 60 * 60);
    // minutes
    unsigned long long m = ending / (1000000ULL * 60);
    message << std::setw(2) << std::setfill('0')
      << m << "min:";
    ending %= (1000000ULL * 60);
    // seconds
    unsigned long long s = ending / (1000000ULL);
    message << std::setw(2) << std::setfill('0')
      << s << "s:";
    ending %= (1000000ULL);
    // miliseconds
    message << std::setw(6) << std::setfill('0') << ending;
    return message.str();
  }

  // TODO: considering that every values is already defined
  Answer process() {
    // hash of the execution
    const string key = params.sha1;
    // control the number of status prints
    unsigned long long mask = 1;
    while(mask<= params.initialization) mask<<=1;
    mask=(mask>>4)-1;
    // begin transaction
    params.database.open_transaction();
    params.save();
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
      save_iteration(now,key);
      for(unsigned int iter = 1; iter <= params.maximum_iteration; ++iter) {
        if(optimize(now)) {
          save_iteration(now,key);
        } else {
          break;
        }
      }
      // optimize the best result
      best = min(best, now);
      if(init == 1 or init == params.initialization or !(init&mask)) {
        ALERT(status((double)(init)/params.initialization));
      }
    }
    save_best(best, key);
    // end transaction
    params.database.close_transaction();
    // TODO: print result in a tex/pdf file
    return best;
  }

};

#endif
