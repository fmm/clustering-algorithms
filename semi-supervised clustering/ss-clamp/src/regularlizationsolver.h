#ifndef REGULARLIZATION_SOLVER_H_
#define REGULARLIZATION_SOLVER_H_

#include "lib.h"

struct RegularlizationSolver {

  vector<double> term;
  double beta;

  RegularlizationSolver(vector<double> term, double beta) {
    this->term = term;
    this->beta = beta;
  }

  vector<double> partial(double gamma) {
    vector<double> vet;
    for(double a : term) {
      double delta = a*a + 4*(2*beta)*(gamma);
      vet.push_back((sqrt(delta)-a)/2/(2*beta));
    }
    return vet;
  }

  double eval(double gamma) {
    vector<double> fgamma = partial(gamma);
    double product = 1.0;
    for(double value : fgamma) {
      product *= value;
    }
    return product;
  }

  bool itsok(double gamma) {
    return Util::cmp(eval(gamma),1.0) <= 0;
  }

  vector<double> process() {
    double step = 1<<30, gamma = 0;
    // this precision should be good enough
    for(unsigned int i = 0; i < 100; i++) {
      double new_gamma = gamma + step;
      if(itsok(new_gamma)) {
        gamma = new_gamma;
      }
      step /= 2;
    }
    return partial(gamma);
  }
};

#endif
