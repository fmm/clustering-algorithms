#ifndef NORMALIZED_METRICS_H_
#define NORMALIZED_METRICS_H_

#include "lib.h"

namespace Metrics {

  struct NormalizedMetrics {
    virtual double operator()(Row p, Row q) const = 0;
  };

  const struct L1 : public NormalizedMetrics {
    virtual double operator()(Row p, Row q) const {
      double d = 0;
      for(unsigned int i = 0; i < p.size(); ++i) 
        d += fabs(p[i] - q[i]);
      return d / 2.0;
    }
  } L1;
};

#endif
