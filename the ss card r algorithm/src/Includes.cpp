#include "Includes.h"

int cmp(double x, double y) {
  return fabs(x-y) < eps ? 0 : x < y ? -1 : +1;
}

double sqr(double x) {
  return (x * x);
}
