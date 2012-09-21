#include "Includes.h"

int cmp(double x, double y) {
  return (x <= y + eps) ? (x + eps < y) ? -1 : 0 : 1;
}

double sqr(double x) {
  return (x * x);
}

