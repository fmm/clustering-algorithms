#ifndef IMPRIME_H_
#define IMPRIME_H_

#include <iostream>
#include <iomanip>

using namespace std;

template<class _Tp> class Array;

void printLinha1(ostream&);

void printLinha2(ostream&);

void printInicializacaoHeader(ostream&, int);

void print2(ostream&);

void print3(ostream&, int, double, double);

void printLinha1(ostream&, int);

void printLinha2(ostream&, int);

void printRelacaoHeader(ostream&, int);

void print2(ostream&, int, int, Array<double>&, int, int);

#endif /* IMPRIME_H_ */
