#ifndef ALGORITMO_H_
#define ALGORITMO_H_

#include "Includes.h"
#include "Repositorio.h"
#include "Resultado.h"
#include "MinCostFlow.h"
#include "Imprime.h"

struct Algoritmo {

  static MinCostFlow mcf;

  static string saida;

  const size_t inicializacoes;
  const size_t clusters;
  const size_t prototipos;
  const Repositorio& repositorio;
  ostream& out;

  size_t individuos;
  size_t limite;

  Resultado melhor, atual;

  Algoritmo(size_t, size_t, size_t, size_t, double, const Repositorio&, ostream&);
  virtual ~Algoritmo();

  void executar();
  bool inicializacao();
  void etapa1();
  bool etapa2();

  operator string() const;

  void printAcessFile(size_t []);
  void imprimirMatriz(ostream &);

  friend ostream& operator<<(ostream&, Algoritmo&);
};

#endif /* ALGORITMO_H_ */
