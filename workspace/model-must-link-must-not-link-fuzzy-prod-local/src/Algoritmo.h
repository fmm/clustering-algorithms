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

  const size_t inicializacoes; // numero de inicializacoes
  const size_t clusters; // numero de clusters
  const size_t prototipos; // numero de prototipos
  const Repositorio& repositorio; // dados
  ostream& out; // arquivo de saida

  size_t individuos; // numero de individuos
  size_t limite; // limite de iterações

  Resultado melhor, atual;

  Algoritmo(size_t, size_t, size_t,size_t, size_t, double,const Repositorio&, ostream&);
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
