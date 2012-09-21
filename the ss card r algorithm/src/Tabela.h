#ifndef TABELA_H_
#define TABELA_H_

////////////////////////////////////
//Variáveis Globais
extern int numeroClassesPriori; //Número de classes a priori
extern int numeroClasses; //Número de classes na partição
extern double m; //Parâmetro de nebulosidade (maior que 1)
extern int limiteIteracao; //Número máximo de iterações
extern double epsilon; //Diferença mínima entre os critérios W
extern int numeroVariaveis; //Número de variáveis
extern int totalPadroes; //Número total de padrões
extern int numeroInicializacoes; //Número de inicializações da partição para cada conjunto
extern int numeroMonteCarlo; //Número de ciclos do Monte Carlo
extern int numeroClassesEscolhidas;
////////////////////////////////////

#include "Includes.h"
#include "Array.h"
#include "../sodas/Pattern.h"

struct Tabela {
  int n;
  Array< Array<double> > matriz;

  Tabela(size_t = 0);
  Tabela(string);
  Tabela(Pattern *, Funcao = ABS);
  virtual ~Tabela();

  double& operator()(int, int);
  const double& operator()(int, int) const;
  operator string() const;

  friend istream& operator>>(istream&, Tabela&);
  friend ostream& operator<<(ostream&, const Tabela&);
};
#endif /* TABELA_H_ */
