#ifndef TABELA_H_
#define TABELA_H_

////////////////////////////////////
//Vari�veis Globais
extern int numeroClassesPriori; //N�mero de classes a priori
extern int numeroClasses; //N�mero de classes na parti��o
extern double m; //Par�metro de nebulosidade (maior que 1)
extern int limiteIteracao; //N�mero m�ximo de itera��es
extern double epsilon; //Diferen�a m�nima entre os crit�rios W
extern int numeroVariaveis; //N�mero de vari�veis
extern int totalPadroes; //N�mero total de padr�es
extern int numeroInicializacoes; //N�mero de inicializa��es da parti��o para cada conjunto
extern int numeroMonteCarlo; //N�mero de ciclos do Monte Carlo
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
