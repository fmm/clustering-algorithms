// Pattern.h: interface for the Pattern class.
//
//////////////////////////////////////////////////////////////////////

#ifndef PATTERN_H_
#define PATTERN_H_

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

#include "../src/Includes.h"

class Pattern {
private:
	int classePriori; //Classe verdadeira do padr�o
	int classe; //Classe atribuida pelo m�todo

	double *grauPertinencia; //Grau de pertin�ncia do padr�o em rela��o �s parti��es

	double *xL; //Vetor de variaveis que descreve o limite inferior do padr�o
	double *xU; //Vetor de variaveis que descreve o limite superior do padr�o

	char * nome; //variavel que armazenara o nome do padrao

public:

	Pattern();
	virtual ~Pattern();

	int getClassePriori();
	void setClassePriori(int cPriori);
	int getClasse();
	void setClasse(int c);
	double getGrauPertinencia(int i);
	double * getGrauPertinencia();
	void setGrauPertinencia(int i, double valor);

	double * getXL();
	double getXL(int j);
	double * getXU();
	double getXU(int j);

	char * getNome();
	void setNome(char * n);

	friend double dissimilaridade(Pattern& , Pattern& , Funcao = ABS);
};

#endif
