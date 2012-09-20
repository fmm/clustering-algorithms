// Pattern.h: interface for the Pattern class.
//
//////////////////////////////////////////////////////////////////////

#ifndef PATTERN_H_
#define PATTERN_H_

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

#include "../src/Includes.h"

class Pattern {
private:
	int classePriori; //Classe verdadeira do padrão
	int classe; //Classe atribuida pelo método

	double *grauPertinencia; //Grau de pertinência do padrão em relação às partições

	double *xL; //Vetor de variaveis que descreve o limite inferior do padrão
	double *xU; //Vetor de variaveis que descreve o limite superior do padrão

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
