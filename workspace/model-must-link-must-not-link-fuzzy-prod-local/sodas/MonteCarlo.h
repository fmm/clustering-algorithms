#ifndef MONTECARLO_H_
#define MONTECARLO_H_

#include "../src/Includes.h"
#include "Pattern.h"

/***************************************/
extern int numeroClassesPriori;
extern int numeroClasses;
extern double m;
extern int limiteIteracao;
extern double epsilon;
extern int numeroVariaveis;
extern int totalPadroes;
extern int numeroInicializacoes;
extern int numeroMonteCarlo;
extern int numeroClassesEscolhidas;
/***************************************/

struct AmostraNormalBi {
	double x;
	double y;
};

struct MonteCarlo {
	MonteCarlo();
	virtual ~MonteCarlo();

	void executar();
	void interfaceMonteCarlo();
	
	double **bestClassError;
	char *nomeDados;
	double *iteracoesMedio;
	double calculaDesvioPadrao(double *vetor, double media, int tamanho);
	double calculaMedia(double *vetor, int tamanho);
	double *bestGlobalError;
	double *bestDC;
	double *bestCR;
	double *bestB;
	double *bestT;
	double *bestR;
	double *bestJ;
	char *nomeSaida;
	void imprimeResultado(ostream& saida);
	Pattern *tabela;
	int metodo;
	int *numeroPattern;
	double **medias;
	double **desvioPadrao;
	double *ro;
	int *limitesGama;
	int b;
	int c;
	int tipoDeInicializacao;
	int tipoDeGama;
	double geraGamaX(AmostraNormalBi a);
	double geraGamaY(AmostraNormalBi a);
	double geraGama();

	AmostraNormalBi geraNormalBivariada(int i, double **medias, double **desvioPadrao, double *ro);
	void geraDados();
	void geraNormalMultivariada();
};

#endif /* MONTECARLO_H_ */
