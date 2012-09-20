// MonteCarlo.h: interface for the MonteCarlo class.
//
//////////////////////////////////////////////////////////////////////

#ifndef MONTECARLO_H_
#define MONTECARLO_H_

#include "../src/Includes.h"
#include "Pattern.h"

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

//Representa uma amostra de uma distribui��o normal bivariada
struct AmostraNormalBi {
	double x;
	double y;
};

////////////////////////////////////

class MonteCarlo {
public:
	MonteCarlo();
	virtual ~MonteCarlo();

	void executar();
	void interfaceMonteCarlo();

private:
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
	Pattern *tabela; //Tabela de dados
	int metodo; //Indica o m�todo escolhido
	int *numeroPattern; //Indica o n�mero de objetos de cada classe

	double **medias; //Vetor de m�dias das vari�veis das parti��es
	double **desvioPadrao; //Vetor de desvio padr�o das vari�veis das parti��es
	double *ro; //Este vetor armazena as covari�ncias das vari�veis para o caso bivariado
	int *limitesGama;
	int b;
	int c;

protected:
	int tipoDeInicializacao;
	int tipoDeGama;
	double geraGamaX(AmostraNormalBi a);
	double geraGamaY(AmostraNormalBi a);
	double geraGama();

	AmostraNormalBi geraNormalBivariada(int i, double **medias, double **desvioPadrao, double *ro);
	void geraDados();

	void geraNormalMultivariada(); //Esta fun��o precisa ser construida
};

#endif
