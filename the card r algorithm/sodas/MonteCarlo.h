// MonteCarlo.h: interface for the MonteCarlo class.
//
//////////////////////////////////////////////////////////////////////

#ifndef MONTECARLO_H_
#define MONTECARLO_H_

#include "../src/Includes.h"
#include "Pattern.h"

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

//Representa uma amostra de uma distribuição normal bivariada
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
	int metodo; //Indica o método escolhido
	int *numeroPattern; //Indica o número de objetos de cada classe

	double **medias; //Vetor de médias das variáveis das partições
	double **desvioPadrao; //Vetor de desvio padrão das variáveis das partições
	double *ro; //Este vetor armazena as covariâncias das variáveis para o caso bivariado
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

	void geraNormalMultivariada(); //Esta função precisa ser construida
};

#endif
