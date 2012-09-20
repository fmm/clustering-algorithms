// MonteCarlo.cpp: implementation of the MonteCarlo class.
//
//////////////////////////////////////////////////////////////////////

#include "MonteCarlo.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

MonteCarlo::MonteCarlo() {

}

MonteCarlo::~MonteCarlo() {

}

//////////////////////////////////////////////////////////////////////
//Esta função é responsável pela interface
void MonteCarlo::interfaceMonteCarlo() {
	// colocar aqui a interface...
	/*do
	 {
	 cout<<"\nEscolha o metodo desejado:\n";
	 cout<<"0 - Metodo Cobertura\n";

	 cin>>metodo;
	 cout<<"\n";
	 }
	 while(metodo != 0);*/

	// TODO
	///////////////////////////////////////////////////
	//Leitura do nome do arquivo de saída
	//cout << "\nDigite o nome do arquivo de saida: ";
	nomeSaida = new char[60];
	strcpy(nomeSaida, "TESTE");
	//cin >> nomeSaida;


	// TODO
	//Leitura da configuração dos dados a partir de arquivos
	char *nomeEntrada = new char[60];
	//cout << "\nDigite o nome do arquivo que contem a configuracao dos dados\n(vide template_simulacao.txt): ";
	//cout.flush();
	//cin >> nomeEntrada;
	strcpy(nomeEntrada, "fmm_simu.txt");

	ifstream fin;
	fin.open(nomeEntrada);
	assert(fin);

	nomeDados = new char[60];

	//Lê nome dos dados
	fin.getline(nomeDados, 60);
	//Lê número de ciclos monte carlo
	fin >> numeroMonteCarlo;
	//Lê numero de inicializações
	fin >> numeroInicializacoes;
	//Lê limite de iterações
	fin >> limiteIteracao;

	//Lê número de variáveis
	fin >> numeroVariaveis;
	//Lê número de classes
	fin >> numeroClasses;

	//Lê número de padrões
	fin >> totalPadroes;

	//Lê o tipo de inicialização
	fin >> tipoDeInicializacao;

	//Lê o tipo de geração dos gamas
	fin >> tipoDeGama;

	if (tipoDeGama == 1) {
		//Lê os limites de gama
		limitesGama = new int[2];
		fin >> limitesGama[0];
		fin >> limitesGama[1];
	} else {
		fin >> b;
		fin >> c;
	}

	numeroClassesPriori = numeroClasses;

	//Preparação para a leitura dos dados de cada classe
	numeroPattern = new int[numeroClasses];

	medias = new double*[numeroClasses];

	for (int i = 0; i < numeroClasses; i++) {
		medias[i] = new double[numeroVariaveis];
	}

	desvioPadrao = new double*[numeroClasses];

	for (int i = 0; i < numeroClasses; i++) {
		desvioPadrao[i] = new double[numeroVariaveis];
	}

	if (numeroVariaveis == 2) {
		ro = new double[numeroClasses];
	}

	//Leitura dos dados de cada classe
	for (int i = 0; i < numeroClasses; i++) {
		fin >> numeroPattern[i];

		for (int j = 0; j < numeroVariaveis; j++) {
			fin >> medias[i][j];
			fin >> desvioPadrao[i][j];
		}

		if (numeroVariaveis == 2) {
			fin >> ro[i];
		}
	}

	fin.close();

	delete nomeEntrada;
}

//////////////////////////////////////////////////////////////////////
void MonteCarlo::executar() {

	/*

	time_t timetamp;
	time(&timetamp);
	srand(timetamp);

	/////////////////////////////
	//Criação do Arquivo de Saída
	ofstream saida;
	char extensao[5];
	strcpy(extensao, ".wri");
	char *arquivoSimulacao = new char[65];
	strcpy(arquivoSimulacao, nomeSaida);
	strncat(arquivoSimulacao, extensao, 4);
	saida.open(arquivoSimulacao);
	/////////////////////////////

	////////////////////////////
	//Impressão
	saida << "Results: \n\n";
	saida << "\n" << setw(25) << "Input Data: " << nomeDados;
	saida << "\n" << setw(25) << "Parameter m: " << m;
	saida << "\n" << setw(25) << "Clusters: " << numeroClasses;
	saida << "\n" << setw(25) << "Initialization Number: " << numeroInicializacoes;
	saida << "\n" << setw(25) << "Iteration Limit: " << limiteIteracao;
	////////////////////////////

	bestJ = new double[numeroMonteCarlo];
	bestT = new double[numeroMonteCarlo];
	bestR = new double[numeroMonteCarlo];
	bestB = new double[numeroMonteCarlo];
	bestCR = new double[numeroMonteCarlo];
	bestDC = new double[numeroMonteCarlo];
	bestGlobalError = new double[numeroMonteCarlo];

	bestClassError = new double*[numeroClasses];
	for (int i = 0; i < numeroClasses; i++) {
		bestClassError[i] = new double[numeroMonteCarlo];
	}

	iteracoesMedio = new double[numeroMonteCarlo];

	//saida << "\n\n --- x --- \n\n";
	// geraDados();

	// Tabela tab(tabela);

	// saida << tab;

	vector<Tabela> v;
	for (int mc = 0; mc < numeroMonteCarlo; mc++) {
		//print("gerando dados...");
		geraDados();
		//print("... dados gerados.\ncriando repositorio...");
		v.push_back(Tabela(tabela));
		//print("funfou!!");
		//print("criando repositorio...");
		Repositorio repositorio(v);
		//print("repositorio criado...");
		//print("inicializando algoritmo...");
		Algoritmo algoritmo(numeroInicializacoes, numeroClasses, 1, repositorio, cout);
		//print("preparando para executar...");
		//print("executando...");
		algoritmo.executar();
		//print("terminado...");


		bestJ[mc] = algoritmo.melhor.J;
		bestT[mc] = 0.0;
		bestR[mc] = 0.0;
		bestB[mc] = 0.0;
		bestCR[mc] = algoritmo.melhor.CR;
		bestDC[mc] = 0.0;
		bestGlobalError[mc] = 0.0;

		for (int i = 0; i < numeroClasses; i++) {
			bestClassError[i][mc] = 0.0;
		}

		for (int j = 0; j < totalPadroes; j++) {
			tabela[j].~Pattern();
		}

		cout << "\n\nReplication: " << mc + 1 << "\nJ = " << bestJ[mc] << "\nCR = " << bestCR[mc];
		cout.flush();
	}

	imprimeResultado(saida);

	*/

	/*for (int i = 0; i < totalPadroes; i++) {
	 saida << "individuo[" << i << "] == ";
	 saida << "(" << tabela[i].getXL(0) << " : " << tabela[i].getXU(0) << ")";
	 saida << ", (" << tabela[i].getXL(1) << " : " << tabela[i].getXU(1) << ")\n";
	 }*/

	//bool rotulado = true;
	//switch (metodo) {

	/*case 11:
	 {
	 ///////////////////
	 //Impressão
	 saida<<"\n"<<setw(25)<<"Algorithm: "<<"Weighted-Single Hard C-Means\n";
	 ///////////////////

	 Metodo11 * metodo11;

	 for(int mc=0;mc<numeroMonteCarlo;mc++)
	 {
	 geraDados();

	 metodo11 = new Metodo11();
	 metodo11 -> gerarObjetosPadroes(tabela, true);

	 metodo11 -> executarSimulacao(nomeDados,nomeSaida,tipoDeInicializacao);

	 bestJ[mc] = metodo11 -> getMelhorJ();
	 bestT[mc] = 0.0;
	 bestR[mc] = 0.0;
	 bestB[mc] = 0.0;
	 bestCR[mc] = metodo11 -> getMelhorCR();
	 bestDC[mc] = 0.0;
	 bestGlobalError[mc] = metodo11 -> calcularErrosClassificacao(nomeSaida);

	 for(int i=0;i<numeroClasses;i++)
	 {
	 bestClassError[i][mc] = metodo11 -> getMelhorErroClasse(i);
	 }

	 iteracoesMedio[mc] = metodo11 -> getMediaIteracoes();

	 for(int j = 0; j < totalPadroes; j++)
	 {
	 tabela[j].~Pattern();
	 }

	 cout<<"\n\nReplication: "<<mc+1<<"\nJ = "<<bestJ[mc]<<"\nCR = "<<bestCR[mc];
	 cout.flush();

	 }
	 imprimeResultado(saida);
	 }
	 break;*/

	//}

	////////////////////////////////////////
	// TESTE
	/*_dbg(saida,m);
	 _dbg(saida,numeroParticao);
	 _dbg(saida,epsilon);*/
	/*
	 extern int numeroVariaveis;
	 extern int limiteIteracao;
	 extern int totalPadroes;
	 extern int numeroInicializacoes;
	 extern int numeroMonteCarlo;
	 extern int numeroClassesPriori;*/
	////////////////////////////////////////

	//imprimeResultado(saida);

	/*delete bestJ;
	delete bestT;
	delete bestR;
	delete bestB;
	delete bestCR;
	delete bestDC;
	delete bestGlobalError;

	for (int i = 0; i < numeroClasses; i++) {
		delete bestClassError[i];
	}
	delete bestClassError;

	delete iteracoesMedio;*/
}

//////////////////////////////////////////////////////////////////////
//Esta função gera os dados de entrada
void MonteCarlo::geraDados() {
	tabela = new Pattern[totalPadroes];

	//Se os dados forem bivariados
	if (numeroVariaveis == 2) {
		if (tipoDeGama == 1) {
			for (int i = 0, j = 0; i < numeroClasses; i++) {
				for (int k = 0; k < numeroPattern[i]; k++) {
					double *xL = tabela[j].getXL();
					double *xU = tabela[j].getXU();

					AmostraNormalBi a = geraNormalBivariada(i, medias, desvioPadrao, ro);

					double gama1 = geraGama();
					double gama2 = geraGama();

					xL[0] = a.x - (gama1 / 2.0);
					xU[0] = a.x + (gama1 / 2.0);//(a*xL)+b

					xL[1] = a.y - (gama2 / 2.0);
					xU[1] = a.y + (gama2 / 2.0);

					tabela[j].setClassePriori(i);

					j++;

				}
			}
		} else {
			for (int i = 0, j = 0; i < numeroClasses; i++) {
				for (int k = 0; k < numeroPattern[i]; k++) {
					double *xL = tabela[j].getXL();
					double *xU = tabela[j].getXU();

					AmostraNormalBi a = geraNormalBivariada(i, medias, desvioPadrao, ro);

					double gama1 = geraGamaX(a);
					double gama2 = geraGamaY(a);

					xL[0] = a.x - (gama1 / 2.0);
					xU[0] = a.x + (gama1 / 2.0);

					xL[1] = a.y - (gama2 / 2.0);
					xU[1] = a.y + (gama2 / 2.0);

					tabela[j].setClassePriori(i);

					j++;
				}
			}
		}
	} else { // TODO Se os dados não forem bivariados, precisa fazer
		return;
	}
}

//////////////////////////////////////////////////////////////////////
void MonteCarlo::geraNormalMultivariada() {

}

//////////////////////////////////////////////////////////////////////
//Esta função gera uma amostra uma normal bivariada
AmostraNormalBi MonteCarlo::geraNormalBivariada(int i, double **medias, double **desvioPadrao, double *ro) {
	//////////////////////////////
	//Algoritmo N2
	double u1 = (rand() % 100) / 100.0;

	double u2 = 0;
	do {
		u2 = (rand() % 100) / 100.0;
	} while (u2 == 0.0);

	double e = -log(u2);

	double z1 = (sqrt(2 * e)) * cos((2 * pi * u1));
	double z2 = (sqrt(2 * e)) * sin((2 * pi * u1));

	//////////////////////////////

	double y = medias[i][1] + desvioPadrao[i][1] * z2;

	double mediaX_Condicional = medias[i][0] + (ro[i] * desvioPadrao[i][0] / desvioPadrao[i][1]) * (y - medias[i][1]);

	double desvioPadraoX_Condicional = sqrt(pow(desvioPadrao[i][0], 2) * (1 - pow(ro[i], 2)));

	double xCondicional = mediaX_Condicional + z1 * desvioPadraoX_Condicional;

	AmostraNormalBi a;

	a.x = xCondicional;
	a.y = y;

	return a;

}

//////////////////////////////////////////////////////////////////////
void MonteCarlo::imprimeResultado(ostream& saida) {
	for (int mc = 0; mc < numeroMonteCarlo; mc++) {
		saida << "\n\nReplication " << mc + 1;
		saida << "\n" << setw(15) << "J" << setw(15) << "T" << setw(15) << "R" << setw(15) << "B" << setw(15) << "CR" << setw(15) << "DC" << setw(15) << "Global Error";

		for (int i = 0; i < numeroClasses; i++) {
			saida << setw(7) << "Class" << i + 1 << setw(7) << "Error";
		}

		saida << setw(25) << "Medium Iteration Number";

		saida << "\n----------------------------------------------------------------------------------------------------------------------------------";
		saida << setw(15 * numeroClasses + 1) << setfill('-') << " " << setfill(' ') << "\n";
		saida << setw(15) << bestJ[mc] << setw(15) << bestT[mc] << setw(15) << bestR[mc] << setw(15) << bestB[mc] << setw(15) << bestCR[mc] << setw(15) << bestDC[mc] << setw(15) << 100 * bestGlobalError[mc];

		for (int i = 0; i < numeroClasses; i++) {
			saida << setw(15) << 100 * bestClassError[i][mc];
		}

		saida << setw(25) << iteracoesMedio[mc];
	}

	int indice = 0;

	for (int mc = 1; mc < numeroMonteCarlo; mc++) {
		if (bestJ[mc] < bestJ[indice]) {
			indice = mc;
		}
	}

	saida << "\n\n";

	double media;

	saida << "\n\n--------------------------------------------\n";
	saida << "\nInformation about J\n";
	media = calculaMedia(bestJ, numeroMonteCarlo);
	saida << setw(20) << "Avg" << setw(20) << "Std";
	saida << "\n";
	saida << setw(20) << media << setw(20) << calculaDesvioPadrao(bestJ, media, numeroMonteCarlo);

	saida << "\n\n--------------------------------------------\n";
	saida << "\nInformation about Adjusted Rand Index (RC)\n";
	media = calculaMedia(bestCR, numeroMonteCarlo);
	saida << setw(20) << "Avg" << setw(20) << "Std";
	saida << "\n";
	saida << setw(20) << media << setw(20) << calculaDesvioPadrao(bestCR, media, numeroMonteCarlo);

	saida << "\n\n--------------------------------------------\n";
	saida << "\nInformation about DC\n";
	media = calculaMedia(bestDC, numeroMonteCarlo);
	saida << setw(20) << "Avg" << setw(20) << "Std";
	saida << "\n";
	saida << setw(20) << media << setw(20) << calculaDesvioPadrao(bestDC, media, numeroMonteCarlo);

	saida << "\n\n--------------------------------------------\n";
	saida << "\nInformation about T\n";
	media = calculaMedia(bestT, numeroMonteCarlo);
	saida << setw(20) << "Avg" << setw(20) << "Std";
	saida << "\n";
	saida << setw(20) << media << setw(20) << calculaDesvioPadrao(bestT, media, numeroMonteCarlo);

	saida << "\n\n--------------------------------------------\n";
	saida << "\nInformation about R\n";
	media = calculaMedia(bestR, numeroMonteCarlo);
	saida << setw(20) << "Avg" << setw(20) << "Std";
	saida << "\n";
	saida << setw(20) << media << setw(20) << calculaDesvioPadrao(bestR, media, numeroMonteCarlo);

	saida << "\n\n--------------------------------------------\n";
	saida << "\nInformation about B\n";
	media = calculaMedia(bestB, numeroMonteCarlo);
	saida << setw(20) << "Avg" << setw(20) << "Std";
	saida << "\n";
	saida << setw(20) << media << setw(20) << calculaDesvioPadrao(bestB, media, numeroMonteCarlo);

	saida << "\n\n--------------------------------------------\n";
	saida << "\nInformation about the Classification Error\n";
	saida << setw(7) << "CLUSTER" << setw(20) << "Avg" << setw(20) << "Std";
	saida << "\n";
	for (int i = 0; i < numeroClasses; i++) {
		media = calculaMedia(bestClassError[i], numeroMonteCarlo);
		saida << setw(7) << i + 1 << setw(20) << media * 100 << setw(20) << 100 * calculaDesvioPadrao(bestClassError[i], media, numeroMonteCarlo);
		saida << "\n";
	}
	media = calculaMedia(bestGlobalError, numeroMonteCarlo);
	saida << setw(7) << "Global" << setw(20) << media * 100 << setw(20) << 100 * calculaDesvioPadrao(bestGlobalError, media, numeroMonteCarlo);

	saida << "\n\n--------------------------------------------\n";
	saida << "\nInformation about the Medium Iteration Number\n";
	media = calculaMedia(iteracoesMedio, numeroMonteCarlo);
	saida << setw(20) << "Avg" << setw(20) << "Std";
	saida << "\n";
	saida << setw(20) << media << setw(20) << calculaDesvioPadrao(iteracoesMedio, media, numeroMonteCarlo);

}

//////////////////////////////////////////////////////////////////////
double MonteCarlo::calculaMedia(double *vetor, int tamanho) {
	double soma = 0.0;
	for (int i = 0; i < tamanho; i++) {
		soma += vetor[i];
	}

	double media = soma / tamanho;

	return media;
}

//////////////////////////////////////////////////////////////////////
double MonteCarlo::calculaDesvioPadrao(double *vetor, double media, int tamanho) {
	double soma = 0.0;
	for (int i = 0; i < tamanho; i++) {
		soma += pow((vetor[i] - media), 2);
	}
	double desvioPadrao = sqrt(soma / (tamanho - 1));

	return desvioPadrao;
}
//////////////////////////////////////////////////////////////////////
double MonteCarlo::geraGama() {
	double semente = (rand() % 100) / 100.0;

	double gama = limitesGama[0] + (limitesGama[1] - limitesGama[0]) * semente;

	return gama;
}

//////////////////////////////////////////////////////////////////////
double MonteCarlo::geraGamaX(AmostraNormalBi a) {
	double gama = ((b * abs(a.x)) + c);

	return gama;
}

//////////////////////////////////////////////////////////////////////
double MonteCarlo::geraGamaY(AmostraNormalBi a) {
	double gama = ((b * abs(a.y)) + c);

	return gama;
}

//////////////////////////////////////////////////////////////////////
