/***************************************/
int numeroClassesPriori;
int numeroClasses;
double m = 2;
int limiteIteracao = 350;
double epsilon = 1e-10;
int numeroVariaveis;
int totalPadroes;
int numeroInicializacoes = 60;
int numeroMonteCarlo = 60;
int numeroClassesEscolhidas = 0;
/***************************************/

#include "src/Algoritmo.h"

vector<string> v;
string arq, saida;
int numCluster, numPrototipo;
int numInicializacao;
int var_classe;
const int parametro_m = 2;
int numIteracoes;

double label;
double alpha;

void read(const char *nome) {
  ifstream in(nome, ios::in);
  string s;
  bool input = false, output = false;
  while (in >> s) {
    if (s.find("(numCluster)") != string::npos) {
      in >> s;
      sscanf(s.c_str(), "%d", &numCluster);
      input = output = false;
    }
    else if (s.find("(numPrototipos)") != string::npos) {
      in >> s;
      sscanf(s.c_str(), "%d", &numPrototipo);
      input = output = false;
    }
    else if (s.find("(numInicializacao)") != string::npos) {
      in >> s;
      sscanf(s.c_str(), "%d", &numInicializacao);
      input = output = false;
    }
    else if (s.find("(variavel_classe)") != string::npos) {
      in >> s;
      sscanf(s.c_str(), "%d", &var_classe);
      Dados::var_classe = var_classe;
      input = output = false;
    }
    else if (s.find("(label)") != string::npos) {
      in >> s;
      sscanf(s.c_str(), "%lf", &label);
      input = output = false;
    }
    else if (s.find("(alfa)") != string::npos) {
      in >> s;
      sscanf(s.c_str(), "%lf", &alpha);
      input = output = false;
    }
    else if (s.find("(numIteracoes)") != string::npos) {
      in >> s;
      sscanf(s.c_str(), "%d", &numIteracoes);
      input = output = false;
    }
    else if (s.find("(input)") != string::npos) {
      input = true;
      output = false;
    } 
    else if (s.find("(output)") != string::npos) {
      input = false;
      output = true;
    }
    else if (input) {
      v.push_back(s);
    }
    else if (output) {
      saida = s;
      Algoritmo::saida = s;
    }
    else {
      assert(0);
    }
  }
  Dbg(numCluster);
  Dbg(numInicializacao);
  Dbg(numIteracoes);
  Dbg(var_classe);
  Print("input:");
  rep(i,v.size()) Print((i+1) << " : " << v[i]);
  Print("output:");
  Print(saida);
}

void wait(int seconds) {
  clock_t endwait;
  endwait = clock() + seconds * CLOCKS_PER_SEC;
  while (clock() < endwait);
}

int main(int argc, char **argv) {

  if(argc > 1) {
    arq = string(argv[1]);
  }
  else {
    cout << "Digite o nome do arquivo de configuracao: " << endl;
    cin >> arq;
  }

  read(arq.c_str());

  srand(7);
  //srand(time(NULL));

  try {
    Repositorio repositorio(v, label);
    ofstream out(saida.c_str(), ios::out);
    Algoritmo algoritmo(numInicializacao, numCluster, numPrototipo, parametro_m, numIteracoes, alpha, repositorio, out);
    algoritmo.executar();
  } catch (exception &error) {
    cout << error.what() << endl;
  } catch (...) {
    cout << "eita esse erro foi muito do mau!!\n";
  }

  cout.flush();
  cerr << ("ending...") << endl;

  wait(1);

  return 0;
}
