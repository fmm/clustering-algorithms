#include "Tabela.h"

Tabela::Tabela(size_t n) :
  n(n) {
    this->matriz.resize(n);
    for (size_t i = 0; i < n; i++) {
      this->matriz[i].resize(n);
    }
  }

Tabela::Tabela(string arquivo) {
  ifstream in(arquivo.c_str(), ios::in);
  in >> *this;
  in.close();
}

Tabela::Tabela(Pattern *pattern, Funcao funcao) {
  *this = Tabela(totalPadroes);
  for (int i = 0; i < n; i++) {
    this->matriz[i][i] = 0.0;
    for (int j = i + 1; j < n; j++) {
      this->matriz[i][j] = this->matriz[j][i] = dissimilaridade(pattern[i], pattern[j], funcao);
    }
  }
}

Tabela::~Tabela() {
}

double& Tabela::operator()(int i, int j) {
  return this->matriz[i][j];
}

const double& Tabela::operator()(int i, int j) const {
  return this->matriz[i][j];
}

Tabela::operator string() const {
  stringstream out;
  out << "N = " << this->n;
  out << fixed << setprecision(2);
  for (int i = 0; i < this->n; i++) {
    out << "\n| ";
    for (int j = 0; j < this->n; j++) {
      out << this->matriz[i][j] << " ";
    }
    out << "|";
  }
  return out.str();
}

istream& operator>>(istream &in, Tabela& t) {
  string linha, pedaco;
  int n;

  if (!in) {
    __throw_ios_failure("Arquivo inexistente");
  }

  // ler a quantidade de individuos
  // ler "indiv_nb"
  do {
    if (!getline(in, linha)) {
      __throw_ios_failure("Arquivo invalido");
    }
  } while (linha.find("indiv_nb") == string::npos);
  if (sscanf(linha.c_str(), " indiv_nb = %d ", &n) != 1) {
    __throw_ios_failure("Arquivo invalido");
  }

  // cria a tabela
  t = Tabela(n);

  // ler tabela de dissimilaridades
  // ler "DIST_MATRIX= ("
  do {
    if (!getline(in, linha)) {
      __throw_ios_failure("Arquivo invalido");
    }
  } while (linha.find("DIST_MATRIX= (") == string::npos);

  char c;
  double d;

  size_t pos;
  for (int i = 0; i < n; i++) {
    linha.clear();
    while (in >> c && c != '(')
      ;
    linha.push_back(c);
    while (in >> c && c != ')') {
      linha.push_back(c);
    }
    linha.push_back(c);
    while (in >> c && c != ',') {
      linha.push_back(c);
    }

    linha[linha.find("(")] = ' ';
    linha[linha.find(")")] = ' ';

    while ((pos = linha.find(",")) != string::npos)
      linha[pos] = ' ';

    stringstream str(linha);
    for (int j = 0; j <= i; j++) {
      if (!(str >> d)) {
        __throw_ios_failure("Arquivo invalido!!");
      }
      t.matriz[i][j] = t.matriz[j][i] = d;
    }
  }

  /*for (int i = 0; i < t.n; i++) {
    if (!getline(in, linha)) {
    __throw_ios_failure("Arquivo invalido");
    }

    linha[linha.find("(")] = ' ';
    linha[linha.find(")")] = ',';
    is str(linha);

    for (int j = 0; j <= i; j++) {
    if (!(str >> pedaco)) {
    __throw_ios_failure("Arquivo invalido");
    }
    if (sscanf(pedaco.c_str(), " %lf, ", &d) != 1) {
    __throw_ios_failure("Arquivo invalido");
    }
    t.matriz[i][j] = t.matriz[j][i] = d;
    }
    }*/

  return in;
}

ostream& operator<<(ostream& out, const Tabela& t) {
  return out << ((string) t);
}
