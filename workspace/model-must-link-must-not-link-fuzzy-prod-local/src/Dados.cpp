#include "Dados.h"

int Dados::var_classe = 0;

Dados::Dados() {

}

Dados::Dados(string arquivo) {
  ifstream in(arquivo.c_str(), ios::in);
  in >> *this;
  in.close();
}

Dados::~Dados() {
}

#define comb(x) (((x) * ((x) - 1)) / 2)
double Dados::calculaCR(Array<Cluster>& cluster) const {
  size_t individuos = 0;
  for (size_t i = 0; i < cluster.size(); i++) {
    individuos += cluster[i].size();
  }

  size_t k = cluster.size();
  size_t p = this->prioriCluster.size();

  size_t table[k + 1][p + 1];

  // table[i][j] = cluster[i] && prioriCluster[j]
  for (size_t i = 0; i < k; i++) {
    for (size_t j = 0; j < p; j++) {
      table[i][j] = (cluster[i] && this->prioriCluster[j]).size();
    }
  }

  // table[i][p] = cluster[i].tamanho
  for (size_t i = 0; i < k; i++) {
    table[i][p] = cluster[i].size();
  }

  // table[k][i] = prioriCluster[i].tamanho
  for (size_t i = 0; i < p; i++) {
    table[k][i] = this->prioriCluster[i].size();
  }

  double termo[4], temp[2];
  double pot = std::pow(comb(individuos), -1.0);

  // termo 1
  termo[1] = 0;
  for (size_t i = 0; i < k; i++) {
    for (size_t j = 0; j < p; j++) {
      termo[1] += comb(table[i][j]);
    }
  }

  // temps
  temp[0] = 0;
  for (size_t i = 0; i < k; i++) {
    temp[0] += comb(table[i][p]);
  }
  temp[1] = 0;
  for (size_t i = 0; i < p; i++) {
    temp[1] += comb(table[k][i]);
  }

  // termo 2
  termo[2] = pot * (temp[0] * temp[1]);

  // termo 3
  termo[3] = 0.5 * (temp[0] + temp[1]);

  return (termo[1] - termo[2]) / (termo[3] - termo[2]);
}

Cluster& Dados::operator[](size_t __n) {
  return this->prioriCluster[__n];
}

const Cluster& Dados::operator[](size_t __n) const {
  return this->prioriCluster[__n];
}

void Dados::generateTables(double taxe) {
  map<int,int> label;
  for(unsigned int i = 0; i < prioriCluster.size(); i++) {
    Cluster &c = prioriCluster[i];
    for(Cluster::iterator j = c.begin(); j != c.end(); j++) {
      label[*j] = i;
    }
  }
  int n = label.size();
  mustLink.clear();
  mustNotLink.clear();
  int m = int(n * taxe + 0.7777777);
  vector<int> object(n);
  rep(i,n) object[i] = i;
  random_shuffle(object.begin(),object.end());
  assert(m <= n);
  rep(i,m) rep(j,m) {
    int a = object[i], b = object[j];
    if(label[a] == label[b]) {
      mustLink.insert( mp(a,b) );
    }
    else {
      mustNotLink.insert( mp(a,b) );
    }
  }

}

Dados::operator string() const {
  stringstream out;
  _dbg(out,prioriCluster.size());
  for (size_t i = 0; i < prioriCluster.size(); i++) {
    _dbg(out,i);
    _print(out,prioriCluster[i]);
  }
  return out.str();
}

istream& operator>>(istream& in, Dados& a) {
  string linha, pedaco;
  int n;
  map<string, int> hash;
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

  // ler tabela de descricao de variaveis
  // ler "RECTANGLE_MATRIX = ("
  do {
    if (!getline(in, linha)) {
      __throw_ios_failure("Arquivo invalido");
    }
  } while (linha.find("RECTANGLE_MATRIX = (") == string::npos);

  char c = '?';
  int ct;
  size_t priori = 0;
  size_t pos;
  vector<string> v;
  string piece;

  for (int i = 0; i < n; i++) {

    linha.clear();

    ct = 0;
    while (!ct) {
      if (!(in >> c))
        __throw_ios_failure("Arquivo invalido");
      if (c == '(')
        ct++;
    }

    while (ct) {
      if (!(in >> c))
        __throw_ios_failure("Arquivo invalido");
      if (c == '(')
        ct++;
      if (c == ')')
        ct--;
      if (ct)
        linha.push_back(c);
    }

    while ((pos = linha.find(' ')) != string::npos)
      linha.erase(pos, 1);

    ct = 0;
    piece.clear();
    v.clear();
    for (size_t j = 0; j < linha.size(); j++) {
      if (linha[j] == '(')
        ct++;
      if (linha[j] == ')')
        ct--;
      piece.push_back(linha[j]);
      if (!ct) {
        if (piece[0] != ',')
          v.push_back(piece);
        piece.clear();
      }
    }

    piece = v[Dados::var_classe - 1];
    assert(sscanf(piece.c_str(),"%u",&priori) == 1);
    assert(priori> 0);

    for (size_t j = a.prioriCluster.size(); j < priori; j++) {
      a.prioriCluster.push_back(Cluster(0));
    }
    a.prioriCluster[priori - 1].insert(i);
  }

  return in;
}

ostream& operator<<(ostream& out, const Dados& a) {
  return out << ((string) a);
}
