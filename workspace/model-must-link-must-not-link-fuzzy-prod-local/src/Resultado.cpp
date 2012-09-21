#include "Resultado.h"

Resultado::Resultado(double alpha, size_t k, size_t p, size_t t, size_t n, size_t m, int id, double J, double CR) {
  this->alpha = alpha;
  this->updateAlpha = alpha < 0;
  for (size_t i = 0; i < k; i++) {
    this->cluster.push_back(Cluster(p));
  }
  this->m = m;
  this->U.resize(n);
  for (size_t i = 0; i < n; i++) {
    this->U[i].resize(k);
  }
  this->coeficiente.resize(k);
  for (size_t i = 0; i < k; i++) {
    this->coeficiente[i].resize(t);
  }
  init(id, J, CR);
}

Resultado::~Resultado() {
}

void Resultado::init(size_t id, double J, double CR) {
  this->id = id;
  this->J = J;
  this->CR = CR;
}

void Resultado::clear() {
  for (size_t i = 0; i < this->cluster.size(); i++) {
    this->cluster[i].clear();
  }
}

void Resultado::srand(const Repositorio& repositorio) {
  const Array<Tabela> &tabela = repositorio.tabela;
  size_t individuos = tabela[0].n;
  size_t prototipos = this->cluster[0].prototipo.size();

  if (this->cluster.size() * prototipos > individuos) {
    __throw_invalid_argument("Quantidade invalida de prototipos");
  }
  map<int, int> repete;
  for (size_t i = 0; i < this->cluster.size(); i++) {
    this->cluster[i].prototipo.clear();
    for (size_t j = 0; j < prototipos; j++) {
      int novoPrototipo;
      do {
        novoPrototipo = (int) (((double) individuos * rand()) / (RAND_MAX+1.0));
        novoPrototipo %= individuos;
      } while (repete[novoPrototipo]++);
      this->cluster[i].prototipo.push_back(novoPrototipo);
    }
  }

  // inicializa coeficientes
  for (size_t i = 0; i < this->coeficiente.size(); i++) {
    for (size_t j = 0; j < this->coeficiente[i].size(); j++) {
      this->coeficiente[i][j] = 1.0;
    }
  }

  // inicializa matriz U
  this->atualizarUInicializacao(repositorio);

  // inicializa J
  this->J = DBL_MAX;
  this->atualizaJ(repositorio);

  // inicializa CR
  this->CR = 0;
  this->atualizaCR(repositorio);
}

void Resultado::atualizaCluster(const Repositorio& repositorio) {
  const Array<Tabela> &tabela = repositorio.tabela;
  size_t individuos = tabela[0].n;
  for (size_t i = 0; i < this->cluster.size(); i++) {
    vector<pair<double, int> > v;
    for (size_t j = 0; j < individuos; j++) {
      double d = 0;
      for (size_t k = 0; k < individuos; k++) {
        for (size_t t = 0; t < tabela.size(); t++) {
          d += (tabela[t](k, j) * pow(this->U[k][i], m) * this->coeficiente[i][t]); // TODO
        }
      }
      v.push_back(mp(d,j));
    }
    sort(v.begin(), v.end(), less<pair<double, int> > ());
    for (size_t j = 0; j < this->cluster[i].prototipo.size(); j++) {
      this->cluster[i].prototipo[j] = v[j].second;
    }
  }
}

bool Resultado::atualizaCoeficiente(const Repositorio& repositorio) {
  const Array<Tabela> &tabela = repositorio.tabela;
  size_t individuos = tabela[0].n;
  double num = 1, den[tabela.size()];

  for(size_t i = 0; i < this->cluster.size(); i++) {
    num = 1;
    for(size_t j = 0; j < tabela.size(); j++) {
      den[j] = 0;
      for(size_t k = 0; k < individuos; k++) {
        den[j] += cluster[i].distancia(k,tabela[j]) * pow(this->U[k][i],m);
      }
      num *= std::pow(den[j], 1.0 / tabela.size());
    }

    if( cmp(num) ) {
      for (size_t j = 0; j < tabela.size(); j++) {
        assert(den[i]);
        this->coeficiente[i][j] = num / den[j];
      }
    }
    else {
      Alert("division by zero during coefficient calculation");
      return false;
    }
  }

  return true;
}

void Resultado::atualizarUInicializacao(const Repositorio& repositorio) {
  const Array<Tabela> &tabela = repositorio.tabela;
  size_t individuos = tabela[0].n;
  double dist[this->cluster.size()];
  for (size_t i = 0; i < individuos; i++) {
    for (size_t j = 0; j < this->cluster.size(); j++) {
      dist[j] = 0;
      for (size_t t = 0; t < tabela.size(); t++) {
        dist[j] += (this->cluster[j].distancia(i, tabela[t]) * this->coeficiente[j][t]); // TODO
      }
    }

    vector<size_t> v;
    for (size_t j = 0; j < this->cluster.size(); j++) {
      if (cmp(dist[j]) == 0) {
        v.push_back(j);
      }
    }

    if (v.size()) {
      for (size_t j = 0; j < this->cluster.size(); j++) {
        this->U[i][j] = 0;
      }
      for (size_t j = 0; j < v.size(); j++) {
        this->U[i][v[j]] = (1.0 / v.size());
      }
    } else {
      for (size_t j = 0; j < this->cluster.size(); j++) {
        this->U[i][j] = 0;
        for (size_t k = 0; k < this->cluster.size(); k++) {
          this->U[i][j] += pow(dist[j] / dist[k], 1.0 / (m - 1));
        }
        assert(cmp(this->U[i][j]));
        this->U[i][j] = 1.0 / this->U[i][j];
      }
    }
  }

  this->clear();
  for (size_t i = 0; i < individuos; i++) {
    vector<pair<double, int> > v;
    for (size_t j = 0; j < this->cluster.size(); j++) {
      v.push_back(mp(this->U[i][j],j));
    }
    sort(v.begin(), v.end(), greater<pair<double, int> > ());
    this->cluster[v[0].second].insert(i);
  }
}

void Resultado::atualizarU(const Repositorio& repositorio) {
  const Array<Tabela> &tabela = repositorio.tabela;
  int individuos = tabela[0].n;
  assert(repositorio.rotulado);
  const set< pair<int,int> > &mustLink = repositorio.dados.mustLink;
  const set< pair<int,int> > &mustNotLink = repositorio.dados.mustNotLink;
  if(updateAlpha) {
    // FIXME
  }
  double dist[this->cluster.size()];
  Array< Array<double> > newU = this->U;
  for (int i = 0; i < individuos; i++) {
    // MFCMdd-RWL-P
    for (size_t j = 0; j < this->cluster.size(); j++) {
      dist[j] = 0;
      for (size_t t = 0; t < tabela.size(); t++) {
        dist[j] += (this->cluster[j].distancia(i, tabela[t]) * this->coeficiente[j][t]); // TODO
      }
    }
    vector<size_t> v;
    for (size_t j = 0; j < this->cluster.size(); j++) {
      if (cmp(dist[j]) == 0) {
        v.push_back(j);
      }
    }
    if (v.size()) {
      for (size_t j = 0; j < this->cluster.size(); j++) {
        newU[i][j] = 0;
      }
      for (size_t j = 0; j < v.size(); j++) {
        newU[i][v[j]] = (1.0 / v.size());
      }
    } else {
      for (size_t j = 0; j < this->cluster.size(); j++) {
        newU[i][j] = 0;
        for (size_t k = 0; k < this->cluster.size(); k++) {
          newU[i][j] += pow(dist[j] / dist[k], 1.0 / (m - 1));
        }
        assert(cmp(newU[i][j]));
        newU[i][j] = 1.0 / newU[i][j];
      }
      // CONST
      vector<double> D(this->cluster.size(),0.0);
      for(size_t j = 0; j < this->cluster.size(); j++) {
        for(size_t t = 0; t < tabela.size(); t++) {
          D[j] += this->cluster[j].distancia(i, tabela[t]) * this->coeficiente[j][t];  
        }
      }
      for(size_t j = 0; j < this->cluster.size(); j++) {
        double num = 0, den = 0;
        for(size_t h = 0; h < this->cluster.size(); h++) {
          double tmp = 0;
          // must link
          for(set< pair<int,int> >::const_iterator x = mustLink.begin(); x != mustLink.end(); x++) {
            int l = x->first, m = x->second;
            if(l != i) continue;
            tmp += this->U[m][j] - this->U[m][h];
          }
          // must not link
          for(set< pair<int,int> >::const_iterator x = mustNotLink.begin(); x != mustNotLink.end(); x++) {
            int l = x->first, m = x->second;
            if(l != i) continue;
            tmp += this->U[m][h] - this->U[m][j];
          }
          assert( cmp(D[h]) );
          num += tmp / 2 / D[h];
          den += D[j] / D[h];
        }
        assert(cmp(den));
        newU[i][j] += this->alpha * num / den;
      }
    }
  }

  for(int i = 0; i < individuos; i++) {
    for(size_t j = 0; j < this->cluster.size(); j++) {
      this->U[i][j] = newU[i][j];
    }
  }

  bool normalize = false;
  rep(i,U.size()) rep(j,U[i].size()) if(U[i][j] < 0) U[i][j] = 0, normalize = true; else if(U[i][j] > 1.0) U[i][j] = 1.0, normalize = true;

  if(normalize) {
    Alert("Fixing U membership.");
    rep(t,individuos) {
      double sum = 0;
      rep(c,this->cluster.size()) sum += U[t][c];
      assert(cmp(sum));
      rep(c,this->cluster.size()) U[t][c] /= sum;
    }
  }

  this->clear();
  for (int i = 0; i < individuos; i++) {
    vector<pair<double, int> > v;
    for (size_t j = 0; j < this->cluster.size(); j++) {
      v.push_back(mp(this->U[i][j],j));
    }
    sort(v.begin(), v.end(), greater<pair<double, int> > ());
    this->cluster[v[0].second].insert(i);
  }

}

double Resultado::atualizaJ(const Repositorio& repositorio) {
  const Array<Tabela> &tabela = repositorio.tabela;
  size_t individuos = tabela[0].n;

  double J1 = 0;
  for (size_t i = 0; i < this->cluster.size(); i++) {
    for (size_t j = 0; j < individuos; j++) {
      double dist = 0;
      for (size_t t = 0; t < tabela.size(); t++) {
        dist += (this->cluster[i].distancia(j, tabela[t]) * this->coeficiente[i][t]); // TODO
      }
      J1 += ((dist * pow(this->U[j][i], m)));
    }
  }

  assert(repositorio.rotulado);
  const set< pair<int,int> > &mustLink = repositorio.dados.mustLink;
  const set< pair<int,int> > &mustNotLink = repositorio.dados.mustNotLink;

  double J2 = 0;
  // must link
  for(set< pair<int,int> >::const_iterator x = mustLink.begin(); x != mustLink.end(); x++) {
    int l = x->first, m = x->second;
    for(unsigned int r=0;r<this->cluster.size();r++) for(unsigned int s=0;s<this->cluster.size();s++) if(r!=s) {
      J2 += this->U[l][r] * this->U[m][s];
    }
  }
  // must not link
  for(set< pair<int,int> >::const_iterator x = mustNotLink.begin(); x != mustNotLink.end(); x++) {
    int l = x->first, m = x->second;
    for(unsigned int r=0;r<this->cluster.size();r++) {
      J2 += this->U[l][r] * this->U[m][r];
    }
  }

  double novoJ = J1 + this->alpha * J2;

  if( cmp(novoJ,this->J) > 0 ) {
    Alert("problem updating function J, it didn't decrease");
    Msg("Old J: " << this->J);
    Msg("New J: " << novoJ);
  }

  if(false) assert(cmp(novoJ,this->J) <= 0);  return (this->J = novoJ);
  
  return (this->J = novoJ);
}

double Resultado::atualizaCR(const Repositorio& repositorio) {
  if (repositorio.rotulado) {
    this->CR = repositorio.dados.calculaCR(this->cluster);
  }
  return this->CR;
}

Cluster& Resultado::operator[](size_t __n) {
  return this->cluster[__n];
}

const Cluster& Resultado::operator[](size_t __n) const {
  return this->cluster[__n];
}

double& Resultado::operator()(size_t __n, size_t __m) {
  return this->U[__n][__m];
}

const double& Resultado::operator()(size_t __n, size_t __m) const {
  return this->U[__n][__m];
}

bool Resultado::operator<(const Resultado& r) const {
  return this->J < r.J;
}

Resultado::operator string() const {
  stringstream out;
  _dbg(out,id);
  _dbg(out,J);
  _dbg(out,CR);
  _dbg(out,cluster.size());
  for (size_t i = 0; i < cluster.size(); i++) {
    _dbg(out,i);
    _print(out,cluster[i]);
  }
  return out.str();
}

ostream& operator<<(ostream& out, const Resultado& r) {
  return out << (string(r));
}
