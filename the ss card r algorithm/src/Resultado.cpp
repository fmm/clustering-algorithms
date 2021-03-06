#include "Resultado.h"

Resultado::Resultado(double alpha, size_t c, size_t n, size_t q, int id, double J, double CR) {
  this->alpha = alpha;
  for (size_t i = 0; i < c; i++) {
    this->cluster.push_back(Cluster(0));
  }
  this->m = 2;
  this->q = q;

  this->U = Array< Array<double> >(c, Array<double>(n,0.0));

  init(id, J, CR);
  updateAlpha = alpha < 0;
  if(updateAlpha)
    Msg("alpha will be update automatically in each iteration");
}

Resultado::~Resultado() {
}

void Resultado::srand(const Repositorio& repositorio) {
  const Array<Tabela> &tabela = repositorio.tabela;
  size_t n = tabela[0].n;
  size_t c = this->cluster.size();
  size_t s = tabela.size();

  rep(i,U.size()) rep(j,U[i].size()) U[i][j]=rand();

  rep(j,n) {
    double sum = 0;
    rep(cc,c) sum += U[cc][j];
    rep(cc,c) U[cc][j] /= sum;
  }

  rep(cc,c) cluster[cc].clear();

  rep(i,n) {
    int x = -1;
    rep(cc,c) {
      if(x == -1 || U[cc][i] > U[x][i]) x = cc;
    }
    cluster[x].insert(i);
  }

  this->beta = 0;
  this->W = Array< Array<double> >(c, Array<double>(s,1.0/s));

  atualizaJ(repositorio);
  atualizaCR(repositorio);
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

bool Resultado::atualizarUW(const Repositorio& repositorio) {
  const Array<Tabela> &tabela = repositorio.tabela;
  size_t n = tabela[0].n;
  size_t c = this->cluster.size();
  size_t s = tabela.size();

  assert(repositorio.rotulado);
  const set< pair<int,int> > &mustLink = repositorio.dados.mustLink;
  const set< pair<int,int> > &mustNotLink = repositorio.dados.mustNotLink;

  // update alpha
  if(updateAlpha) {
    double J1 = 0;
    rep(cc,c) {
      double num = 0;
      rep(j,n) rep(k,n) {
        double val = 0;
        rep(ss,s) val += pow(W[cc][ss],q) * tabela[ss](j,k);
        num += pow(U[cc][j],2.0) * pow(U[cc][k],2.0) * val;
      }
      double den = 0;
      rep(j,n) den += pow(U[cc][j],2.0);
      den *= 2;
      assert(cmp(den));
      J1 += num / den;
    }
    double J2 = 0;
    // must link
    for(set< pair<int,int> >::const_iterator x = mustLink.begin(); x != mustLink.end(); x++) {
      int l = x->first, m = x->second;
      for(unsigned int r=0;r<this->cluster.size();r++) for(unsigned int s=0;s<this->cluster.size();s++) if(r!=s) {
        J2 += this->U[r][l] * this->U[s][m];
      }
    }
    // must not link
    for(set< pair<int,int> >::const_iterator x = mustNotLink.begin(); x != mustNotLink.end(); x++) {
      int l = x->first, m = x->second;
      for(unsigned int r=0;r<this->cluster.size();r++) {
        J2 += this->U[r][l] * this->U[r][m];
      }
    }
    //assert(cmp(J2));
    if(!cmp(J2)) {
      Alert("attempt of division by zero, alpha now is 1.0");
      alpha = 1;
    }
    else {
      alpha = J1 / J2;
    }
    Msg("alpha now is: " << alpha);
  }

  typedef Array< double > V;
  typedef Array< V > M;

  __typeof(U) Upow2 = U;
  rep(i,Upow2.size()) rep(j,Upow2[i].size()) Upow2[i][j] *= U[i][j];

  // equation [5]
  Array< M > R(c, M(n,V(n,0)));
  rep(i,n) rep(j,n) rep(cc,c) rep(ss,s) {
    R[cc][i][j] += pow(W[cc][ss],q) * tabela[ss](i,j);
  }

  // equation [8]
  Array< M > Rbeta(c, M(n,V(n,0)));
  rep(cc,c) rep(i,n) rep(j,n) {
    Rbeta[cc][i][j] += R[cc][i][j] + beta * ((i != j) - (i == j));
  }

  // equation [6]
  Array< V > v(c, V(n,0));
  rep(cc,c) {
    double den=0;
    rep(j,n) v[cc][j] = Upow2[cc][j], den += v[cc][j];
    assert(cmp(den));
    rep(j,n) v[cc][j] /= den;
  }

  // equation [5]
  M d(c,V(n,0));
  rep(cc,c) {
    V A(n,0);
    rep(i,n) rep(j,n) A[i] += Rbeta[cc][i][j] * v[cc][j];
    double B=0;
    rep(i,n) B += v[cc][i] * A[i] / 2;
    rep(i,n) d[cc][i] = A[i] - B;
  }

  double delta=0;
  rep(cc,c) rep(j,n) if(cmp(d[cc][j]) < 0) {
    double norm = 0;
    rep(i,n) norm += pow(v[cc][i] - (j==i),2);
    delta = max(delta, -2 * d[cc][j] / norm);
  }

  beta += delta;
  delta /= 2;

  if(cmp(delta)>0){
    rep(cc,c) rep(j,n) {
      double norm = 0;
      rep(i,n) norm += pow(v[cc][i] - (j==i),2);
      d[cc][j] += delta * norm;
    }
  }

  // equation [7]
  rep(cc,c) rep(j,n) {
    if(cmp(d[cc][j]) > 0) {
      double val = 0.0;
      rep(k,c) {
        if(cmp(d[k][j]) == 0) {
          Alert("attempt of division by zero (d[k][j] was zero),  fucntion \"atualizarUW\" is being finished");
          return false;
        }
        val += d[cc][j] / d[k][j];
      }
      assert(cmp(val) > 0);
      U[cc][j] = 1.0 / val;
    }
    else {
      U[cc][j] = 0.0;
    }
  }

  __typeof(U) Ucte = U, a = U, C = U;

  rep(v,a.size()) rep(t,a[v].size()) {
    double num = 0, den = 0;
    rep(k,n) {
      double tmp = 0;
      rep(ss,s) tmp += pow(W[v][ss],q) * tabela[ss](t,k);
      num += Upow2[v][k] * tmp;
      den += Upow2[v][k];
    }
    assert( cmp(den) );
    a[v][t] = 2 * num / den;
  }

  rep(v,C.size()) rep(t,C[v].size()) {
    C[v][t] = 0;
    for(set< pair<int,int> >::const_iterator x = mustLink.begin(); x != mustLink.end(); x++) {
      int tt = x->first, k = x->second;
      if(tt != t) continue;
      rep(l,c) if(l != v) C[v][t] += sqrt(Upow2[l][k]);
    }
    for(set< pair<int,int> >::const_iterator x = mustNotLink.begin(); x != mustNotLink.end(); x++) {
      int tt = x->first, k = x->second;
      if(tt != t) continue;
      C[v][t] += sqrt(Upow2[v][k]);
    }
  }

  Array< double > Ct(n,0.0);

  rep(t,Ct.size()) {
    double num = 0, den = 0;
    rep(i,c) {
      assert( cmp(a[i][t]) );
      num += C[i][t] / a[i][t];
      den += 1.0 / a[i][t];
    }
    assert( cmp(den) );
    Ct[t] = num / den;
  }

  rep(v,Ucte.size()) rep(t,Ucte[v].size()) {
    assert( cmp(a[v][t]) );
    Ucte[v][t] = alpha * (Ct[t] - C[v][t]) / a[v][t];
  }

  rep(i,U.size()) rep(j,U[i].size()) U[i][j] += Ucte[i][j];

  // fix clusters
  rep(cc,c) cluster[cc].clear();

  rep(i,n) {
    int x = -1;
    rep(cc,c) if(x == -1 || U[cc][i] > U[x][i]) x = cc;
    cluster[x].insert(i);
  }

  bool normalize = false;
  rep(i,U.size()) rep(j,U[i].size()) if(U[i][j] < 0) U[i][j] = 0, normalize = true; else if(U[i][j] > 1.0) U[i][j] = 1.0, normalize = true;

  if(normalize) {
    Alert("Fixing U membership.");
    rep(t,n) {
      double sum = 0;
      rep(cc,c) sum += U[cc][t];
      assert(cmp(sum));
      rep(cc,c) U[cc][t] /= sum;
    }
  }

  rep(i,U.size()) rep(j,U[i].size()) Upow2[i][j] = U[i][j] * U[i][j];

  // equation [22]
  M D(c,V(s,0));
  rep(cc,c) rep(ss,s) rep(j,n) rep(k,n) D[cc][ss] += Upow2[cc][j] * Upow2[cc][k] * tabela[ss](j,k);

  rep(cc,c) rep(ss,s) {
    double val = 0.0;
    rep(pp,s) {
      assert(cmp(D[cc][pp]));
      val += pow(D[cc][ss] / D[cc][pp], 1.0 / (q - 1));
    }
    if(cmp(val)) {
      W[cc][ss] = 1.0 / val;
    }
    else {
      Alert("attempt of division by zero (sum of D[cc][ss] was zero),  fucntion \"atualizarUW\" is being finished");
      return false;
    }
  }
  return true;
}

double Resultado::atualizaJ(const Repositorio& repositorio) {
  const Array<Tabela> &tabela = repositorio.tabela;
  size_t n = tabela[0].n;
  size_t c = this->cluster.size();
  size_t s = tabela.size();

  double J1 = 0;
  rep(cc,c) {
    double num = 0;
    rep(j,n) rep(k,n) {
      double val = 0;
      rep(ss,s) val += pow(W[cc][ss],q) * tabela[ss](j,k);
      num += pow(U[cc][j],2.0) * pow(U[cc][k],2.0) * val;
    }
    double den = 0;
    rep(j,n) den += pow(U[cc][j],2.0);
    den *= 2;
    assert(cmp(den));
    J1 += num / den;
  }

  assert(repositorio.rotulado);
  const set< pair<int,int> > &mustLink = repositorio.dados.mustLink;
  const set< pair<int,int> > &mustNotLink = repositorio.dados.mustNotLink;

  double J2 = 0;
  // must link
  for(set< pair<int,int> >::const_iterator x = mustLink.begin(); x != mustLink.end(); x++) {
    int l = x->first, m = x->second;
    for(unsigned int r=0;r<this->cluster.size();r++) for(unsigned int s=0;s<this->cluster.size();s++) if(r!=s) {
      J2 += this->U[r][l] * this->U[s][m];
    }
  }
  // must not link
  for(set< pair<int,int> >::const_iterator x = mustNotLink.begin(); x != mustNotLink.end(); x++) {
    int l = x->first, m = x->second;
    for(unsigned int r=0;r<this->cluster.size();r++) {
      J2 += this->U[r][l] * this->U[r][m];
    }
  }

  double novoJ = J1 + this->alpha * J2;

  if( cmp(novoJ,this->J) > 0 ) {
    Alert("problem updating function J, it didn't decrease");
    Msg("Old J: " << this->J);
    Msg("New J: " << novoJ);
  }

  if(false) assert(cmp(novoJ,this->J) <= 0);

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
