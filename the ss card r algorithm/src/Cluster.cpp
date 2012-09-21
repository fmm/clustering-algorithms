#include "Cluster.h"

Cluster::Cluster(size_t n) :
  prototipo(n) {
  }

Cluster::~Cluster() {
}

double Cluster::distancia(size_t individuo, const Tabela& t) const {
  double d = 0;
  for (size_t i = 0; i < this->prototipo.size(); i++) {
    d += t(individuo, prototipo[i]);
  }
  return d;
}

Cluster Cluster::operator&&(const Cluster& c) const {
  Cluster intersecao;
  tr ((*this), iter) {
    if (c.count(*iter)) {
      intersecao.insert(*iter);
    }
  }
  return intersecao;
}

Cluster::operator string() const {
  stringstream out;
  if(this->prototipo.size()) {
    out << this->prototipo << "\n";
  }
  out << "Elemento(s): {";
  for (Cluster::iterator iter = begin(); iter != end(); iter++) {
    if (iter != begin()) {
      out << ",";
    }
    out << " " << (*iter);
  }
  out << " }";
  return out.str();
}

ostream& operator<<(ostream& out, const Cluster& c) {
  return out << (string(c));
}
