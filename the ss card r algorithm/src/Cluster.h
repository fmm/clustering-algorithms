#ifndef CLUSTER_H_
#define CLUSTER_H_

#include "Includes.h"
#include "Array.h"
#include "Tabela.h"

struct Cluster: public set<int, less<int> > {

  struct Prototipo: public Array<int> {
    Prototipo(size_t n = 0) {
      resize(n);
    }

    operator string() const {
      stringstream out;
      out << "Prototipo(s): {";
      for (size_t i = 0; i < this->size(); i++) {
        if(i > 0) {
          out << ",";
        }
        out << " " << (*this)[i];
      }
      out << " }";
      return out.str();
    }

    friend ostream& operator<<(ostream& out, const Prototipo& p) {
      return out << ((string) p);
    }
  } prototipo;

  Cluster(size_t = 0);
  virtual ~Cluster();

  double distancia(size_t, const Tabela&) const;

  Cluster operator&&(const Cluster&) const;
  operator string() const;

  friend ostream& operator<<(ostream&, const Cluster&);
};

#endif /* CLUSTER_H_ */
