#ifndef DADOS_H_
#define DADOS_H_

#include "Includes.h"
#include "Array.h"
#include "Cluster.h"

struct Dados {

  static int var_classe;

  Array<Cluster> prioriCluster;
  set< pair<int,int> > mustLink, mustNotLink;

  Dados();
  Dados(string);
  virtual ~Dados();

  double calculaCR(Array<Cluster>&) const;

  Cluster& operator[](size_t);
  const Cluster& operator[](size_t) const;
  operator string() const;

  void generateTables(double);

  friend istream& operator>>(istream&, Dados&);
  friend ostream& operator<<(ostream&, const Dados&);
};

#endif /* DADOS_H_ */
