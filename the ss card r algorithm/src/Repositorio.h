#ifndef REPOSITORIO_H_
#define REPOSITORIO_H_

#include "Includes.h"
#include "Tabela.h"
#include "Array.h"
#include "Dados.h"

struct Repositorio {

  vector<string> arquivos;
  Array<Tabela> tabela;

  bool rotulado;
  Dados dados;

  Repositorio(vector<string>, double);
  Repositorio(vector<Tabela>&);
  virtual ~Repositorio();

  Tabela& operator[](size_t);
  const Tabela& operator[](size_t) const;
  operator string() const;

  friend ostream& operator<<(ostream&, const Repositorio&);
};

#endif /* REPOSITORIO_H_ */
