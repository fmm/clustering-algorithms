#include "Repositorio.h"

Repositorio::Repositorio(vector<string> arquivos, double label) : arquivos(arquivos) {
  this->rotulado = (Dados::var_classe > 0);
  for (size_t i = 0; i < arquivos.size(); i++) {
    this->tabela.push_back(Tabela(arquivos[i]));
  }
  if (this->rotulado) {
    try {
      this->dados = Dados(arquivos[0]);
      this->dados.generateTables(label);
    } catch (...) {
      this->dados = Dados();
      this->rotulado = false;
    }
  }
}

Repositorio::Repositorio(vector<Tabela>& tabelas) {
  for (size_t i = 0; i < tabelas.size(); i++) {
    stringstream out;
    out << "tabela(" << i << ")";
    this->arquivos.push_back(out.str());
    this->tabela.push_back(tabelas[i]);
  }
}

Repositorio::~Repositorio() {
}

Tabela& Repositorio::operator[](size_t __n) {
  return this->tabela[__n];
}

const Tabela& Repositorio::operator[](size_t __n) const {
  return this->tabela[__n];
}

Repositorio::operator string() const {
  stringstream out;
  out << "Arquivos[" << this->arquivos.size() << "] = {";
  for (size_t i = 0; i < this->arquivos.size(); i++) {
    if (i > 0) {
      out << ",";
    }
    out << " " << this->arquivos[i];
  }
  out << " }";
  return out.str();
}

ostream& operator<<(ostream& out, const Repositorio& r) {
  return out << ((string) r);
}
