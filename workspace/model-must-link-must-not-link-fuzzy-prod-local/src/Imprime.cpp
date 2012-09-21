#include "Imprime.h"

#define DISTANCIA 5
#define DIVISA '|'
#define TRACO '-'
void printLinha1(ostream &out) {
  out << DIVISA;
  for (int i = 0; i < 38; i++) {
    out << TRACO;
  }
  for (int i = 6; i < 6 * DISTANCIA; i++) {
    out << TRACO;
  }
  out << DIVISA << '\n';
}

void printLinha2(ostream &out) {
  out << DIVISA;
  for (int i = 0; i < 11; i++) {
    out << TRACO;
  }
  for (int i = 2; i < 2 * DISTANCIA; i++) {
    out << TRACO;
  }
  out << DIVISA;
  for (int i = 0; i < 15; i++) {
    out << TRACO;
  }
  for (int i = 2; i < 2 * DISTANCIA; i++) {
    out << TRACO;
  }
  out << DIVISA;
  for (int i = 0; i < 10; i++) {
    out << TRACO;
  }
  for (int i = 2; i < 2 * DISTANCIA; i++) {
    out << TRACO;
  }
  out << DIVISA << '\n';
}

void printInicializacaoHeader(ostream &out, int inicializacao) {
  out << DIVISA;
  out << " # Initialization: ";
  out << setw(19 + (DISTANCIA - 1) * 6) << left << inicializacao;
  out << DIVISA << '\n';
}
void print2(ostream &out) {
  out << DIVISA;
  out << setw(DISTANCIA) << " ";

  out << "Iteration";

  out << setw(DISTANCIA) << " ";
  out << DIVISA;
  out << setw(DISTANCIA) << " ";

  out << setw(6) << " ";
  out << "J";
  out << setw(6) << " ";

  out << setw(DISTANCIA) << " ";
  out << DIVISA;
  out << setw(DISTANCIA) << " ";

  out << setw(3) << " ";
  out << "CR";
  out << setw(3) << " ";

  out << setw(DISTANCIA) << " ";
  out << DIVISA << '\n';
}
void print3(ostream &out, int iteracao, double J, double CR) {


  out << DIVISA;
  out << setw(DISTANCIA) << " ";

  out << setfill('0');
  out << setw(9) << right << iteracao;
  out << setfill(' ');

  out << setw(DISTANCIA) << " ";
  out << DIVISA;
  out << setw(DISTANCIA) << " ";

  out << scientific << J;

  out << setw(DISTANCIA) << " ";
  out << DIVISA;
  out << setw(DISTANCIA) << " ";

  out << fixed << setprecision(6) << CR;

  out << setw(DISTANCIA) << " ";
  out << DIVISA << '\n';

}
#undef DISTANCIA
#define DISTANCIA 2

/*
 * Funcoes para imprimir a tabela 2
 */

void printLinha1(ostream &out, int clus) {
  out << DIVISA;
  for (int i = 0; i < 43; i++) {
    out << TRACO;
  }
  for (int i = 0; i < 15 * clus; i++) {
    out << TRACO;
  }
  int espaco = DISTANCIA * (3 + clus);
  for (int i = 3 + clus; i < espaco; i++) {
    out << TRACO << TRACO;
  }
  out << DIVISA << '\n';
}

void printLinha2(ostream &out, int clus) {
  out << DIVISA;
  for (int i = 0; i < 9; i++) {
    out << TRACO;
  }
  for (int i = 2; i < 2 * DISTANCIA; i++) {
    out << TRACO;
  }
  out << DIVISA;

  for (int j = 0; j < clus; j++) {
    for (int i = 0; i < 14; i++) {
      out << TRACO;
    }
    for (int i = 2; i < 2 * DISTANCIA; i++) {
      out << TRACO;
    }
    out << DIVISA;
  }

  for (int i = 0; i < 14; i++) {
    out << TRACO;
  }
  for (int i = 2; i < 2 * DISTANCIA; i++) {
    out << TRACO;
  }
  out << DIVISA;

  for (int i = 0; i < 18; i++) {
    out << TRACO;
  }
  for (int i = 2; i < 2 * DISTANCIA; i++) {
    out << TRACO;
  }
  out << DIVISA << '\n';
}

void printRelacaoHeader(ostream &out, int clus) {
  out << DIVISA;
  out << setw(DISTANCIA) << " ";

  out << "Pattern";

  out << setw(DISTANCIA) << " ";
  out << DIVISA;

  for (int i = 0; i < clus; i++) {
    out << setw(DISTANCIA) << " ";
    out << "cluster[";
    out << setfill('0');
    out << setw(3) << (i + 1);
    out << setfill(' ');
    out << "]";
    out << setw(DISTANCIA) << " ";
    out << DIVISA;
  }

  out << setw(DISTANCIA) << " ";

  out << "Hard Cluster";

  out << setw(DISTANCIA) << " ";

  out << DIVISA;
  out << setw(DISTANCIA) << " ";

  out << "A Priori Cluster";

  out << setw(DISTANCIA) << " ";

  out << DIVISA << '\n';
}

#include "Array.h"
void print2(ostream &out, int p, int clus, Array<double>& a, int select, int priori) {
  out << DIVISA;
  out << setw(DISTANCIA) << " ";

  out << setfill('0');
  out << setw(7) << right << p;
  out << setfill(' ');

  out << setw(DISTANCIA) << " ";

  for (int i = 0; i < clus; i++) {
    out << DIVISA;
    out << setw(DISTANCIA) << " ";

    out << setw(3) << " ";
    out << fixed << setprecision(4) << a[i];
    out << setw(3) << " ";

    out << setw(DISTANCIA) << " ";
  }

  out << DIVISA;
  out << setw(DISTANCIA) << " ";

  out << setw(5) << " ";
  out << setfill('0');
  out << setw(3) << (select + 1);
  out << setfill(' ');
  out << setw(4) << " ";

  out << setw(DISTANCIA) << " ";
  out << DIVISA;

  out << setw(DISTANCIA) << " ";

  out << setw(7) << " ";
  out << setfill('0');
  out << setw(3) << (priori + 1);
  out << setfill(' ');
  out << setw(6) << " ";

  out << setw(DISTANCIA) << " ";
  out << DIVISA << '\n';
}
#undef DISTANCIA
#undef DIVISA
#undef TRACO
