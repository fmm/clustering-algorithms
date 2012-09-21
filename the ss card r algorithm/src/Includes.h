#ifndef INCLUDES_H_
#define INCLUDES_H_

#ifdef WINDOWS
#pragma comment(linker, "/STACK:16777216")
#endif

using namespace std;

#include <bits/stdc++.h>

#define rep(x,n) for(int x = 0; x < int(n); ++x)
#define tr(container,it) for(typeof(container.begin()) it = container.begin(); it != container.end(); it++)

#define _print(y,x) y << x << endl
#define print(x) _print(cout,x)
#define _dbg(y,x) y << #x << " == " << x << endl
#define dbg(x) _dbg(cerr,x)
#define _ << " , " <<

#define Print(x) _print(cerr,x)
#define Dbg(x) _dbg(cerr,x)

#define Msg(x) Print("*** " << x);
#define Alert(x) Print("!!! " << x << " !!!");

#define mp make_pair
#define x first
#define y second

#define pi M_PI

#define cast(x,t) *({stringstream ss;static t __ret;ss<<x,ss>>__ret;&__ret;})

typedef long long ll;
typedef unsigned long long ull;

const int inf = ~0U>>1;
const ll INF = ~0ULL>>1;
const double eps = 1e-7;

int cmp(double x, double y = 0);
double sqr(double x);

enum Funcao {
  ABS, POW
};

enum Relacao {
  INDIVIDUAL, GRUPO
};

enum Coeficiente {
  INDIVIDUO, TABELA, TABELA_CLUSTER
};

#endif /* INCLUDES_H_ */
