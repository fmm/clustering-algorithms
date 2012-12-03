#if 0
#define EXTREMELY_DETAILED
#define DETAILED
#define PROTOTYPE_VARIANCE
#define MAGIC_VARIANCE
#define RUN_ANOTHER
#endif

#ifdef DETAILED
#define EXTREMELY_DETAILED
#endif

#include <bits/stdc++.h>
#include <sys/time.h>
#define rep(x,n) for(int x = 0; x < int(n); ++x)
#define dbg(x) cerr << #x << " = " << x << endl
#define _ << " , " <<
using namespace std;

typedef long double long_double;
#define double long_double

const double eps = 1e-7;
const int inf = ~0u>>1;
const unsigned long long TLE = 60 * 60000000ULL; // minutes
const double max_alpha = 1e3;

typedef vector< double > Row;
typedef vector< Row > Matrix;
typedef set< int, less< int > > Cluster;
typedef vector< int > Prototype;
typedef pair< int, int > Pair;

namespace MersenneTwister {
  const int MT_N = 624;
  const int MT_M = 397;
  const unsigned int MT_MSB = 0x80000000U;
  const unsigned int MT_LS31B = 0x7FFFFFFFU;
  const unsigned int MT_A = 2567483615U;
  unsigned int twistory[MT_N];
  int pos;
  void build(uint seed=0) {
    twistory[0]=seed;
    for (int i=1;i<MT_N;i++) twistory[i]=1812433253U*(twistory[i-1]^(twistory[i-1]>>30))+i;
    pos=0;
  }
  void generate(void) {
    unsigned int tmp, i;
    for (i=0;i<MT_N-MT_M;i++) {
      tmp=(twistory[i]&MT_MSB)+(twistory[i+1]&MT_LS31B);
      twistory[i]=twistory[i+MT_M]^(tmp>>1)^(MT_A&-(tmp&1));
    }
    for (;i<MT_N-1;i++) {
      tmp=(twistory[i]&MT_MSB)+(twistory[i+1]&MT_LS31B);
      twistory[i]=twistory[i+MT_M-MT_N]^(tmp>>1)^(MT_A&-(tmp&1));
    }
    tmp=(twistory[i]&MT_MSB)+(twistory[0]&MT_LS31B);
    twistory[i]=twistory[MT_M-1]^(tmp>>1)^(MT_A&-(tmp&1));
  }
  unsigned int rand_unsigned() {
    if (pos==0) generate();
    uint ans=twistory[pos++];
    pos&=-(pos!=624);
    ans^=ans>>11;
    ans^=(ans<<7)&2636928640U;
    ans^=(ans<<15)&4022730752U;
    ans^=ans>>18;
    return ans;
  }
  int rand_signed() {
    return rand_unsigned()>>1;
  }
  int next_int(int n) {
    return rand_unsigned()%n;
  }
  int next_int(int a,int b) {
    return rand_unsigned()%(b-a+1)+a;
  }
};

#define get_random() MersenneTwister::rand_unsigned()

namespace Util {
  int cmp(double x, double y = 0) {
    if(fabs(x-y) <= eps) return 0; else return x < y ? -1 : +1;
  }
  double f_rand(double f_min, double f_max) {
    double f = double(get_random()) / UINT_MAX;
    return f_min + f * (f_max - f_min);
  }
  template<class T> T square(T x) {
    return x * x;
  }
  template<class T> void pv(T a, T b, ostream& out = cerr) {
    for (T i = a; i != b; ++i) out << *i << " ";
    out << endl;
  }
  template<class S, class T> S cast(T a) {
    stringstream s;
    s << a;
    S b;
    s >> b;
    return b;
  }
  void randomize(vector<int>& x) {
    for (int i=int(x.size())-1; i>0; --i) {
      swap (x[i],x[get_random()%(i+1)]);
    }
  }
};

namespace Counter {
  struct timeval start, now;
  void start_timer() {
    gettimeofday(&start, NULL);
  }
  unsigned long long elapsed() {
    gettimeofday(&now, NULL);
    return (now.tv_sec-start.tv_sec)*1000000ULL+(now.tv_usec-start.tv_usec);
  }
};

namespace Analysis {
  pair< double, double > confidence_interval(double mean, double deviation, int n, double magic = 1.959964) { // 95%
    double delta = magic * deviation / sqrt(n);
    return make_pair(mean - delta, mean + delta);
  }
  pair< double, double > build(vector< double > sample) {
    double mean = 0, deviation = 0;
    rep(i,sample.size()) mean += sample[i];
    mean /= sample.size();
    rep(i,sample.size()) deviation += Util::square(sample[i]-mean);
    deviation = sqrt(deviation/sample.size());
    return confidence_interval(mean, deviation, sample.size());
  }
};

namespace MinimumCostMaximumFlow {
  vector< int > from, to, w, cap, ant, adj;
  int n, m;
  void add(int x, int y, int z, int c) {
    while(x+1 > n) n++, adj.push_back(-1); while(y+1 > n) n++, adj.push_back(-1);
    from.push_back(x), to.push_back(y), w.push_back(z), cap.push_back(c), ant.push_back(adj[x]), adj[x] = m++;
    swap(x,y), z = -z, c = 0;
    from.push_back(x), to.push_back(y), w.push_back(z), cap.push_back(c), ant.push_back(adj[x]), adj[x] = m++;
  }
  Pair mcmf(int source, int sink) {
    Pair p(0,0); // cost,flow
    for(;;) {
      vector< int > pot(n,inf), pai(n,-1);
      pot[source] = 0;
      rep(k,n) rep(i,n) for(int j = adj[i]; j >= 0; j = ant[j]) {
        if(cap[j] && pot[i] + w[j] < pot[to[j]]) {
          pot[to[j]] = pot[i] + w[j];
          pai[to[j]] = j;
        }
      }
      if(pot[sink] < inf) {
        int cost = 0, flow = inf;
        for(int x = sink; x != source; x = from[pai[x]]) if(cap[pai[x]] < flow) flow = cap[pai[x]];
        for(int x = sink; x != source; x = from[pai[x]]) cap[pai[x]] -= flow, cap[pai[x]^1] += flow, cost += w[pai[x]]*flow;
        p.first += cost, p.second += flow;
      }
      else break;
    }
    return p;
  }
  void init() {
    n = m = 0;
    adj.clear(), from.clear(), to.clear(), w.clear(), cap.clear(), ant.clear();
  }
};

namespace Validation {
  pair< double, vector< int> > global_error(int n, vector< vector<int> >& tab) {
    int N = tab.size() - 1, M = tab[0].size() - 1, source = 0, sink = N + M + 1;
    vector< int > match(N, -1);
    MinimumCostMaximumFlow::init();
    rep(i,N) rep(j,M) MinimumCostMaximumFlow::add(i + 1, j + N + 1, -tab[i][j], 1);
    rep(i,N) MinimumCostMaximumFlow::add(source, i + 1, 0, 1);
    int cflow = (N == M ? 1 : N);
    rep(i,M) MinimumCostMaximumFlow::add(i + N + 1, sink, 0, cflow);
    Pair ret = MinimumCostMaximumFlow::mcmf(source, sink);
    rep(i,N) rep(j,M) if(MinimumCostMaximumFlow::cap[(i*M+j)<<1] != cflow) match[i] = j;
    assert(ret.second == N);
    return make_pair(1.0 + double(ret.first) / n, match);
  }
  int comb(int x) {
    return x * (x - 1) / 2;
  }	
  double adjusted_rand_index(int n, vector< vector<int> >& tab) {
    double term[4] = {0}, temp[2] = {0};
    double pot = 1.0 / comb(n);
    int N = tab.size() - 1, M = tab[0].size() - 1;
    rep(i,N) rep(j,M) term[1] += comb(tab[i][j]);
    rep(i,N) temp[0] += comb(tab[i][M]);
    rep(i,M) temp[1] += comb(tab[N][i]);
    term[2] = pot * (temp[0] * temp[1]);
    term[3] = 0.5 * (temp[0] + temp[1]);
    return (term[1]-term[2])/(term[3]-term[2]);
  }
  double f_measure(int n, vector< vector<int> >& tab) {
    double F = 0;
    int N = tab.size() - 1, M = tab[0].size() - 1;
    rep(j,M) {
      double vmax = 0;
      rep(i,N) if(tab[i][j] != 0) {
        double rappel = double(tab[i][j]) / tab[N][j];
        double precision = double(tab[i][j]) / tab[i][M];
        vmax = max(vmax, 2 * rappel * precision / (rappel + precision));
      }
      F += vmax * tab[N][j];
    }
    return F / n;
  }
};

struct Gauss {
  vector< vector< double > > A;
  vector< double > B;
  int N, M;
  Gauss(int N) {
    this->N = N;
    this->M = 0;
    A.clear();
    B.clear();
  }
  void create_row() {
    A.push_back(vector< double >(N,0.0));
    B.push_back(0.0);
    M++;
  }
  int solve() {
    int row = 0;
    rep(pivot,N) {
      int pos = -1;
      double vmax = 0;
      for(int i = row; i < M; i++) {
        if(Util::cmp(A[i][pivot]) && Util::cmp(fabs(A[i][pivot]),vmax) > 0) {
          pos = i;
          vmax = fabs(A[i][pivot]);
        }
      }
      if(pos == -1) continue;
      if(pos != row) swap(A[row],A[pos]), swap(B[row],B[pos]);
      double F = A[row][pivot];
      for(int i = row; i < N; i++) A[row][i] /= F;
      B[row] /= F;
      for(int i = row + 1; i < M; i++) {
        F = A[i][pivot];
        if(Util::cmp(F) == 0) continue;
        for(int j = pivot; j < N; j++) A[i][j] = (A[row][j] * F) - (A[i][j] * A[row][pivot]);
        B[i] = (B[row] * F) - (B[i] * A[row][pivot]);
      }
      row++;
    }
    for(int i = row; i < M; i++) if(Util::cmp(B[i])) return 0;
    if(row < N) return 0;
    for(int i = row-1; i >= 0; i--) {
      if(Util::cmp(A[i][i]) == 0) return 0;
      B[i] /= (A[i][i]);
      rep(j,i) B[j] -= (A[j][i] * B[i]);
    }
    return 1;
  }
  double& operator[](size_t __n) {
    return __n == 0 ? B[M - 1] : A[M - 1][__n - 1];
  }
  double& operator()(size_t __n) {
    return B[__n - 1];
  }
};

struct Latex {
  ostream& out;
  Latex(ostream& out = cout): out(out) {	
    assert(out);
  }
  void begin() {
    out<<"\\documentclass[en,msc,oneside,onehalfspacing]{risethesis}\n"
      "\\usepackage{natbib}\n"
      "\\usepackage{babel}\n"
      "\\usepackage{supertabular}\n"
      "\\usepackage{longtable}\n"
      "\\usepackage{fancybox}\n"
      "\\usepackage{acronym}\n"
      "\\usepackage{moreverb}\n"
      "\\usepackage{multirow}\n"
      "\\usepackage{subfig}\n"
      "\\usepackage{enumitem}\n"
      "\\setlist{nolistsep}\n"
      "\\usepackage[linkcolor=black,citecolor=black,urlcolor=black,colorlinks,pdfpagelabels,pdftitle={Semi-Supervised Clustering Algorithm Analysis},pdfauthor={Filipe Martins de Melo}]{hyperref}\n"
      "\\address{Recife}\n"
      "\\universitypt{Universidade Federal de Pernambuco}\n"
      "\\universityen{Federal University of Pernambuco}\n"
      "\\departmentpt{Centro de Informática}\n"
      "\\departmenten{Center for Informatics}\n"
      "\\programpt{Mestrado em Ciência da Computação}\n"
      "\\programen{Master Degree in Computer Science}\n"
      "\\majorfieldpt{Ciência da Computação}\n"
      "\\majorfielden{Computer Science}\n"
      "\\title{Semi-Supervised Clustering Algorithm Analysis}\n"
      "\\date{\\today}\n"
      "\\author{Filipe Martins de Melo}\n"
      "\\adviser{Prof. Dr. Francisco de Assis Tenório Carvalho}\n"
      "\\begin{document}\n";
    out<<endl;
  }
  void end() {
    out<<"\\end{document}";
    out<<endl;
    out.flush();
  }
  void print_name(string name) {
    out<<"%% here it goes %%\n"
      "\\centering\n"
      "\\paragraph{\\textit{Fuzzy Analysis @ "<<name<<"}}\n";
  }
  void print_description(const map< string, vector< string > > &M) {
    out<<"\\small\n"
      "\\begin{description}\n";
    for(map< string, vector< string > >::const_iterator i = M.begin(); i != M.end(); i++) {
      string x = i->first;
      rep(j,x.size()) if(x[j] == '_') x[j] = ' ';
      out<<"\\item "<<x;
      vector< string > y = i->second;
      rep(j,y.size()) {
        x = string("");
        rep(k,y[j].size()) if(y[j][k] == '_') x += "\\_"; else x += y[j][k];
        out<<x<<"\n";
      }
    }
    out<<"\\end{description}\n";
  }
  void print_itemize(vector< string > elem) {
    out<<"\\begin{itemize}\n";
    rep(i,elem.size()){
      out<<"\\item "<<elem[i]<<"\n";
    }
    out<<"\\end{itemize}\n";
  }
  void print_itemize(vector< pair< string, vector< string > >  > elem) {
    out<<"\\begin{itemize}\n";
    rep(i,elem.size()){
      out<<"\\item "<<elem[i].first<<"\n";
      print_itemize(elem[i].second);
    }
    out<<"\\end{itemize}\n";
  }
  void print_table(int n, string caption, vector< string > elem) {
    out<<"\\begin{table}[h]\n"
      "\\caption{"<<caption<<"}\n"<<
      "\\centering\n"
      "\\begin{tabular}{c";rep(i,n-1) out<<" c"; out<<"}\n"
      "\\hline\n";
    rep(i,elem.size()) {
      if(i == n) out<<"\\hline\n";
      out<<elem[i]<<((i+1)%n?"&":"\\\\\n");
    }
    out<<"\\hline\n"
      "\\end{tabular}\n"
      "\\end{table}\n";
  }
  void print_long_table(int n, string caption, vector< string > elem) {
    out<<"\\begin{longtable}{c";rep(i,n-1) out<<" c"; out<<"}\n";
    out<<"\\caption{"<<caption<<"}\n"
      "\\centering\n"
      "\\tabularnewline\n\\hline\n";
    rep(i,elem.size()) {
      if(i == n) out<<"\\hline\n";
      out<<elem[i]<<((i+1)%n?"&":"\\\\\n");
    }
    out<<"\\hline\n"
      "\\end{longtable}\n";
  }
  Latex& operator<<(string str) {
    out << str;
    return *this;
  }
  Latex& operator<<(ostream& out2) {
    out << out2;
    return *this;
  }
};

// data
unsigned long long seed;
map< string, vector< string> > M;
int class_variable;
int initialization_number;
int maximum_iteration_number;
int N, C, P, T;
double alpha, label;
vector< Cluster > priori_cluster;
set< Pair > must_link, must_not_link;
vector< string > file;
vector< Matrix > table;
bool automatic_alpha;

struct Scanner {  
  Scanner(string config) {
    M.clear();
    ifstream in(config.c_str(), ios::in);
    for(string kind, value; in>>value;){
      if(value[0] == '[' && value[value.size() - 1] == ']') kind = value;
      else M[kind].push_back(value);
    }
    in.close();
    class_variable = Util::cast<int>(M["[class_variable]"][0]);
    initialization_number = Util::cast<int>(M["[initialization_number]"][0]);
    maximum_iteration_number = Util::cast<int>(M["[maximum_iteration_number]"][0]);
    C = Util::cast<int>(M["[cluster_number]"][0]);
    P = Util::cast<int>(M["[prototype_number]"][0]);
    alpha = Util::cast<double>(M["[alpha]"][0]);
    label = Util::cast<double>(M["[label]"][0]);
    file = M["[input]"];
    T = file.size();
    read_individual_number_and_priori_cluster(file[0]);
    M["[individual_number]"].push_back(Util::cast<string>(N));
    table = vector< Matrix >(T, Matrix(N, Row(N,0.0)));
    rep(i,T) read_table(file[i], table[i]);
    automatic_alpha = M.count("[automatic_alpha]") && M["[automatic_alpha]"][0] == string("True");
  }
  void read_individual_number_and_priori_cluster(string file_name) {
    ifstream in(file_name.c_str(), ios::in);
    assert(in);
    string line;
    do if(!getline(in,line)) assert(0); while(line.find("indiv_nb") == string::npos);
    assert(sscanf(line.c_str()," indiv_nb = %d", &N) == 1);
    do if(!getline(in,line)) assert(0); while(line.find("RECTANGLE_MATRIX = (") == string::npos);
    priori_cluster.clear();
    rep(i,N) {
      line.clear();
      char c;
      int ct = 0;
      while(ct == 0) {
        assert(in>>c);
        if(c == '(') ct++;
      }
      while(ct != 0) {
        assert(in>>c);
        if(c == '(') ct++; else if(c == ')') ct--;
        if(ct != 0) line.push_back(c);
      }
      line.push_back(',');
      ct = 0;
      rep(j,line.size()) if(line[j] != ' ') line[ct++] = line[j];
      line.resize(ct);
      string piece;
      vector< string > v;
      ct = 0;
      rep(j,line.size()) {
        if(line[j] == '(') ct++; else if(line[j] == ')') ct--;
        piece.push_back(line[j]);
        if(ct == 0) {
          if(piece[0] == ',') piece.clear();
          else if(line[j+1] == ',') v.push_back(piece);
        }
      }
      piece = v[class_variable - 1];
      unsigned int priori;
      assert(sscanf(piece.c_str(),"%u",&priori) == 1);
      while(priori > priori_cluster.size()) priori_cluster.push_back(Cluster());
      priori_cluster[priori - 1].insert(i);
    }
    in.close();
  }
  void read_table(string file_name, Matrix& M) {
    ifstream in(file_name.c_str(), ios::in);
    assert(in);
    string line;
    do if(!getline(in,line)) assert(0); while(line.find("DIST_MATRIX= (") == string::npos);
    rep(i,N) {
      line.clear();
      char c;
      for(c = '?';c != '('; assert(in>>c));
      for(assert(in>>c);c != ')'; assert(in>>c)) line.push_back(c);
      rep(j,line.size()) if(line[j] == '(' || line[j] == ')' || line[j] == ',') line[j] = ' ';
      stringstream s(line);
      rep(j,i+1) {
        double d;
        assert(s>>d);
        M[i][j] = M[j][i] = d;
      }
    }
    in.close();
  }
  void generate_restrictions() {
    must_link.clear();
    must_not_link.clear();
    map< int, int > mapping;
    rep(k,priori_cluster.size()) for(Cluster::const_iterator x = priori_cluster[k].begin(); x != priori_cluster[k].end(); x++) mapping[*x] = k;
    int n = (N * label + 0.7777777);
    assert(n <= N);
    vector< int > object(N);
    rep(i,N) object[i] = i;
    Util::randomize(object);
    rep(i,n) rep(j,n) if(i != j) {
      int a = object[i], b = object[j];
      if(mapping[a] == mapping[b]) must_link.insert( make_pair(a,b) );
      else must_not_link.insert( make_pair(a,b) );
    }
  }
  void manual(vector< int > object) {
    must_link.clear();
    must_not_link.clear();
    map< int, int > mapping;
    rep(k,priori_cluster.size()) for(Cluster::const_iterator x = priori_cluster[k].begin(); x != priori_cluster[k].end(); x++) mapping[*x] = k;
    int n = (N * label + 0.7777777);
    assert(n <= N);
    rep(i,n) rep(j,n) if(i != j) {
      int a = object[i], b = object[j];
      if(mapping[a] == mapping[b]) {
        must_link.insert( make_pair(a,b) );
      }
      else {
        must_not_link.insert( make_pair(a,b) );
      }
    }
  }
};

namespace Algorithm {
  bool update_coefficient_table = false;
  bool update_coefficient_cluster_table = false;
  struct Answer {
    int id;
    double J, restriction;
    vector< double > objective;
    vector< Cluster > cluster;
    vector< Prototype > prototype;
    Matrix U, coefficient;
    Answer(int id) {
      this->id = id;
      this->J = DBL_MAX;
      this->cluster = vector< Cluster >(C);
      this->prototype = vector< Prototype >(C);
      this->U = Matrix(N, Row(C));
      this->coefficient = Matrix(C, Row(T));
    }
    vector< int > get_order() {
      vector< double > doubt(N, 0.0);
      rep(i,N) rep(k,C) doubt[i] = max(doubt[i],U[i][k]);
      vector< pair< double, int > > v(N);
      rep(i,N) v[i] = make_pair(doubt[i],i);
      sort(v.begin(),v.end());
      vector< int > order(N);
      rep(i,N) order[i] = v[i].second;
      return order;
    }
    void srand() {
      if(C * P > N) __throw_invalid_argument("Invalid number of prototypes.");
      vector< int > new_prototype(N);
      rep(i,N) new_prototype[i] = i;
      Util::randomize(new_prototype);
      rep(k,C) rep(p,P) prototype[k].push_back(new_prototype[k*P+p]);
      rep(k,C) rep(t,T) coefficient[k][t] = 1.0;
      J = restriction = DBL_MAX;
      objective.clear();
      rep(i,N) {
        vector< double > dist(C,0.0);
        rep(t,T) rep(k,C) rep(p,P) dist[k] += table[t][i][prototype[k][p]] * coefficient[k][t];
        vector< int > is_zero;
        rep(k,C) if(Util::cmp(dist[k]) <= 0) is_zero.push_back(k);
        if(is_zero.size()) {
          rep(k,C) U[i][k] = 0.0;
          rep(k,is_zero.size()) U[i][is_zero[k]] = 1.0 / is_zero.size();
          continue;
        }
        rep(k,C) {
          U[i][k] = 0;
          rep(h,C) U[i][k] += dist[k] / dist[h];
          assert(Util::cmp(U[i][k]));
          U[i][k] = 1.0 / U[i][k];
        }
      }
      updateCluster();
      updateJ();
    }
    void updateU() {
      for(int turn=0, max_turn=10;turn<max_turn;turn++){
        double diff = 0;
        vector< int > order(N);
        rep(i,N) order[i] = i;
        Util::randomize(order);
        rep(_i,N) {
          int i = order[_i];
          vector< double > previous = U[i], a(C,0.0), b(C,0.0);
          rep(k,C) rep(t,T) rep(p,P) a[k] += table[t][i][prototype[k][p]] * coefficient[k][t];
          rep(k,C) a[k] *= 2.0;
          vector< int > is_zero;
          rep(k,C) if(Util::cmp(a[k]) <= 0) is_zero.push_back(k);
          if(is_zero.size()) {
            rep(k,C) U[i][k] = 0.0;
            rep(k,is_zero.size()) U[i][is_zero[k]] = 1.0 / is_zero.size();
            goto end;
          }
          for(set< Pair >::const_iterator x = must_link.lower_bound(make_pair(i,-1)); x != must_link.end() && x->first == i; x++) { // must link
            int m = x->second;
            rep(r,C) rep(s,C) if(r != s) b[r] += U[m][s];
          }
          for(set< Pair >::const_iterator x = must_not_link.lower_bound(make_pair(i,-1)); x != must_not_link.end() && x->first == i; x++) { // must not link
            int m = x->second;
            rep(r,C) b[r] += U[m][r];
          }
          rep(k,C) b[k] *= alpha;
          is_zero = vector< int >(C,0);
          for(vector< double > psi(C,0.0);;) {
            double gamma = ({
                double num = 0, den = 0;
                rep(k,C) num += (b[k] - psi[k]) / a[k], den += 1.0 / a[k];
                (1.0 + num) / den;
                });
            bool ok = true;
            rep(k,C) {
              U[i][k] = (gamma + psi[k] - b[k]) / a[k];
              if(is_zero[k]) {
                assert(Util::cmp(U[i][k]) == 0);
                U[i][k] = 0.0;
              }
              else if(Util::cmp(U[i][k]) < 0) {
                is_zero[k] = 1;
                ok = false;
              }
            }
            if(ok) break;
            vector< int > mapping(C,0);
            rep(k,C) mapping[k] = is_zero[k] + (k ? mapping[k - 1] : 0);
            Gauss gauss(mapping[C-1]);
            double sum_a = 0, sum_b = 0;
            rep(k,C) sum_a += 1.0 / a[k], sum_b += b[k] / a[k];
            rep(k,C) if(is_zero[k]) {
              gauss.create_row();
              gauss[0] = 1.0 + sum_b - b[k] * sum_a;
              gauss[mapping[k]] = -sum_a;
              rep(h,C) if(is_zero[h]) gauss[mapping[h]] += 1.0 / a[h];
              for(int h = k + 1; h < C; h++) if(is_zero[h]) {
                gauss.create_row();
                gauss[mapping[k]] = +1.0;
                gauss[mapping[h]] = -1.0;
                gauss[0] = b[k] - b[h];
              }
            }
            assert(gauss.solve());
            rep(k,C) if(is_zero[k]) psi[k] = gauss(mapping[k]);
          }
          if(true) {
            rep(k,C) if(U[i][k] < 0) U[i][k] = 0; else if(U[i][k] > 1.0) U[i][k] = 1.0;
            double sum = 0;
            rep(k,C) sum += U[i][k];
            assert(fabs(sum-1) <= eps);
            rep(k,C) U[i][k] /= sum;
          }
end: rep(k,C) diff = max(diff, previous[k] - U[i][k]);
        }
        if(diff <= eps) break;
      }
    }
    void updatePrototypes() {
      rep(k,C) {
        vector< double > dist(N,0.0);
        rep(t,T) rep(i,N) rep(j,N) dist[i] += table[t][i][j] * Util::square(U[j][k]) * coefficient[k][t];
        vector< pair< double, int > > v(N);
        rep(i,N) v[i] = make_pair(dist[i], i);
        nth_element(v.begin(),v.begin()+P,v.end());
        rep(p,P) prototype[k][p] = v[p].second;
      }
    }
    void updateCluster() {
      rep(k,C) cluster[k].clear();
      rep(i,N) {
        vector< pair< double, int > > v(C);
        rep(k,C) v[k] = make_pair(U[i][k],k);
        cluster[max_element(v.begin(),v.end())->second].insert(i);
      }
    }
    void updateCoefficient() {
      if(update_coefficient_table) {
        double num = 1.0;
        vector< double > den(T,0.0);
        rep(t,T) {
          rep(k,C) rep(i,N) rep(p,P) den[t] += table[t][i][prototype[k][p]] * Util::square(U[i][k]);
          assert(den[t]);
          num *= pow(den[t], 1.0 / T);
        }
        rep(k,C) rep(t,T) coefficient[k][t] = num / den[t];
      }
      else if(update_coefficient_cluster_table) rep(k,C) {
        double num = 1.0;
        vector< double > den(T,0.0);
        rep(t,T) rep(i,N) rep(p,P) den[t] += table[t][i][prototype[k][p]] * Util::square(U[i][k]);
        rep(t,T) num *= pow(den[t], 1.0 / T);
        if(Util::cmp(num) == 0) assert(0);
        rep(t,T) assert(den[t]);
        rep(t,T) coefficient[k][t] = num / den[t];
      }
    }
    double updateJ() {
      double J1 = 0;
      Matrix d(N,Row(C,0.0));
      rep(t,T) rep(i,N) rep(k,C) rep(p,P) d[i][k] += table[t][i][prototype[k][p]] * coefficient[k][t];
      rep(i,N) rep(k,C) J1 += Util::square(U[i][k]) * d[i][k];
      double J2 = 0;
      for(set< Pair >::const_iterator x = must_link.begin(); x != must_link.end(); x++) { // must link
        int l = x->first, m = x->second;
        rep(r,C) rep(s,C) if(r != s) J2 += U[l][r] * U[m][s];
      }
      for(set< Pair >::const_iterator x = must_not_link.begin(); x != must_not_link.end(); x++) { // must not link
        int l = x->first, m = x->second;
        rep(r,C) J2 += U[l][r] * U[m][r];
      }
      double new_J = J1 + alpha * J2;
      objective.push_back(J = new_J);
      restriction = J2;
      return J;
    }
    vector< vector< int > > confusing_matrix() {
      int k = C;
      int p = priori_cluster.size();
      vector< vector< int > > table(k+1,vector< int >(p+1,0));
      rep(i,k) rep(j,p) for(Cluster::const_iterator x = cluster[i].begin(); x != cluster[i].end(); x++) table[i][j] += priori_cluster[j].count(*x);
      rep(i,k) table[i][p] = cluster[i].size();
      rep(i,p) table[k][i] = priori_cluster[i].size();
      table[k][p] = N;
      return table;
    }
    bool operator<(const Answer& answer) const {
      return J < answer.J;
    }
  };
  double find_alpha(double start = max_alpha) {
    const double this_alpha = alpha;
    double ans = start, take = start, R= log2(start) + 10;
    rep(i,R) {
      double nans = ans - take; take /= 2;
      if(nans < 0) continue; else alpha = nans;
      int ok = 0;
      rep(j,5) {
        Answer x(-1);
        x.srand();
        rep(r,350) {
          x.updatePrototypes();
          x.updateCoefficient();
          x.updateU();
          x.updateCluster();
          double old_J = x.J;
          double new_J = x.updateJ();
          if(fabs(old_J-new_J) <= eps || x.restriction <= eps) break;
        }
        if(x.restriction <= eps) {
          ok = 1;
          break;
        }
      }
      if(ok) ans = nans;
    }
    alpha = this_alpha;
    return ans;
  }
  Answer main(Latex& latex, bool print = true) {
    const double this_alpha = alpha;
    if(automatic_alpha) {
      double alpha1 = find_alpha();
      alpha = alpha1;
#ifdef DETAILED
      dbg(alpha1);
#endif
    }
    Counter::start_timer();
    Answer opt(-1);
    rep(initialization,initialization_number) {
      if(Counter::elapsed()>TLE) {
#ifdef DETAILED
        cerr << "tle at " << initialization << endl;
#endif
        break;
      }
      Answer now(initialization);
      now.srand();
      rep(iteration,maximum_iteration_number) {
#ifdef EXTREMELY_DETAILED
        dbg(iteration _ now.J _ now.restriction);
#endif
        now.updatePrototypes();
        now.updateCoefficient();
        now.updateU();
        now.updateCluster();
        double old_J = now.J;
        double new_J = now.updateJ();
#ifdef DETAILED
        if(new_J>old_J&&fabs(new_J-old_J)>1e-0) {
          fprintf(stderr,"increment (%+.7Lf)\n", new_J - old_J);
        }
#endif
        if(fabs(old_J - new_J) <= eps) break; else opt = min(opt,now);
      }
#ifdef EXTREMELY_DETAILED
      dbg(initialization _ now.J _ now.restriction);
#endif
    }
    if(true) { // optmize
      const double magic = M.count("[magic]") ? Util::cast<double>(M["[magic]"][0]) : 1e9;
      const double inside = 1.0 - 1.0 / (C + magic), outside = 1.0 / (C + magic);
      __typeof(must_link) ml = must_link;
      __typeof(must_not_link) mnl = must_not_link;
      rep(k,C) rep(i,N) {
        if(opt.U[i][k] > inside) rep(j,i) {
          if(opt.U[j][k] > inside) {
            must_link.insert( make_pair(i,j) );
            must_link.insert( make_pair(j,i) );
          }
          else if(opt.U[j][k] < outside) {
            must_not_link.insert( make_pair(i,j) );
            must_not_link.insert( make_pair(j,i) );
          }
        }
        else if(opt.U[i][k] < outside) rep(j,i) {
          if(opt.U[j][k] > inside) {
            must_not_link.insert( make_pair(i,j) );
            must_not_link.insert( make_pair(j,i) );
          }
        }
      }
      opt.J = DBL_MAX;
      if(automatic_alpha) {
        double alpha2 = find_alpha();
        alpha = alpha2;
#ifdef DETAILED
        dbg(alpha2);
#endif
      }
      Counter::start_timer();
      rep(initialization2,initialization_number) {
        if(Counter::elapsed()>TLE) {
#ifdef DETAILED
          cerr << "tle tle at " << initialization2 << endl;
#endif
          break;
        }
        Answer now2(initialization2);
        now2.srand();
        rep(iteration2,maximum_iteration_number) {
#ifdef EXTREMELY_DETAILED
          dbg(iteration2 _ now2.J _ now2.restriction);
#endif
          now2.updatePrototypes();
          now2.updateCoefficient();
          now2.updateU();
          now2.updateCluster();
          double old_J2 = now2.J;
          double new_J2 = now2.updateJ();
#ifdef DETAILED
          if(new_J2>old_J2&&fabs(new_J2-old_J2)>1e-0) {
            fprintf(stderr,"increment2 (%+.7Lf)\n", new_J2 - old_J2);
          }
#endif
          if(fabs(old_J2 - new_J2) <= eps) break; else opt = min(opt,now2);
        }
#ifdef EXTREMELY_DETAILED
        dbg(initialization2 _ now2.J _ now2.restriction);
#endif
      }
      must_link = ml;
      must_not_link = mnl;   
    }
    if(print) {
      if(automatic_alpha) {
        M["[alpha]"] = vector< string >(1,Util::cast<string>(alpha));
      }
      char text[1<<10];
      string caption;
      vector<string> elem;
      latex.print_name("X algorithm"); // name
      latex.print_description(M);
      caption="Initialization \\textbf{"+Util::cast<string>(opt.id+1)+"}"; // best initialization
      elem.clear();
      elem.push_back("Iteration");
      elem.push_back("Objective function");
      elem.push_back("~");
      rep(i,opt.objective.size()) {
        elem.push_back(Util::cast<string>(i+1));
        sprintf(text,i+1==(int)opt.objective.size()?"\\textbf{%.7Lf}":"%.7Lf",opt.objective[i]);
        elem.push_back(text);
        sprintf(text,i == 0 ? "~":"(%.7Lf)",opt.objective[i]-opt.objective[i?i-1:0]);
        elem.push_back(text);
      }
      latex.print_long_table(3,caption,elem);
      vector< vector< int > > table = opt.confusing_matrix();
      pair< double, vector< int > > mcmf = Validation::global_error(N,table);
      caption="Confusing matrix";
      elem.clear();
      elem.push_back("~");
      rep(k,C) elem.push_back(Util::cast<string>(k+1));
      rep(k,C) {
        elem.push_back(Util::cast<string>(k+1)+"'");
        rep(h,C) {
          string s = Util::cast<string>(table[k][h]);
          if(mcmf.second[k] == h) s = "\\textbf{" + s + "}";
          elem.push_back(s);
        }
      }
      latex.print_table(C+1,caption,elem);
      elem.clear();
      sprintf(text,"Global error: \\textbf{%.7Lf\\%%}", mcmf.first * 100.0);
      elem.push_back(text);
      latex.print_itemize(elem);
      elem.clear();
      sprintf(text,"Adjusted rand index: \\textbf{%.7Lf}",Validation::adjusted_rand_index(N,table));
      elem.push_back(text);
      sprintf(text,"F measure: \\textbf{%.7Lf}",Validation::f_measure(N,table));
      elem.push_back(text);
      latex.print_itemize(elem);
      if(update_coefficient_table) {
        elem.clear();
        rep(t,T) {
          sprintf(text,"Table %d: \\textbf{%.7Lf}",t+1,opt.coefficient[0][t]);
          elem.push_back(text);
        }
        latex.print_itemize(elem);
      }
      else if(update_coefficient_cluster_table) {
        vector< pair< string, vector< string > > > v_elem(C);
        rep(k,C) {
          sprintf(text,"Cluster %d:", k+1);
          pair< string, vector< string > >& p_elem = v_elem[k];
          p_elem.first = text;
          rep(t,T) {
            sprintf(text,"Table %d: \\textbf{%.7Lf}", t+1, opt.coefficient[k][t]);
            p_elem.second.push_back(text);
          }
        }
        latex.print_itemize(v_elem);
      }
      if(true) {
        elem.clear();
        rep(k,C) {
          stringstream ss;
          rep(p,P) {
            if(p) ss << ",";
            ss << opt.prototype[k][p]+1;
          }
          sprintf(text,"Cluster %d: \\textbf{%s}",k+1,ss.str().c_str());
          elem.push_back(text);
        }
        latex.print_itemize(elem);
      }
      caption = "Fuzzy clustering";
      elem.clear();
      elem.push_back("~");
      rep(k,C) elem.push_back(Util::cast<string>(k+1));
      rep(i,N) {
        sprintf(text,"%.4d",i+1);
        elem.push_back(text);
        rep(k,C) {
          sprintf(text,opt.cluster[k].count(i)?"\\textbf{%.7Lf}":"%.7Lf",opt.U[i][k]);
          elem.push_back(text);
        }
      }
      latex.print_long_table(C+1,caption,elem);
      if(true) {
        latex<<"\\tiny\n";
        stringstream s;
        s << setprecision(6) << fixed << Counter::elapsed() / 60000000.0 << " minutes, with seed = " << seed << ".\n";
        latex<<s.str();
      }
      latex<<"\\pagebreak\n";
    }
    alpha = this_alpha;
    return opt;
  }
};

namespace Another {
  const int q = 2;
  struct Answer {
    int id;
    double J, beta, restriction;
    vector< double > objective;
    vector< Cluster > cluster;
    Matrix U, W;
    Answer(int id) {
      this->id = id;
      this->J = DBL_MAX;
      this->cluster = vector< Cluster >(C);
      this->U = Matrix(C, Row(N));
      this->W = Matrix(C, Row(T));
    }
    vector< int > get_order() {
      vector< double > doubt(N, 0.0);
      rep(i,N) rep(k,C) doubt[i] = max(doubt[i],U[k][i]);
      vector< pair< double, int > > v(N);
      rep(i,N) v[i] = make_pair(doubt[i],i);
      sort(v.begin(),v.end());
      vector< int > order(N);
      rep(i,N) order[i] = v[i].second;
      return order;
    }
    void srand() {
      rep(i,N) {
        double sum = 0;
        rep(k,C) U[k][i] = get_random(), sum += U[k][i];
        assert(Util::cmp(sum) != 0);
        rep(k,C) U[k][i] /= sum;
      }
      updateCluster();
      beta = 0;
      J = restriction = DBL_MAX;
      objective.clear();
      rep(k,C) rep(t,T) W[k][t] = 1.0 / T;
      updateJ();
    }
    void updateCluster() {
      rep(k,C) cluster[k].clear();
      rep(i,N) {
        vector< pair< double, int > > v(C);
        rep(k,C) v[k] = make_pair(U[k][i],k);
        cluster[max_element(v.begin(),v.end())->second].insert(i);
      }
    }
    bool update() {
      Matrix Upow2(C,Row(N,0)), newU(C,Row(N,0)), newW(C,Row(T,0));
      rep(k,C) rep(i,N) Upow2[k][i] = Util::square(U[k][i]);
      // equation [15]
      vector< Matrix > R(C,Matrix(N,Row(N,0)));
      rep(i,N) rep(j,N) rep(k,C) rep(t,T) R[k][i][j] += pow(W[k][t],q) * table[t][i][j];
      // equation [8]
      vector< Matrix > Rbeta(C,Matrix(N,Row(N,0)));
      rep(k,C) rep(i,N) rep(j,N) Rbeta[k][i][j] += R[k][i][j] + beta * ((i!=j) - (i==j));
      // equation [6]
      vector< Row > v(C,Row(N,0));
      rep(k,C) {
        double den = 0;
        rep(i,N) v[k][i] = Upow2[k][i], den += v[k][i];
        if(Util::cmp(den) <= 0) return false;
        rep(i,N) v[k][i] /= den;
      }
      // equation [5]
      Matrix d(C,Row(N,0));
      rep(k,C) {
        Row A(N,0);
        rep(i,N) rep(j,N) A[i] += Rbeta[k][i][j] * v[k][j];
        double B=0;
        rep(i,N) B += v[k][i] * A[i] / 2.0;
        rep(i,N) d[k][i] = A[i] - B;
      }
      // equation [9]
      double delta=0;
      rep(k,C) rep(i,N) if(Util::cmp(d[k][i]) < 0) {
        double norm = 0;
        rep(j,N) norm += Util::square(v[k][j] - (i==j));
        if(Util::cmp(norm) <= 0) return false;
        delta = max(delta,-2 * d[k][i] / norm);
      }
      beta += delta;
      delta /= 2;
      if(Util::cmp(delta) > 0) rep(k,C) rep(i,N) {
        double norm = 0;
        rep(j,N) norm += Util::square(v[k][j] - (i==j));
        d[k][i] += delta * norm;
      }
      // u-rfcm
      rep(k,C) rep(i,N) assert(Util::cmp(d[k][i]) >= 0);
      rep(k,C) rep(i,N) {
        if(Util::cmp(d[k][i]) > 0) {
          double val = 0;
          rep(h,C) if(Util::cmp(d[h][i]) > 0) val += d[k][i] / d[h][i];
          assert(Util::cmp(val) != 0);
          newU[k][i] = 1.0 / val;
        }
        else {
          newU[k][i] = 0.0;
        }
      }
      // u-const
      Matrix a = d, c(C,Row(N,0));
      rep(k,C) rep(i,N) a[k][i] *= 2;
      rep(k,C) rep(i,N) {
        for(set< Pair >::const_iterator x = must_link.lower_bound(make_pair(i,-1)); x != must_link.end() && x->first == i; x++) { // must link
          int m = x->second;
          rep(r,C) if(r != k) c[k][i] += U[r][m];
        }
        for(set< Pair >::const_iterator x = must_not_link.lower_bound(make_pair(i,-1)); x != must_not_link.end() && x->first == i; x++) { // must not link
          int m = x->second;
          c[k][i] += U[k][m];
        }
      }
      vector< double > ct(N,0.0);
      rep(i,N) {
        double num = 0, den = 0;
        rep(k,C) {
          if(Util::cmp(a[k][i]) == 0) return false;
          num += c[k][i] / a[k][i];
          den += 1.0 / a[k][i];
        }
        assert(Util::cmp(den) != 0);
        ct[i] = num / den;
      }
      rep(k,C) rep(i,N) newU[k][i] += alpha * (ct[i] - c[k][i]) / a[k][i];
      bool normalize = false;
      rep(k,C) rep(i,N) if(newU[k][i] < 0) newU[k][i] = 0, normalize = true; else if(newU[k][i] > 1.0) newU[k][i] = 1.0, normalize = true;
      rep(k,C) rep(i,N) newU[k][i] = max(newU[k][i],eps);
      if(normalize) {
        rep(i,N) {
          double sum = 0;
          rep(k,C) sum += newU[k][i];
          assert(Util::cmp(sum) != 0);
          rep(k,C) newU[k][i] /= sum;
        }
      }
      rep(k,C) rep(i,N) Upow2[k][i] = Util::square(newU[k][i]);
      if(true){
        // equation [22]
        Matrix D(C,Row(T,0));
        rep(k,C) rep(t,T) rep(i,N) rep(j,N) D[k][t] += Upow2[k][i] * Upow2[k][j] * table[t][i][j];
        rep(k,C) rep(t,T) {
          double val = 0;
          rep(p,T) {
            if(Util::cmp(D[k][p]) == 0) return false;
            val += pow(D[k][t] / D[k][p], 1.0 / (q - 1));
          }
          assert(Util::cmp(val) != 0);
          newW[k][t] = 1.0 / val;
        }
      }
      // check
      if(!updateJ(newU,newW)) return false;
      rep(k,C) rep(i,N) U[k][i] = newU[k][i];
      rep(k,C) rep(t,T) W[k][t] = newW[k][t];
      updateCluster();
      return true;
    }
    bool updateJ() {
      return updateJ(U, W);
    }
    bool updateJ(const Matrix &u, const Matrix &w) {
      double J1 = 0;
      rep(k,C) {
        double num = 0;
        rep(i,N) rep(j,N) {
          double val = 0;
          rep(t,T) val += pow(w[k][t],q) * table[t][i][j];
          num += Util::square(u[k][i] * u[k][j]) * val;
        }
        double den = 0;
        rep(i,N) den += Util::square(u[k][i]);
        den *= 2;
        assert(Util::cmp(den) != 0);
        J1 += num / den;
      }
      double J2 = 0;
      for(set< Pair >::const_iterator x = must_link.begin(); x != must_link.end(); x++) { // must link
        int l = x->first, m = x->second;
        rep(r,C) rep(s,C) if(r != s) J2 += u[r][l] * u[s][m];
      }
      for(set< Pair >::const_iterator x = must_not_link.begin(); x != must_not_link.end(); x++) { // must not link
        int l = x->first, m = x->second;
        rep(r,C) J2 += u[r][l] * u[r][m];
      }
      double new_J = J1 + alpha * J2;
      objective.push_back(J = new_J);
      restriction = J2;
      return true;
    }
    vector< vector< int > > confusing_matrix() {
      int k = C;
      int p = priori_cluster.size();
      vector< vector< int > > table(k+1,vector< int >(p+1,0));
      rep(i,k) rep(j,p) for(Cluster::const_iterator x = cluster[i].begin(); x != cluster[i].end(); x++) table[i][j] += priori_cluster[j].count(*x);
      rep(i,k) table[i][p] = cluster[i].size();
      rep(i,p) table[k][i] = priori_cluster[i].size();
      table[k][p] = N;
      return table;
    }
    bool operator<(const Answer& answer) const {
      return J < answer.J;
    }
  };
  double find_alpha(double start = max_alpha) {
    const double this_alpha = alpha;
    double ans = start, take = start, R= log2(start) + 10;
    rep(i,R) {
      double nans = ans - take; take /= 2;
      if(nans < 0) continue; else alpha = nans;
      int ok = 0;
      rep(j,5) {
        Answer x(-1);
        x.srand();
        rep(r,350) {
          double old_J = x.J;
          if(!x.update()) break;
          double new_J = x.updateJ();
          if(fabs(old_J-new_J) <= eps || x.restriction <= eps) break;
        }
        if(x.restriction <= eps) {
          ok = 1;
          break;
        }
      }
      if(ok) ans = nans;
    }
    alpha = this_alpha;
    return ans;
  }
  Answer main(Latex& latex, bool print = true) {
    const double this_alpha = alpha;
    if(automatic_alpha) {
      double alpha3 = find_alpha(); // FIXME
      alpha = alpha3;
#ifdef DETAILED
      dbg(alpha3);
#endif
    }
    Answer opt(-1);
    Counter::start_timer();
    rep(initialization,initialization_number) {
      if(Counter::elapsed()>TLE) {
#ifdef DETAILED
        cerr << "tle at " << initialization << endl;
#endif
        break;
      }
      Answer now(initialization);
      now.srand();
      rep(iteration,maximum_iteration_number) {
        double old_J = now.J;
        if(!now.update()) break;
        double new_J = now.J;
        if(fabs(old_J - new_J) <= 1e-12) break; else opt = min(opt,now);
      }
    }
    if(print){
      if(automatic_alpha) {
        M["[alpha]"] = vector< string >(1,Util::cast<string>(alpha));
      }
      char text[1<<10];
      string caption;
      vector<string> elem;
      latex.print_name("Another algorithm"); // name
      latex.print_description(M);
      caption="Initialization \\textbf{"+Util::cast<string>(opt.id+1)+"}"; // best initialization
      elem.clear();
      elem.push_back("Iteration");
      elem.push_back("Objective function");
      elem.push_back("~");
      rep(i,opt.objective.size()) {
        elem.push_back(Util::cast<string>(i+1));
        sprintf(text,i+1==(int)opt.objective.size()?"\\textbf{%.7Lf}":"%.7Lf",opt.objective[i]);
        elem.push_back(text);
        sprintf(text,i == 0 ? "~":"(%.7Lf)",opt.objective[i]-opt.objective[i?i-1:0]);
        elem.push_back(text);
      }
      latex.print_long_table(3,caption,elem);
      vector< vector< int > > table = opt.confusing_matrix();
      pair< double, vector< int > > mcmf = Validation::global_error(N,table);
      caption="Confusing matrix";
      elem.clear();
      elem.push_back("~");
      rep(k,C) elem.push_back(Util::cast<string>(k+1));
      rep(k,C) {
        elem.push_back(Util::cast<string>(k+1)+"'");
        rep(h,C) {
          string s = Util::cast<string>(table[k][h]);
          if(mcmf.second[k] == h) s = "\\textbf{" + s + "}";
          elem.push_back(s);
        }
      }
      latex.print_table(C+1,caption,elem);
      elem.clear();
      sprintf(text,"Global error: \\textbf{%.7Lf\\%%}", mcmf.first * 100.0);
      elem.push_back(text);
      latex.print_itemize(elem);
      elem.clear();
      sprintf(text,"Adjusted rand index: \\textbf{%.7Lf}",Validation::adjusted_rand_index(N,table));
      elem.push_back(text);
      sprintf(text,"F measure: \\textbf{%.7Lf}",Validation::f_measure(N,table));
      elem.push_back(text);
      latex.print_itemize(elem);
      if(true) {
        vector< pair< string, vector< string > > > v_elem(C);
        rep(k,C) {
          sprintf(text,"Cluster %d:", k+1);
          pair< string, vector< string > >& p_elem = v_elem[k];
          p_elem.first = text;
          rep(t,T) {
            sprintf(text,"Table %d: \\textbf{%.7Lf}", t+1, opt.W[k][t]);
            p_elem.second.push_back(text);
          }
        }
        latex.print_itemize(v_elem);
      }
      caption = "Fuzzy clustering";
      elem.clear();
      elem.push_back("~");
      rep(k,C) elem.push_back(Util::cast<string>(k+1));
      rep(i,N) {
        sprintf(text,"%.4d",i+1);
        elem.push_back(text);
        rep(k,C) {
          sprintf(text,opt.cluster[k].count(i)?"\\textbf{%.7Lf}":"%.7Lf",opt.U[k][i]);
          elem.push_back(text);
        }
      }
      latex.print_long_table(C+1,caption,elem);
      if(true) {
        latex<<"\\tiny\n";
        stringstream s;
        s << setprecision(6) << fixed << Counter::elapsed() / 6000000.0 << " minutes, with seed = " << seed << ".\n\n";
        latex<<s.str();
      }
      latex<<"\\pagebreak\n";
    }
    alpha = this_alpha;
    return opt;
  }
};

namespace Generator {
  const int class_number = 2, variable_number = 2, pattern_number = 300;
  const int class_size[class_number] = 
  {150,150};
  const double mean[class_number][variable_number] = 
    //{{-0.4,0.1},{0.1,32.0}};
  {{-5.06,0.04},{5.08,0.13}};
  const double deviation[class_number][variable_number] =
    //{{sqrt(236.6),sqrt(1.0)},{sqrt(1.0),sqrt(215.2)}}; 
  {{sqrt(4.79),sqrt(0.70)},{sqrt(0.77),sqrt(4.96)}};
  const double rho[class_number] = 
    //{0.6,-0.2};
  {0.0,0.0};
  vector< pair< double, double > > pattern;
  vector< int > where;
  void bivariate_gaussian(double sigma_x, double sigma_y, double mean_x, double mean_y, double rho, double &x, double &y) {
    double u, v, r2, scale;
    do {
      u = -1 + 2 * Util::f_rand (0.0,1.0);
      v = -1 + 2 * Util::f_rand (0.0,1.0);
      r2 = u * u + v * v;
    } while (r2 > 1.0 || r2 == 0);
    scale = sqrt (-2.0 * log (r2) / r2);
    x = mean_x + sigma_x * u * scale;
    y = mean_y + sigma_y * (rho * u + sqrt(1 - rho*rho) * v) * scale;
  }
  void generate_bivariate_normal_distribution() {
    pattern.resize(pattern_number);
    where.resize(pattern_number);
    rep(i,pattern_number) {
      int &my = where[i] = -1;
      if(i < class_size[0]) my = 0;
      else if(i < class_size[0] + class_size[1]) my = 1;
      else my = 2;
      bivariate_gaussian(deviation[my][0], deviation[my][1], mean[my][0], mean[my][1], rho[my], pattern[i].first, pattern[i].second);
    }
  }
  void prepare() {
    // variables
    C = class_number;
    P = 2;
    initialization_number = 100;
    label = 0.10;
    alpha = 0.0;
    maximum_iteration_number = 150;
    N = pattern_number;
    T = 2;
    // put in scanner
    M["[cluster_number]"].push_back(Util::cast<string>(C));
    M["[prototype_number]"].push_back(Util::cast<string>(P));
    M["[initialization_number]"].push_back(Util::cast<string>(initialization_number));
    M["[label]"].push_back(Util::cast<string>(label));
    M["[alpha]"].push_back(Util::cast<string>(alpha));
    M["[maximum_iteration_number]"].push_back(Util::cast<string>(maximum_iteration_number));
    M["[output]"].push_back("main.tex");
    // generate data
    generate_bivariate_normal_distribution();
    priori_cluster = vector< Cluster >(C, Cluster());
    rep(i,N) priori_cluster[where[i]].insert(i);
    // generate table
    table = vector< Matrix >(T, Matrix(N, Row(N,0.0)));
    rep(t,T) rep(i,N) rep(j,N) {
      double var1 = (t == 0 ? pattern[i].first : pattern[i].second);
      double var2 = (t == 0 ? pattern[j].first : pattern[j].second);
      double dist = Util::square(var1 - var2);
      table[t][i][j] = table[t][j][i] = dist;
    }
  }
  void generate_restrictions() {
    must_link.clear();
    must_not_link.clear();
    map< int, int > mapping;
    rep(k,priori_cluster.size()) for(Cluster::const_iterator x = priori_cluster[k].begin(); x != priori_cluster[k].end(); x++) mapping[*x] = k;
    int n = (N * label + 0.7777777);
    assert(n <= N);
    vector< int > object(N);
    rep(i,N) object[i] = i;
    Util::randomize(object);
    rep(i,n) rep(j,n) if(i != j) {
      int a = object[i], b = object[j];
      if(mapping[a] == mapping[b]) must_link.insert( make_pair(a,b) );
      else must_not_link.insert( make_pair(a,b) );
    }
  }
};

int main() {
  seed = time(NULL);
  MersenneTwister::build(seed);
  dbg(seed);
  Scanner scan("config");
  if(false) {
    Generator::prepare();
    if(true){
      ofstream out("data.in");
      rep(i,N) {
        out << Generator::pattern[i].first << "," << Generator::pattern[i].second << endl;
      }
      out.close();
    }
    if(false) return 0;
  }
  string out = M["[output]"][0];
  ofstream out_file(out.c_str(),ios::out);
  Latex latex(out_file);
  latex.begin();
#ifdef PROTOTYPE_VARIANCE
  for(int this_P=Util::cast<int>(M["[prototype_number_begin]"][0]);this_P<=Util::cast<int>(M["[prototype_number_end]"][0]);this_P++){
    P=this_P;
    M["[prototype_number]"][0] = Util::cast<string>(this_P);
#else
    {
#endif
#ifdef MAGIC_VARIANCE
      for(int this_magic=Util::cast<int>(M["[magic_begin]"][0]);this_magic<=Util::cast<int>(M["[magic_end]"][0]);this_magic++){
        M["[magic]"][0] = Util::cast<string>(this_magic);
#else
        {
#endif
          vector< double > E1, E2;
          vector< vector< int > > tab;
          vector< int > v1, v2;
          int do_manual = int(M["[order]"].size());
          if(do_manual) {
            vector< string > vs = M["[order]"];
            vector< int > vi(vs.size());
            rep(i,vs.size()) vi[i] = Util::cast<int>(vs[i]);
            scan.manual(vi);
            M.erase("[order]");
          }
          int repeat = Util::cast<int>(M["[repeat]"][0]);
          rep(R,repeat) {
            if(!do_manual) scan.generate_restrictions();
            Algorithm::update_coefficient_cluster_table = true;
            Algorithm::Answer a = Algorithm::main(latex);
            tab = a.confusing_matrix();
            v1 = a.get_order();
#ifdef DETAILED
            cerr << "! " << a.restriction << endl;
#endif
            E1.push_back( Validation::global_error(N,tab).first );
#ifdef RUN_ANOTHER
            Another::Answer b = Another::main(latex);
            tab = b.confusing_matrix();
            v2 = b.get_order();
#ifdef DETAILED
            cerr << "? " << b.restriction << endl;
#endif
#endif
            E2.push_back( Validation::global_error(N,tab).first );
#ifdef DETAILED
            dbg( R _ E1.back() _ E2.back() );
#else
            cerr<<"("<<E1.back();
#ifdef RUN_ANOTHER
            <<","<<E2.back();
#endif
            cerr<<"),";
#endif
          }
          cerr<<endl;
          pair< double, double > i1 = Analysis::build(E1);
          pair< double, double > i2 = Analysis::build(E2);
          for(double magic = M.count("[magic]") ? Util::cast<double>(M["[magic]"][0]) : 1e9;;){
            dbg(P _ magic _ i1.first _ i1.second);
#ifdef RUN_ANOTHER
            dbg(P _ magic _ i2.first _ i2.second);
#endif
            break;
          }
#ifdef DETAILED
          cerr << "v1: "; Util::pv(v1.begin(),v1.end());
          cerr << "v2: "; Util::pv(v2.begin(),v2.end());
#endif
          if(true) {
            latex<<"\\pagebreak\n";
            latex<<"\\large\n";
            stringstream s;
            s << setprecision(12) << fixed << "my: (" << i1.first << ":" << i1.second << ")\n";
            s << "\\\\\n";
            s << setprecision(12) << fixed << "his: (" << i2.first << ":" << i2.second << ")\n";
            latex<<s.str();
          }
#if 1
        }
#endif
#if 1
      }
#endif
  latex.end();
  out_file.close();
  return 0;
}

