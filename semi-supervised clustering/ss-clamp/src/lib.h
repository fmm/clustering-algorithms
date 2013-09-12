#ifndef LIB_H_
#define LIB_H_ 
using namespace std;

#include <bits/stdc++.h>
#include <sqlite3.h>
#include <sys/time.h>

// stringify
#define STR(x) #x

// assertive
#define ASSERT(condition,message) __assert(condition, message, __FILE__, __LINE__)

void __assert(bool _condition, string _message, const char* _file, int _line) {
  if(!(_condition)) {
    std::cerr << "ASSERTION FAILED!" << std::endl
      << "File = " << _file << std::endl
      << "Line = " << _line << std::endl
      << "Message = " << _message << std::endl;
#ifndef NDEBUG
    std::exit(EXIT_FAILURE);
#endif     
  }
}

// warning
#define WARNING(message) std::cerr << "WARNING: " << message << std::endl

// alert
#define ALERT(message) std::cerr << "[" << message << "]" << std::endl

// for debug
#ifndef NDEBUG
#define _ << ", " <<
#define dbg(x) std::cerr << STR(x) << " = " << x << std::endl
#define pv(x,y) \
  do { \
    for(decltype(y) __xy = (x); __xy != (y); __xy++) { \
      std::cerr << *__xy << " "; \
    } \
    std::cerr << std::endl; \
  } while(false)
#else
#define _ {}
#define dbg(x) {}
#define pv(x,y) {}
#endif

// for better precision
#define double long double
#define DOUBLE "Lf"

// for validation
#define VALIDATE_FILE(condition,file_name) ASSERT(condition,"invalid file:" + file_name)
#define VALIDATE_DENOMINATOR(x) ASSERT(Util::cmp(x), "division by zero")

// for general usage
typedef vector<double> Row;
typedef vector<Row> Matrix;
typedef unordered_set<unsigned int> Cluster; 
typedef vector<unsigned int> Prototype;
typedef pair<unsigned int,unsigned int> Pair;

// hash function for pair of integers
namespace std {
  template<> struct hash<Pair> {
    inline size_t operator()(const Pair& p) const {
      std::hash<unsigned long long> hasher;
      return hasher((unsigned long long) p.first << 32 ^ p.second);
    }
  };
}

// constants
const int inf = ~0u >> 2;
const double INF = 1e100;

// util functions
namespace Util {
  int cmp(double x, double y = 0, double eps = 1e-16) {
    if(fabs(x-y) <= eps) return 0; else return x < y ? -1 : +1;
  }
  template<class T> T square(T x) {
    return x * x;
  }
  template<class S, class T> S cast(T a) {
    stringstream s;
    s << a;
    S b;
    s >> b;
    return b;
  }
  template<class S, class T> vector<S> cast(vector<T> a) {
    vector<S> s;
    for(unsigned int i = 0; i < a.size(); ++i) {
      s.push_back(cast<S>(a[i]));
    }
    return s;
  }
  bool ends_with(string a, string b) {
    return a.size() >= b.size() and a.substr(a.size() - b.size()) == b;
  }
};

#endif
