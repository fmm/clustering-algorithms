#ifndef PARAMETER_H_
#define PARAMETER_H_

#include "database.h"
#include "lib.h"
#include "mersennetwister.h"
#include "sha1.h"

struct Parameter {
  // database to store results
  Database database;
  // hash string for all parameters
  string sha1;
  // summarized info about the execution
  string info;
  // seed used for the random generator
  time_t seed;
  // index of the position of the priori cluster on "RECTANGLE_MATRIX"
  // 0 for unlabeled data
  unsigned int class_variable;
  // number of initizalizations
  unsigned int initialization;
  // maximum number of iterations during one initizalization
  unsigned int maximum_iteration;
  // number of individuals
  unsigned int N;
  // mask of individuals to be considered
  vector<bool> mask;
  // number of clusters
  unsigned int C;
  // number of prototypes
  unsigned int P;
  // number of dissimilarities matrices
  unsigned int T;
  // importance of the supervision
  double alpha;
  // time limit (in minutes)
  unsigned long long time_limit;
  // eps for objective function
  double eps_for_criterion;
  // relevance factor
  double relevance_v;
  // percentage of labeled data
  double label;
  // priori cluster for labeled data
  vector<Cluster> priori_cluster;
  // pairwise constraint restrictions
  unordered_set<Pair> must_link, cannot_link;
  // source files
  vector<string> input;
  string pwc_file;
  // dissimilarities matrices
  vector<Matrix> table;
  // summary of all variables
  map< string, vector<string> > summary;

  // class constructor
  Parameter(string file_name) {
    // load general variables
    load(file_name);
    // calculate unique id
    generate_id();
  }

  void save() {
    // save algorithm execution
    string sql = "INSERT INTO algorithm("
      "sha1,"
      "info,"
      "seed,"
      "time_limit,"
      "class_variable,"
      "label_percentage,"
      "initializations,"
      "clusters,"
      "individuals,"
      "prototypes,"
      "eps,"
      "relevance_v,"
      "iterations)"
      "VALUES(" +
      string("\"" + sha1 + "\"") + "," +
      string("\"" + info + "\"") + "," +
      Util::cast<string>(seed) + "," +
      Util::cast<string>(time_limit) + "," +
      Util::cast<string>(class_variable) + "," +
      Util::cast<string>(label) + "," +
      Util::cast<string>(initialization) + "," +
      Util::cast<string>(C) + "," +
      Util::cast<string>(N) + "," +
      Util::cast<string>(P) + "," +
      Util::cast<string>(eps_for_criterion) + "," +
      Util::cast<string>(relevance_v) + "," +
      Util::cast<string>(maximum_iteration) +
      ");";
    database.execute(sql);
    // save input files
    for(unsigned int i = 0; i < input.size(); ++i) {
      sql = "INSERT INTO input(algorithm_id,file)VALUES(" +
        string("\"" + sha1 + "\"") + "," +
        string("\"" + input[i] + "\"") +
        ");";
      database.execute(sql);
    }
  }
  
  void update_multable() {
    // update number of prototypes
    summary["[PROTOTYPES]"].clear();
    summary["[PROTOTYPES]"].push_back(Util::cast<string>(P));
    // update the mask
    summary["[MASK]"].clear();
    summary["[MASK]"] = Util::cast<string>(mask);
  }

  void generate_id() {
    update_multable();
    sha1 = "";
    for(auto iter : summary) {
      sha1 += iter.first;
      sha1 += accumulate(iter.second.begin(),iter.second.end(),string(":"));
    }
    sha1 = sha1::process(sha1);  
  }

  void generate_pwc(string pwc_file, double known) {
    this->pwc_file = pwc_file;
    this->label = known;
    input.push_back(pwc_file);
    summary["[INPUT]"].push_back(pwc_file);
    // (individual,class) pairs
    vector<unsigned int> objects(N);
    for(unsigned int i = 0; i < N; ++i) {
      objects[i] = i;
    }
    // randomize the pairs
    MersenneTwister random(time(0));
    for(unsigned int i = 1; i < N; ++i) {
      unsigned int j = random.rand_unsigned() % (i + 1);
      swap(objects[i], objects[j]);
    }
    // load label info
    vector<int> class_id(N, -1);
    for(unsigned int p = 0; p < priori_cluster.size(); ++p) {
      for(unsigned int i = 0; i < N; ++i) {
        if(priori_cluster[p].count(i)) {
          ASSERT(class_id[i] == -1, "clusters should be disjointed");
          class_id[i] = p;
        }
      }
    }
    must_link.clear();
    cannot_link.clear();
    unsigned int pairs = N * known / 100.0;
    ASSERT(pairs <= N, "label was invalid");
    for(unsigned int i = 0; i < pairs; ++i) {
      for(unsigned int j = i+1; j < pairs; ++j) {
        unsigned int a = objects[i], b = objects[j];
        if(class_id[a] == class_id[b]) {
          must_link.insert(Pair(a,b));
          must_link.insert(Pair(b,a));
        } else {
          cannot_link.insert(Pair(a,b));
          cannot_link.insert(Pair(b,a));
        }
      }
    }
    dbg(must_link.size());
    dbg(cannot_link.size());
    ofstream out(pwc_file.c_str(), ios::out);
    VALIDATE_FILE(out, pwc_file);
    out << "PAIRWISE_CONSTRAINTS = (" "\n"
      "INFO = (" "\n"
      "\t" "percentage = " << label << "," "\n"
      "\t" "must_link = " << must_link.size() << "," "\n"
      "\t" "cannot_link = " << cannot_link.size() << "," "\n"
      ")," "\n"
      "MUST_LINK = (" "\n";
    for(unsigned int i = 0; i < pairs; ++i) {
      for(unsigned int j = i+1; j < pairs; ++j) {
        unsigned int a = objects[i], b = objects[j];
        if(class_id[a] == class_id[b]) {
          out << "\t" "(" << a << "," << b << ")," "\n";
          out << "\t" "(" << b << "," << a << ")," "\n";
        }
      }
    }     
    out << ")," "\n"
      << "CANNOT_LINK = (" "\n";
    for(unsigned int i = 0; i < pairs; ++i) {
      for(unsigned int j = i+1; j < pairs; ++j) {
        unsigned int a = objects[i], b = objects[j];
        if(class_id[a] != class_id[b]) {
          out << "\t" "(" << a << "," << b << ")," "\n";
          out << "\t" "(" << b << "," << a << ")," "\n";
        }
      }
    }
    out << "))" "\n"
      << "END\n";
    out.close();
  }

  void set_default_values() {
    if(summary.count("[SEED]") == 0) {
      summary["[SEED]"].push_back(Util::cast<string>(time(NULL)));
    }
    if(summary.count("[EPS_FOR_CRITERION]") == 0) {
      summary["[EPS_FOR_CRITERION]"].push_back(Util::cast<string>((double)(1e-12)));
    }
    if(summary.count("[RELEVANCE_V]") == 0) {
      summary["[RELEVANCE_V]"].push_back(Util::cast<string>((double)(1.0)));
    }
  }

  void load(string file_name) {
    ifstream in(file_name.c_str(), ios::in);
    VALIDATE_FILE(in, file_name);
    // set summary
    summary.clear();
    for(string kind,value; in>>value;) {
      if(value[0] == '[' and value[value.size()-1] == ']') {
        kind = value;
      } else {
        summary[kind].push_back(value);
      }
    }
    in.close();
    // set default values
    set_default_values();
    // set dabasase
    database = Database(summary["[DATABASE]"][0]);
    // set info
    info = accumulate(summary["[INFO]"].begin(),summary["[INFO]"].end(),string(""));
    // set seed
    seed = Util::cast<time_t>(summary["[SEED]"][0]);
    // set class variable
    class_variable = Util::cast<unsigned int>(summary["[CLASS_VARIABLE]"][0]);
    // set number of initizalization
    initialization = Util::cast<unsigned int>(summary["[INITIALIZATION]"][0]);
    // set maximum iteration
    maximum_iteration = Util::cast<unsigned int>(summary["[MAXIMUM_ITERATION]"][0]);
    // set number of classes
    C = Util::cast<unsigned int>(summary["[CLUSTERS]"][0]);
    // set number of prototypes
    P = Util::cast<unsigned int>(summary["[PROTOTYPES]"][0]);
    // set alpha
    alpha = Util::cast<double>(summary["[ALPHA]"][0]);
    // set time limit
    // converted to seconds
    time_limit = 60000000ULL * Util::cast<double>(summary["[TIME_LIMIT]"][0]);
    // set eps for criterion
    eps_for_criterion = Util::cast<double>(summary["[EPS_FOR_CRITERION]"][0]);
    // set relevance factor
    relevance_v = Util::cast<double>(summary["[RELEVANCE_V]"][0]);
    // set input files
    input = summary["[INPUT]"];
    N = 0;
    for(unsigned int i = 0; i < input.size(); ++i) {
      if(Util::ends_with(input[i],".sds")) {
        if(N <= 0) {
          // set number of individuals
          read_individual_number(input[i]);
          // set mask
          mask = vector<bool>(N, true);
          if(class_variable > 0) {
            // set priori cluster
            read_priori_cluster(input[i]);
          }
        }
        // set matrices of dissimilarities
        table.push_back(Matrix(N,Row(N,0.0)));
        read_table(input[i], table.back());
      } else if(Util::ends_with(input[i],".pwc")) {
        pwc_file = input[i];
        // set label and pairwise constraints
        read_pairwise_constraints(input[i]);
      } else ASSERT(false, "file not supported: " + input[i]);
    }
    // set number of tables
    T = table.size();
  }

  void read_individual_number(string file_name) {
    ifstream in(file_name.c_str(), ios::in);
    VALIDATE_FILE(in, file_name);
    string line;
    do VALIDATE_FILE(getline(in,line),file_name);
    while(line.find("indiv_nb") == string::npos);
    ASSERT(
        sscanf(line.c_str()," indiv_nb = %d", &N) == 1,
        "failed to read the number of individuals"
        );
    in.close();
  }

  void read_priori_cluster(string file_name) {
    ifstream in(file_name.c_str(), ios::in);
    VALIDATE_FILE(in, file_name);
    string line;
    do VALIDATE_FILE(getline(in,line),file_name);
    while(line.find("RECTANGLE_MATRIX = (") == string::npos);
    priori_cluster.clear();
    for(unsigned int i = 0; i < N; ++i) {
      line.clear();
      char c;
      int ct = 0;
      while(ct == 0) {
        VALIDATE_FILE(in>>c, file_name);
        if(c == '(') ct++;
      }
      while(ct != 0) {
        VALIDATE_FILE(in>>c, file_name);
        if(c == '(') ct++; else if(c == ')') ct--;
        if(ct != 0) line.push_back(c);
      }
      line.push_back(',');
      vector<string> var;
      string piece;
      for(unsigned int j = 0; j < line.size(); ++j) if(line[j] != ' ') {
        if (line[j] == '(') ct++; else if (line[j] == ')') ct--;
        if(!ct && line[j] == ',') {
          var.push_back(piece), piece.clear();
        } else piece.push_back(line[j]);
      }
      unsigned int priori;
      ASSERT(
          sscanf(var[class_variable-1].c_str(),"%u",&priori) == 1,
          "failed to read the priori cluster"
          );
      while(priori > priori_cluster.size()) priori_cluster.push_back(Cluster());
      priori_cluster[priori - 1].insert(i);
    }
    in.close();
  }

  void read_table(string file_name, Matrix& M) {
    ifstream in(file_name.c_str(), ios::in);
    VALIDATE_FILE(in, file_name);
    string line;
    do VALIDATE_FILE(getline(in,line), file_name);
    while(line.find("DIST_MATRIX= (") == string::npos);
    for(unsigned int i = 0; i < N; ++i) {
      line.clear();
      char c;
      for(c = '?';c != '('; VALIDATE_FILE(in>>c, file_name));
      for(VALIDATE_FILE(in>>c, file_name);c != ')'; VALIDATE_FILE(in>>c, file_name)) line.push_back(c);
      for(unsigned int j = 0; j < line.size(); ++j) {
        if(line[j] == '(' || line[j] == ')' || line[j] == ',') {
          line[j] = ' ';
        }
      }
      stringstream s(line);
      for(unsigned int j = 0; j <= i; ++j) {
        double d;
        VALIDATE_FILE(s>>d, file_name);
#ifdef CHECK_NORMALIZED
        ASSERT(
            0  <= Util::cmp(d,0.0) and Util::cmp(d,1.0) <= 0,
            "table in " + file_name + " is not normalized"
            );
#endif
        M[i][j] = M[j][i] = d;
      }
    }
    in.close();
  }

  void read_pairwise_constraints(string file_name) {
    ifstream in(file_name.c_str(), ios::in);
    VALIDATE_FILE(in, file_name);
    string line;
    do VALIDATE_FILE(getline(in,line), file_name);
    while(line.find("INFO = (") == string::npos);
    // set percentage of labeled data
    VALIDATE_FILE(getline(in,line), file_name);
    ASSERT(
        sscanf(line.c_str()," percentage = %" DOUBLE, &label) == 1,
        "failed to read the percentage of labeled data"
        );
    unsigned int n_must_link, n_cannot_link;
    VALIDATE_FILE(getline(in,line), file_name);
    ASSERT(
        sscanf(line.c_str()," must_link = %u", &n_must_link) == 1,
        "failed to read the number of must-link pairs"
        );
    VALIDATE_FILE(getline(in,line), file_name);
    ASSERT(
        sscanf(line.c_str()," cannot_link = %u", &n_cannot_link) == 1,
        "failed to read the number of cannot-link pairs"
        );
    // set must-link constraints
    do VALIDATE_FILE(getline(in,line), file_name);
    while(line.find("MUST_LINK = (") == string::npos);
    must_link.clear();
    for(unsigned int i = 0; i < n_must_link; ++i) {
      VALIDATE_FILE(getline(in,line), file_name);
      unsigned int a, b;
      ASSERT(
          sscanf(line.c_str(), " (%u,%u)", &a,&b) == 2,
          "failed to read the must-link constraint: " + line
          );
      must_link.insert(make_pair(a,b));
      must_link.insert(make_pair(b,a));
    }
    // set cannot-link constraints
    do VALIDATE_FILE(getline(in,line), file_name);
    while(line.find("CANNOT_LINK = (") == string::npos);
    cannot_link.clear();
    for(unsigned int i = 0; i < n_cannot_link; ++i) {
      VALIDATE_FILE(getline(in,line), file_name);
      unsigned int a, b;
      ASSERT(
          sscanf(line.c_str(), " (%u,%u)", &a,&b) == 2,
          "failed to read the cannot-link constraint: " + line
          );
      cannot_link.insert(make_pair(a,b));
      cannot_link.insert(make_pair(b,a));
    }
    in.close();
  }

};

#endif
