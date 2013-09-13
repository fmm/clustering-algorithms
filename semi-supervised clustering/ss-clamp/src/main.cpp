// everything
#include "lib.h"

// mine
#include "ssclampproductglobal.h" // prg
#include "ssclampproductlocal.h" // prl 
#include "ssclampsumglobal.h" // smg
#include "ssclampsumlocal.h" // sml

// another
#include "sscard.h" // crd

int main(int argc, char *argv[]) {
  string config;
  double known = 0;
  string method;
  string output;

  for(int opt;(opt = getopt(argc, argv, "h:g:m:o:")) != -1;) {
    switch(opt) {
      case 'h':
        // TODO: help
        break;
      case 'g':
        known = Util::cast<double>(optarg);
        dbg(known);
        break;
      case 'm':
        method = Util::cast<string>(optarg);
        break;
      case 'o':
        output = Util::cast<string>(optarg);
        break;
      default:
        ASSERT(false,"correct usage is: config [-h] [-g known] [-m method] [-o output]");
    }
  }

  ASSERT(optind < argc,"expected argument after options");

  config = Util::cast<string>(argv[optind]);
  dbg(config);
  
  Parameter params(config);
  
  if(known > 0) {
    ASSERT(output.size(),"one output file must be specified");
    string pwc = output;
    if(!Util::ends_with(pwc,".pwc")) {
      pwc = pwc + ".pwc";
    }
    dbg(pwc);
    params.generate_pwc(pwc, known);    
    params.generate_id();
  }

  if(method.size()) {
    Method *machine = NULL;

    if(method.find("prg") != string::npos) {
      assert(machine == NULL);
      machine = new SSClampProductGlobal(params);
    }

    if(method.find("prl") != string::npos) {
      assert(machine == NULL);
      machine = new SSClampProductLocal(params);
    }

    if(method.find("smg") != string::npos) {
      assert(machine == NULL);
      machine = new SSClampSumGlobal(params);
    }

    if(method.find("sml") != string::npos) {
      assert(machine == NULL);
      machine = new SSClampSumLocal(params);
    }

    if(method.find("crd") != string::npos) {
      assert(machine == NULL);
      machine = new SSCARD(params);
    }

    assert(machine != NULL);

    Method::Answer answer = machine->process();
    string tex = output + ".tex";
    dbg(tex);

    vector< vector<unsigned int> > confusing_matrix = machine->compute_confusion_matrix(answer);
    Matrix priori_matrix = machine->compute_priori_matrix();
    
    if(true) for(unsigned int i = 0; i < confusing_matrix.size(); ++i) {
      for(unsigned int j = 0; j < confusing_matrix[i].size(); ++j) {
        cout << setw(3) << confusing_matrix[i][j] << " ";    
      }
      cout << endl;
    }

    dbg( answer.criterion );
    dbg( Validation::accuracy(confusing_matrix).first );
    dbg( Validation::adjusted_rand_index(confusing_matrix) );
    dbg( Validation::f_measure(confusing_matrix) );
    dbg( Validation::fuzzy_rand_index_campello(answer.U,priori_matrix) );
    dbg( Validation::fuzzy_rand_index_hullermeier(answer.U,priori_matrix) );
  
    if(machine != NULL) {
      delete machine;
      machine = NULL;
    }
  }
  
  return 0;
}

