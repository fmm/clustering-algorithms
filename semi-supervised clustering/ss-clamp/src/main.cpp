#include "lib.h"

#include "ssclampproductglobal.h"
#include "ssclampproductlocal.h"
#include "ssclampsumglobal.h"
#include "ssclampsumlocal.h"

int main(int argc, char *argv[]) {

  string config;
  unsigned int known = 0;
  string method;
  string output;

  for(int opt;(opt = getopt(argc, argv, "h:g:m:o:")) != -1;) {
    switch(opt) {
      case 'h':
        // TODO: help
        break;
      case 'g':
        known = Util::cast<int>(optarg);
        dbg(known);
        break;
      case 'm':
        // TODO: choose method here
        break;
      case 'o':
        output = Util::cast<string>(optarg);
        dbg(output);
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
    params.generate_pwc(pwc, known);    
    params.generate_id();
  }

  dbg(params.N);
  dbg(params.sha1);
  dbg(params.seed);

  output = params.sha1;
  
  string tex = output + ".tex";
  dbg(tex);

  //return 0;

  //SSClamp algo(p);
  //SSClampProductGlobal algo(p);
  SSClampProductLocal algo(params);
  //SSClampSumGlobal algo(p); // best weight
  //SSClampSumLocal algo(p);
  
  Method::Answer ans = algo.process();

  vector< vector<unsigned int> > C = algo.compute_confusion_matrix(ans);
  
  if(true) for(unsigned int i = 0; i < C.size(); ++i) {
    for(unsigned int j = 0; j < C[i].size(); ++j) {
      cout << setw(3) << C[i][j] << " ";    
    }
    cout << endl;
  }
  
  dbg(Validation::accuracy(C).first);
  
  //dbg( Validation::fuzzy_rand_index_hullermeier(ans.U,ans.U) );

  return 0;


	/*
	MersenneTwister rnd(1028);

	dbg(rnd.rand_signed());

  Database db("test.db");

	db.execute("select * from config;");

	double d = (double)(0);

	dbg(d);

	//dbg( sha1::process("") );
	*/

  //start_timer();

	/*MinimumCostMaximumFlow mcmf;

	mcmf.add(1,2,1,1);

	Pair p = mcmf.process(1,2);

	dbg(p.first _ p.second);*/
	
	/*time_t t = -1;

  dbg(t);
  
  var v=var(0);
  
  dbg(v);
  
  string x=str(0+1+...+n);
  dbg(x);*/
  
  //Parameter p("x"); 
  
  //int x=12;
  
  //ASSERT(x == 13, string("x must be 12 ") + "OK?!" );
  
  //dbg(x);
  
  /*
  
  Parameter p("config.txt");
  dbg(p.N);
  dbg(p.sha1);
  dbg(p.seed);
  
  */
  
  //
  
  /*dbg(ans.cluster.size());
  for(unsigned int i = 0; i < p.N; ++i) {
    int clus=0;
    while(!ans.cluster[clus].count(i)) ++clus;
    dbg(i _ clus);
  }
  dbg(ans.criterion _ ans.restriction);*/
  
  //Method::load_parameters(p);
  //Method::Answer a(1);
  
  /*;*/
  
  //A a;
  //A::my m;
  
	return 0;
}

