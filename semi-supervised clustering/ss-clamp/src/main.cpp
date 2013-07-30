#include "lib.h"
#include "database.h"
#include "mersennetwister.h"
#include "mersennetwister.h"
#include "minimumcostmaximumflow.h"
#include "method.h"
#include "parameter.h"

#include "ssclamp.h"
#include "validation.h"

#include "ssclampproductglobal.h"
#include "ssclampproductlocal.h"

int main() {

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
  
  
  Parameter p("config.txt");
  dbg(p.N);
  dbg(p.sha1);
  
  //SSClamp algo(p);
  SSClampProductGlobal algo(p);
  //SSClampProductLocal algo(p);
  
  Method::Answer ans = algo.process();
  
  /*dbg(ans.cluster.size());
  for(unsigned int i = 0; i < p.N; ++i) {
    int clus=0;
    while(!ans.cluster[clus].count(i)) ++clus;
    dbg(i _ clus);
  }
  dbg(ans.criterion _ ans.restriction);*/
  
  //Method::load_parameters(p);
  //Method::Answer a(1);
  
  vector< vector<unsigned int> > C = algo.compute_confusion_matrix(ans);
  
  if(true) for(unsigned int i = 0; i < C.size(); ++i) {
    for(unsigned int j = 0; j < C[i].size(); ++j) {
      cout << setw(3) << C[i][j] << " ";    
    }
    cout << endl;
  }
  
  dbg(Validation::accuracy(C).first);
  
  //A a;
  //A::my m;
  
	return 0;
}

