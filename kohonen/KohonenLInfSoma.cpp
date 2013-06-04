using namespace std;
#include <bits/stdc++.h>
#define show(msg) cerr << "ATTENTION! report: " << msg << ", in line(" << __LINE__ << ")" << endl
#define bug(msg) cerr << "bug report: " << msg << ", in line(" << __LINE__ << ")" << endl, exit(1)
#define _ << " _ " <<
#define dbg(x) cerr << #x << " == " << x << endl
#define mp(x,y) make_pair((x),(y))
#define pv(x,y) {for(typeof(y) z=(x);z!=(y);z++)cerr<<*z<<" ";cerr<<endl;}
#define rep(x,y) for(int(x)=(0);(x)<int(y);++(x))
#define x first
#define y second

namespace Util {
	int cmp(double x, double y = 0, double eps = 1e-7) {
		if(fabs(x-y) <= eps) return 0; else return x < y ? -1 : +1;
	}
};

namespace MinCostMaxFlow {
	vector< int > from, to, w, cap, ant;
	vector< int > adj;
	int n, m;

	void add(int x, int y, int z, int c) {
		#define pb push_back
		while(x+1 > n) n++, adj.pb(-1); while(y+1 > n) n++, adj.pb(-1);
		from.pb(x), to.pb(y), w.pb(z), cap.pb(c), ant.pb(adj[x]), adj[x] = m++;
		swap(x,y), z = -z, c = 0;
		from.pb(x), to.pb(y), w.pb(z), cap.pb(c), ant.pb(adj[x]), adj[x] = m++;
		#undef pb
	}

	const int inf = 0x3f3f3f3f;
	vector< int > dist, pot, pai;
	set< pair<int,int> > heap;

	void update(int no, int ndist, int p) {
		if(ndist >= dist[no]) return;
		if(dist[no] < inf) heap.erase(pair<int,int>(dist[no],no));
		dist[no] = ndist, pai[no] = p;
		heap.insert(pair<int,int>(dist[no],no));
	}

	pair<int,int> top() {
		pair<int,int> ret = *heap.begin();
		heap.erase(heap.begin());
		return ret;
	}

	int djikstra(int source, int sink) {
		heap.clear(), dist = vector< int >(n,inf), pai = vector< int >(n,-1);
		update(source,0,-1);
		while(heap.size()) {
			pair<int,int> p = top();
			for(int i = adj[p.second]; i >= 0; i = ant[i]) if(cap[i]) {
				update(to[i], p.first + w[i] + pot[p.second] - pot[to[i]], i);
			}
		}
		return dist[sink] < inf;
	}

	pair<int,int> mcmf(int source, int sink) {
		pair<int,int> p(0,0); // cost,flow
		pot = vector< int >(n,inf); pot[source] = 0;
		for(int k = 0; k < n; k++) for(int i = 0; i < n; i++) for(int j = adj[i]; j >= 0; j = ant[j]) if(cap[j]) {
			pot[to[j]] = min(pot[to[j]], pot[i] + w[j]);
		};
		while(djikstra(source,sink)) {
		  int cost = 0, flow = inf;
		  for(int x = sink; x != source; x = from[pai[x]]) if(cap[pai[x]] < flow) flow = cap[pai[x]];
		  for(int x = sink; x != source; x = from[pai[x]]) cap[pai[x]] -= flow, cap[pai[x]^1] += flow, cost += w[pai[x]]*flow;
		  for(int x = 0; x < n; x++) pot[x] += dist[x]; p.first += cost, p.second += flow;
		}
		return p;
	}

	void init() {
		n = m = 0;
		adj.clear(), from.clear(), to.clear(), w.clear(), cap.clear(), ant.clear();
	}
};

namespace Validation {
	double globalError(int n, vector< vector<int> >& tab) {
		int N = tab.size() - 1, M = tab[0].size() - 1, source = 0, sink = N + M + 1;
		MinCostMaxFlow::init();
		rep(i,N) rep(j,M) MinCostMaxFlow::add(i + 1, j + N + 1, -tab[i][j], 1);
		rep(i,N) MinCostMaxFlow::add(source, i + 1, 0, 1);
		rep(i,M) MinCostMaxFlow::add(i + N + 1, sink, 0, N);;
		pair<int, int> ret = MinCostMaxFlow::mcmf(source, sink);
		assert(ret.second == N);
		return 1 + double(ret.first) / n;
	}
	
	double adjustedRandIndex(int n, vector< vector<int> >& tab) {
		#define comb(x) ((x)*((x)-1)/2)
		double term[4] = {0}, temp[2] = {0};
		double pot = 1.0 / comb(n);
		int N = tab.size() - 1, M = tab[0].size() - 1;
		rep(i,N) rep(j,M) term[1] += comb(tab[i][j]);
		rep(i,N) temp[0] += comb(tab[i][M]);
		rep(i,M) temp[1] += comb(tab[N][i]);
		term[2] = pot * (temp[0] * temp[1]);
		term[3] = 0.5 * (temp[0] + temp[1]);
		#undef comb
		return (term[1]-term[2])/(term[3]-term[2]);
	}
	
	double fMeasure(int n, vector< vector<int> >& tab) {
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
		return F/n;
	}
};

namespace Kohonen {
	struct Item {
		vector< pair<double,double> > v; // #p
		pair<double,double>& operator[](int n) {
			 return v[n];
		}
		void add(string s) {
			double x, y; assert(sscanf(s.c_str(),"(%lf:%lf)",&x,&y) == 2);
			v.push_back( pair<double,double>(x,y) );
		}	
	};

	// var
	vector< Item > individual; // #n
	vector< Item > prototype;  // #k
	vector< vector< double > > lambda; // #k * #p
	int n, k, p;
	int niter, ninit;
	int dx, dy;
	int classVariable;
	double tmin, tmax;
	int sumFactor;

	vector< int > who, from;
	vector< pair<int,int> > localization;
	vector< vector< double > > neighborhood;	
	vector< vector< int > > table;
	
	double opt;
	int which;
	vector< int > res;
	vector< vector< double > > relevance;
	vector< Item > what;
	
	ofstream out;
	string input, output;

	// read the data
	int read(string file) {
		ifstream config(file.c_str());
		map<string,string> mapa;
	
		for(string label, val; config >> label >> val;) mapa[label] = val;
	
		config.close();
	
		assert(sscanf(mapa["(Number_of_initializations)"].c_str(),"%d",&ninit) == 1);
		assert(sscanf(mapa["(Class_variable)"].c_str(),"%d",&classVariable) == 1); classVariable--;
		assert(sscanf(mapa["(Number_of_individuals)"].c_str(),"%d",&n) == 1);
		assert(sscanf(mapa["(Number_of_iterations)"].c_str(),"%d",&niter) == 1);
		if(sscanf(mapa["(Tmin)"].c_str(),"%lf",&tmin) != 1 || sscanf(mapa["(Tmax)"].c_str(),"%lf",&tmax) != 1) tmin = tmax = 0.0;
		assert(sscanf(mapa["(Sum_factor)"].c_str(),"%d",&sumFactor) == 1);
		assert(sscanf(mapa["(Grid_dimension)"].c_str(),"%dx%d",&dx,&dy) == 2); k = dx*dy;
		input = mapa["(Input_file)"];
		output = mapa["(Output_file)"];
	
		ifstream data(input.c_str());
		string line;
	
		do {
			if (!getline(data, line)) bug("Invalid file!");
		} while (line.find("RECTANGLE_MATRIX = (") == string::npos);
	
		individual.clear();
		
		from = vector< int >(n, -1);
	
		rep(i,n) {
			line.clear();
			char c = '?';
			int ct = 0;
			while(!ct) {
				if(!(data >> c)) bug("Invalid file!");
				if(c == '(') ct++;
			}
			while(ct) {
				if(!(data >> c)) bug("Invalid file!");
				if (c == '(') ct++; else if (c == ')') ct--;
				if (ct) line.push_back(c);
			}
			line.push_back(',');
			vector<string> var;
			string piece;
			rep(j,line.size()) if(line[j] != ' ') {
				if (line[j] == '(') ct++; else if (line[j] == ')') ct--;
				if(!ct && line[j] == ',') {
					var.push_back(piece), piece.clear();
				} else piece.push_back(line[j]);
			}
			Item item;
			rep(j,var.size()) if(j != classVariable) item.add(var[j]); else assert(sscanf(var[j].c_str(),"%d",&from[i]) == 1);
			individual.push_back(item);
		}
	
		data.close();
	
		p = individual[0].v.size();
		
		return 1;
	}
	
	void open() {
		out.open(output.c_str());
	}
	
	void close() {
		out.close();
	}

	// neighborhood function between local 'x' and local 'y' in a given temperature
	double kernel(int x, int y, double temperature) {
		return exp(-neighborhood[x][y] / temperature / temperature / 2.0);
	}

	// dissimilarity between individual 'x' and prototype 'y'
	double func(int x, int y) {
		double dist = 0;
		rep(j,p) {
			double mxj = (individual[x][j].second + individual[x][j].first) / 2.0;
			double rxj = (individual[x][j].second - individual[x][j].first) / 2.0;
			double myj = (prototype[y][j].second + prototype[y][j].first) / 2.0;
			double ryj = (prototype[y][j].second - prototype[y][j].first) / 2.0;
			double d = fabs(mxj - myj) + fabs(rxj - ryj);
			dist += lambda[y][j] * d;
		}
		return dist;
	}

	// the criterion function
	double energy(double temperature) {
		double criterion = 0;
		rep(i,n) rep(j,k) {
			criterion += kernel(who[i],j,temperature) * func(i,j);
		}
		return criterion;
	}

	// initial phase
	void initialization() {
	
		lambda = vector< vector< double > >(k, vector< double >(p,1.0));
	
		prototype = individual;
		random_shuffle(prototype.begin(),prototype.end());
		prototype.resize(k);
			
		localization.clear();
		rep(i,dx) rep(j,dy) {
			localization.push_back( pair<int,int>(i,j) );
		}
	
		neighborhood.clear();
		rep(i,k) {
			neighborhood.push_back( vector<double>() );
			rep(j,k) {
				neighborhood[i].push_back(pow(localization[i].first - localization[j].first, 2.0) + pow(localization[i].second - localization[j].second, 2.0));
			}
		}
	
		// auto values for tmin and tmax
		if(tmin == tmax) {
			tmin = 0.190239866550812569689909992121, tmax = 0;
			rep(i,k) rep(j,k) tmax = max(tmax, neighborhood[i][j]);
			tmax = min(100.0, sqrt(tmax / 2.0 / 1e-6));
		}
		
		dbg(tmin _ tmax);
		
		who.clear();
		rep(i,n) {
			int place = -1;
			double best = 1e111;
			rep(j,k) {
				double nbest = 0;
				rep(c,k) nbest += kernel(j,c,tmax) * func(i,c); 
				if(nbest < best) best = nbest, place = j;
			}
			who.push_back(place);
		}
	}

	// auxiliar function
	double ternarySearch(vector< double > &a, vector< double >& b, int maxIter = 200, double inf = 1e9, double eps = 1e-12) { // a[i] - b[i]*x
		long double low = -inf, high = +inf;
		rep(iter,maxIter) {
			long double med1 = (low * 2 + high) / 3, f1 = 0;
			rep(i,n) f1 += fabs(a[i] - b[i] * med1);
			long double med2 = (low + high * 2) / 3, f2 = 0;
			rep(i,n) f2 += fabs(a[i] - b[i] * med2);
			if(f1 >= f2) low = med1; else high = med2;
		}
		assert(fabs(low-high) <= eps);
		return (low + high) / 2.0;
	}

	// step 1
	void update(double temperature) {
		vector<double> a(n), b(n);
		rep(c,k) rep(j,p) {			
			rep(i,n) b[i] = kernel(who[i],c,temperature);
			rep(i,n) {
				double mxi = (individual[i][j].second + individual[i][j].first) / 2.0;
				a[i] = mxi * b[i];
			}
			double rcj = ternarySearch(a,b); // mi
			rep(i,n) {
				double rxi = (individual[i][j].second - individual[i][j].first) / 2.0;
				a[i] = rxi * b[i];
			}
			double mcj = ternarySearch(a,b); // p
			prototype[c][j].first = rcj - mcj;
			prototype[c][j].second = rcj + mcj;			
		}
	}
	
	// step2
	void fixLocal(double temperature) {
		vector< vector< double > > f = vector< vector< double > >(k, vector< double > (p, 0));
		rep(c,k) rep(j,p) rep(i,n) {
			double mxj = (individual[i][j].second + individual[i][j].first) / 2.0;
			double rxj = (individual[i][j].second - individual[i][j].first) / 2.0;
			double myj = (prototype[c][j].second + prototype[c][j].first) / 2.0;
			double ryj = (prototype[c][j].second - prototype[c][j].first) / 2.0;
			double d = fabs(mxj - myj) + fabs(rxj - ryj);
			f[c][j] += kernel(who[i],c,temperature) * d;
		}
		vector< vector< double > > nlambda = lambda;
		rep(c,k) rep(j,p) {
			nlambda[c][j] = 0.0;
			rep(h,p) {
				if(Util::cmp(f[c][h]) == 0) {
					show("fixLocal: division by zero, returning...");
					return;
				}
				nlambda[c][j] += pow(f[c][j] / f[c][h], 1.0 / (sumFactor-1));
			}
			if(Util::cmp(nlambda[c][j]) == 0) {
				show("fixLocal: division by zero, returning...");
				return;
			}
			nlambda[c][j] = 1.0 / nlambda[c][j];
		}
		rep(c,k) rep(j,p) lambda[c][j] = pow(nlambda[c][j], sumFactor);
	}
	
	void fixGlobal(double temperature) {
		vector< double > f(p, 0);
		rep(c,k) rep(j,p) rep(i,n) {
			double mxj = (individual[i][j].second + individual[i][j].first) / 2.0;
			double rxj = (individual[i][j].second - individual[i][j].first) / 2.0;
			double myj = (prototype[c][j].second + prototype[c][j].first) / 2.0;
			double ryj = (prototype[c][j].second - prototype[c][j].first) / 2.0;
			double d = fabs(mxj - myj) + fabs(rxj - ryj);
			f[j] += kernel(who[i],c,temperature) * d;
		}
		vector< vector< double > > nlambda = lambda;
		rep(c,k) rep(j,p) {
			nlambda[c][j] = 0.0;
			rep(h,p) {
				if(Util::cmp(f[h]) == 0) {
					show("fixLocal: division by zero, returning...");
					return;
				}
				nlambda[c][j] += pow(f[j] / f[h], 1.0 / (sumFactor-1));
			}
			if(Util::cmp(nlambda[c][j]) == 0) {
				show("fixLocal: division by zero, returning...");
				return;
			}
			nlambda[c][j] = 1.0 / nlambda[c][j];
		}
		rep(c,k) rep(j,p) lambda[c][j] = pow(nlambda[c][j], sumFactor);
	}

	// step3
	int partition(double temperature) {
		int change = 0;
		rep(i,n) {
			int place = who[i];
			double best = 0;
			rep(c,k) best += kernel(place,c,temperature) * func(i,c);
			rep(j,k) {
				double nbest = 0;
				rep(c,k) nbest += kernel(j,c,temperature) * func(i,c);
				if(Util::cmp(nbest,best) < 0) best = nbest, place = j;
			}
			if(place != who[i]) change = 1, who[i] = place;
		}
		//vector<int> freq(k,0); rep(i,n) freq[who[i]]++; pv(freq.begin(),freq.end());
		return change;
	}
	
	void makeTable() {
		int N = k, M = 0;
		map<int,int> mapa; rep(i,n) mapa[from[i]] = 1;
		for(map<int,int>::iterator iter = mapa.begin(); iter != mapa.end(); iter++) iter->second = M++;
		table = vector< vector< int > >(N+1, vector< int >(M+1,0));
		rep(i,N) rep(j,M) rep(k,n) table[i][j] += (who[k] == i && mapa[from[k]] == j);
		rep(i,N) rep(j,M) table[i][M] += table[i][j];
		rep(j,M) rep(i,N) table[N][j] += table[i][j];
		table[N][M] = n;
	}

	void optimize(double temperature, int nround = 1e9) {
    do {
			update(temperature);
			#ifdef LOCAL
				fixLocal(temperature);
			#elif GLOBAL
				fixGlobal(temperature);
			#endif
			cerr<<".";
    } while(partition(temperature) && nround--);
    cerr<<endl;
  }

	double process(int init, ostream& fout) {
		const static string space(1<<6,'-');
		#ifdef DETAILED
			fout << fixed << setprecision(3);
			fout << "INITIALIZATION #" << init << "\n" << space << endl;
		#endif
		double temperature = tmax, skip = pow(tmin / tmax, 1.0 / niter), ant = 1e111, prox;
		initialization();
		#ifdef OPTIMIZE
		  optimize(temperature);
		#endif
		ant = prox = energy(temperature);
		for(int iteration = 1; iteration <= niter; iteration++) {
			dbg(iteration);
			#ifdef DETAILED
				makeTable();
				fout << "ITERATION " << setw(3) << iteration << " ::" << endl;
				fout << "* Criterion: " << energy(temperature) << endl;
				fout << "* Global error: " << Validation::globalError(n,table) << endl;
				fout << "* Adjusted rand index: " << Validation::adjustedRandIndex(n,table) << endl;
				fout << "* F measure: " << Validation::fMeasure(n,table) << endl;
				fout << space << endl;
			#endif
			temperature *= skip;
			#ifdef OPTIMIZE
			  optimize(temperature);
			#else
			  optimize(temperature,0);
			#endif
			prox = energy(temperature);
			if(Util::cmp(prox,ant) > 0) dbg(prox _ ant _ fabs(prox-ant));
		}
		return prox;
	}
	
	void solve(ostream& fout = out) {
		const static string space(1<<6,'#');
		opt = 1e111;
		rep(init,ninit) {
			dbg(init);
			double nopt = process(init+1,fout);
			if( Util::cmp(nopt,opt) < 0 ) {
				opt = nopt, res = who, relevance = lambda, which = init+1, what = prototype;
			}
			#ifdef DETAILED
				fout << "\n" << space << "\n" << endl;
			#endif
		}
		who = res;
		prototype = what;
		lambda = relevance;
		makeTable();
		fout << "(Better solution in initialization #" << which << ")" << endl;
		fout << "* Confusion's table:";
		int nd = ceil(log10(n)) + 1;
		rep(i,(fout << "\n",table.size())) rep(j,table[i].size()) fout << setw(nd) << table[i][j];
		fout << "* Criterion: " << opt << endl;
		fout << "* Global error: " << Validation::globalError(n,table) << endl;
		fout << "* Adjusted rand index: " << Validation::adjustedRandIndex(n,table) << endl;
		fout << "* F measure: " << Validation::fMeasure(n,table) << endl;
		#ifdef LOCAL
			fout << "* Releance for each prototype:";
			rep(i,(fout << "\n",relevance.size())) rep(j,relevance[i].size()) fout << relevance[i][j] << "  ";
		#elif GLOBAL
			fout << "* Releance for each variable:\n";
			rep(j,relevance[0].size()) fout << relevance[0][j] << "\n";
		#endif
		fout << "\n" << space << endl;
		fout << "* Cluster for each individual:" << endl;
		rep(i,n) fout << i+1 << " : " << who[i] << endl;
		fout << "* Individual who are more alike to a cluster:" << endl; {
			vector< pair<double,int> > minDist(k, pair<double,int>(1e111,-1) );
			rep(i,n) {
				double newDist = func(i, who[i]);
				if(Util::cmp(newDist,minDist[who[i]].first) < 0) minDist[who[i]] = pair<double,int>(newDist,i);
			}
			rep(i,k) {
				fout << "Cluster " << i+1 << ":" << endl;
				stringstream s1, s2;
				s1 << " prototype:";
				s2 << "individual more similar [" << (minDist[i].second + 1) << "]:";
				string ss1 = s1.str(), ss2 = s2.str();
				while(ss1.size() < ss2.size()) ss1 = string(">") + ss1;
				fout << "\t>>"; fout << ss1;
				rep(j,p) fout << " (" << what[i][j].first << "," << what[i][j].second << ")"; fout << endl;
				fout << "\t< "; fout << ss2;
				if(minDist[i].second != -1) rep(j,p) fout << " (" << individual[minDist[i].second][j].first << "," << individual[minDist[i].second][j].second << ")"; fout << endl;
			}
		}
		fout << "* Grid:"; {
			vector< vector< int > > cnt(dx, vector< int >(dy, 0) );
			rep(i,n) {
				pair<int,int> where = localization[ who[i] ];
				cnt[where.first][where.second]++;
			}
			rep(i,(fout << "\n",dx)) rep(j,dy) fout << " " << setw(3) << cnt[i][j];
		}
		fout << "Topology Error (simple): "; {
			map< int, map<int,int> > adjs;
			rep(i,localization.size()) rep(j,localization.size()) {
				int di = localization[i].first - localization[j].first, dj = localization[i].second - localization[j].second;
				int dij = max(di,-di) + max(dj,-dj);
				if(dij <= 1) adjs[i][j] = adjs[j][i] = 1;
			}
			int wrong = 0;
			rep(i,n) {
				vector< pair<double,int> > my(k);
				rep(j,k) my[j] = pair<double,int>(func(i,j),j);
				sort(my.begin(),my.end());
				int ok1 = my.size() < 1, ok2 = my.size() < 2;
				if(!ok1) ok1=adjs[ who[i] ][ my[0].second ];
				for(int j=1;!ok2&&j<my.size()&&fabs(my[j].first-my[1].first)<=1e-10;j++) ok2=adjs[ who[i] ][ my[j].second ];
				if(!ok1) wrong++; else if(!ok2) wrong++;
			}
			fout << wrong * 100.0 / n << "%" << endl;
		}
		fout << "Topology Error (generalized distance): "; {
			map< int, map<int,int> > adjs;
			rep(i,localization.size()) rep(j,localization.size()) {
				int di = localization[i].first - localization[j].first, dj = localization[i].second - localization[j].second;
				int dij = max(di,-di) + max(dj,-dj);
				if(dij <= 1) adjs[i][j] = adjs[j][i] = 1;
			}
			int wrong = 0;
			rep(i,n) {
				vector< pair<double,int> > my(k);
				rep(j,k) {
					double nbest = 0;
					rep(c,k) nbest += kernel(j,c,tmin) * func(i,c);
					my[j] = pair<double,int>(nbest, j);
				}
				sort(my.begin(),my.end());
				int ok1 = my.size() < 1, ok2 = my.size() < 2;
				if(!ok1) ok1=adjs[ who[i] ][ my[0].second ];
				for(int j=1;!ok2&&j<my.size()&&fabs(my[j].first-my[1].first)<=1e-10;j++) ok2=adjs[ who[i] ][ my[j].second ];
				if(!ok1) wrong++; else if(!ok2) wrong++;
			}
			fout << wrong * 100.0 / n << "%" << endl;
		}
	}
};

int main(int argc, char **argv) {
	
	string method;
	
	#ifdef LOCAL
		method = "L1-Local";
	#elif GLOBAL
		method = "L1-Global";
	#else
		method = "L1-Normal";
	#endif
	
	cout << "now running \"" << method << "\"" << endl;
	
	method = string("Config/") + method + string("-Config.txt");
	
	if(argc > 1) method = string(argv[1]);
	
	dbg(method);
	
	srand(time(NULL));
	Kohonen::read(method);
	Kohonen::open();
	Kohonen::solve();
	Kohonen::close();

	return 0;
}
