main(){/*
if(false) {
	vector< int > new_prototype(N);
	rep(i,N) new_prototype[i] = i;
	Util::randomize(new_prototype);
	rep(k,C) rep(p,P) prototype[k].push_back(new_prototype[k*P+p]);
}
else { // explore
	rep(i,C) prototype[i]= vector<int>(P,-1);
  vector< bool > used(N,false);
	used[ (prototype[0][0] = get_random() % N) ] = true;
  rep(i,C-1) {
    double dd=0; int ii=-1;
    rep(j,N) {
      if(used[j]) continue;
      double ndd=0;
      rep(k,i+1) rep(t,T) ndd += table[t][ prototype[k][0] ][j];
      if(ndd>dd)dd=ndd,ii=j;
    }
    used[ (prototype[i+1][0] = ii) ] = true;
  }
  rep(p,P-1) rep(i,C) {
    double dd=0; int ii=-1;
    rep(j,N) {
      if(used[j]) continue;
      double ndd=0;
      rep(k,C) if(k!=i) rep(pp,p+1) rep(t,T) ndd += table[t][prototype[k][pp]][j];
      if(ndd>dd)dd=ndd,ii=j;
    }
    used[ (prototype[i][p+1] = ii) ] = true;
  }
}
*/}
