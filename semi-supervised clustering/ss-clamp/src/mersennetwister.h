#ifndef MERSENNE_TWISTER_H_
#define MERSENNE_TWISTER_H_
#define MT_N 624
#define MT_M 397
#define MT_MSB 0x80000000U
#define MT_LS31B 0x7FFFFFFFU
#define MT_A 2567483615U

struct MersenneTwister {
  unsigned int twistory[MT_N];
	int pos;
	MersenneTwister(unsigned int seed=0) {
    twistory[0]=seed;
    for (int i=1;i<MT_N;i++) twistory[i]=1812433253U*(twistory[i-1]^(twistory[i-1]>>30))+i;
    pos=0;
  }
  void generate() {
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
    unsigned int ans=twistory[pos++];
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

#undef MT_N
#undef MT_M
#undef MT_MSB
#undef MT_LS31B
#undef MT_A
#endif
