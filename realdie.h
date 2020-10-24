/*============================================================================
 * Daniel J. Greenhoe
 * header file for routines for real die routines
 *   'A'--> die face value 1
 *   'B'--> die face value 2
 *   'C'--> die face value 3
 *   'D'--> die face value 4
 *   'E'--> die face value 5
 *   'F'--> die face value 6
 *============================================================================*/
class rdieseq: public dieseq {
  public:
    rdieseq(const long M) : dieseq(M){};
    rdieseq(const long M,const unsigned seed) : dieseq(M,seed){};
    void operator=(dieseq y);     //x=y
    void  operator>>=(long n){shiftR(n);} //shift rdieseq n elements to the right
    void  operator<<=(long n){shiftL(n);} //shift rdieseq n elements to the left
    int metrictbl(void);
    int Rxx(const seqR1 *Rxx, const int showcount);
    int Rxx(const seqR1 *Rxy, const int showcount,const long N, const long M, const long start, const long finish);
    int Rxxo(const seqR1 *rxx, const int showcount);
    double Rxx(const long m);
  };

extern rdieseq rdie_R1todie_euclid(seqR1 xyz);
extern rdieseq rdie_R3todie_larc(seqR3 xyz);
extern rdieseq rdie_R3todie0_larc(seqR3 xyz);
extern rdieseq rdie_R3todie_euclid(seqR3 xyz);
extern rdieseq rdie_R3todie0_euclid(seqR3 xyz);
extern vectR3 rdie_dietoR3(char c);
extern int     rdie_dietoR1(char c);
extern int     rdie_domain(char c);//check if value is in the domain of rdie
extern double  rdie_metric(char a, char b);
extern double  rdie_metric(rdieseq x, rdieseq y); //metric for two sequences

