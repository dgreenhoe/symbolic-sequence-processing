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
class fdieseq: public dieseq {
  public:
    fdieseq(const long M) : dieseq(M){};
    fdieseq(const long M,const unsigned seed) : dieseq(M,seed){};
    void operator=(dieseq y);     //x=y
    int metrictbl(void);
    double Rxx(const long m);
    int Rxx(const seqR1 *Rxx, const int showcount);
    int Rxx(const seqR1 *Rxy, const int showcount,const long N, const long M, const long start, const long finish);
    int Rxxo(const seqR1 *rxx, const int showcount);
    seqR6 dietoR6(void);         //map die face values to R^6
  };

extern fdieseq fdie_R1todie_euclid(seqR1 xyz);
extern fdieseq fdie_R6todie_ae(seqR6 xyz);
extern fdieseq fdie_R6todie_larc(seqR6 xyz);
extern fdieseq fdie_R6todie0_ae(seqR6 xyz);
extern fdieseq fdie_R6todie_euclid(seqR6 xyz);
extern fdieseq fdie_R6todie0_euclid(seqR6 xyz);
extern double  fdie_metric(char a, char b);
extern double  fdie_metric(fdieseq x, fdieseq y); //metric for two sequences

