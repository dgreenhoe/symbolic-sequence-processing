/*============================================================================
 * Daniel J. Greenhoe
 * header file for routines for spinner routines
 *   'A'--> spin face value 1
 *   'B'--> spin face value 2
 *   'C'--> spin face value 3
 *   'D'--> spin face value 4
 *   'E'--> spin face value 5
 *   'F'--> spin face value 6
 *============================================================================*/
class spinseq: public dieseq {
  public:
    spinseq(const long M) : dieseq(M){};
    spinseq(const long M,const unsigned seed) : dieseq(M,seed){};
    seqR1   spintoR1(void);         //map spin face values to R^1
    seqR2   spintoR2(void);         //map spin face values to R^2
    void    operator=(spinseq y);     //x=y
    int     metrictbl(void);
    double  Rxx (const long m);
    int     Rxx (seqR1 *Rxx, const int showcount);
    int     Rxxo(seqR1 *rxx, const int showcount);
    spinseq downsample(const int factor);//downsample by a factor of <factor>
  };

extern int spin_domain(char c);//check if value is in the domain of rspin
extern spinseq spin_R1tospin_euclid( const seqR1 xy);
extern spinseq spin_R2tospin_larc(   const seqR2 xy);
extern spinseq spin_R2tospin0_larc(  const seqR2 xy);
extern spinseq spin_R2tospin_euclid( const seqR2 xy);
extern spinseq spin_R2tospin0_euclid(const seqR2 xy);
extern vectR2  spin_spintoR2(const char c);
extern double  spin_spintoR1(const char c);
extern double spin_metric(const char a, const char b);
extern double spin_metric(const spinseq x, const spinseq y); //metric for two sequences
//extern seqR1 spin_correlation(spinseq x, spinseq y, int showcount);//correlation
//extern seqR1 spin_correlation(spinseq x, spinseq y){return spin_correlation(x,y,0);}//correlation


