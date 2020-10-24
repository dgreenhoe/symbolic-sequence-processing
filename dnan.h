/*============================================================================
 * Daniel J. Greenhoe
 * header file for routines for DNA routines with symbols
 * A, T, C, G, and N
 *============================================================================*/
class dnanseq: public symseq {
  public:
    dnanseq(long M) : symseq(M){};  //constructor initializing to '.'
    void   seed(unsigned seed){srand(seed);}    //
    void   randomize(void){symseq::randomize("ATCGN");}    //
    void   randomize(unsigned seed){srand(seed); randomize();}
    int    randomize(long start, long end, int wA,int wT,int wC,int wG,int wN);
    int    randomize(unsigned seed,int wA,int wT,int wC,int wG,int wN){srand(seed);return randomize(0,getN()-1,wA,wT,wC,wG,wN);}
    int    randomize(int wA,int wT,int wC,int wG, int wN){return randomize(0,getN()-1,wA,wT,wC,wG,wN);}
    int    randomize(long start, long end, unsigned seed,int wA,int wT,int wC,int wG,int wN){srand(seed);return randomize(start,end,wA,wT,wC,wG,wN);}
    char   get(long n){return symseq::get(n,"ATCGN");}     //get a value from x at location n
    void   put(long n, char symbol){symseq::put(n,symbol);}
    void   put(dnanseq *y, const long n, const char symbol);
    void   put(dnanseq *y, long n){return put(y,n,'.');}
    void   put(const long start, const long end, char c);//put a value <c> at locations start to end
    seqR2  dnantoR2(void);         //map gsp face values to R^2
    dnanseq downsample(int factor);//downsample by a factor of <factor>
    seqR1  dnantoR1(void);         //map dnan sequence to R^1
    seqC1  dnantoC1(void);         //map dnan sequence to R^1
    seqR1  dnantoR1pam(void);      //map dnan sequence to R^1 using PAM scheme (symmetric about zero)
    seqR1  dnantoR1bin(void);    //map dnan sequence to R^1 using AT-CG binary scheme
    seqR4  dnantoR4(void);         //map dnan sequence to R^4 
    double Rxx (const long m);
    int    Rxx (const seqR1 *Rxx, const int showcount);
    int    Rxxo(const seqR1 *rxx, const int showcount);
    seqR1  histogram(const long start, const long end, int display, FILE *fptr);//compute, display, and write histogram 
    seqR1  histogram(){return histogram(0,getN()-1,0,NULL);}//compute histogram
    seqR1  histogram(const long start, const long end){return histogram(start,end,0,NULL);}//compute histogram
    seqR1  histogram(int display,FILE *fptr){return histogram(0,getN()-1,1,fptr);}//print histogram to file
    seqR1  histogram(FILE *fptr){return histogram(0,getN()-1,0,fptr);}//print histogram to file
    void operator=(dnanseq y);     //x=y
  };

extern long numsym_fasta_file(const char *filename);
extern int read_fasta_file(const char *filename, char *description, const dnanseq *x);
extern int dnan_domain(const char c);
extern vectR2 dnan_dnatoR2(const char c);
extern double dnan_dnatoR1(const char c);
extern seqR1 dnan_correlation(dnanseq x, dnanseq y, int showcount);//correlation
seqR1 dnan_correlation(dnanseq x, dnanseq y){return dnan_correlation(x,y,0);}//correlation
extern double dnan_metric(const char a, const char b);

