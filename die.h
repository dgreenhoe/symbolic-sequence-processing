/*============================================================================
 * Daniel J. Greenhoe
 * header file for routines for die routines
 *   'A'--> die face value 1
 *   'B'--> die face value 2
 *   'C'--> die face value 3
 *   'D'--> die face value 4
 *   'E'--> die face value 5
 *   'F'--> die face value 6
 *============================================================================*/
class dieseq: public symseq {
  public:
    dieseq(const long M) : symseq(M){};//constructor initializing to '.'
    dieseq(const long M,const unsigned seed) : symseq(M,seed,"ABCDEF"){};//constructor initializing random values
    void  randomize(void){symseq::randomize("ABCDEF");}    //
    void  randomize(unsigned seed){srand(seed); randomize();}
    int   randomize(long start, long end, int wA,int wB,int wC,int wD,int wE,int wF);
    int   randomize(unsigned seed,int wA,int wB,int wC,int wD,int wE,int wF){srand(seed);return randomize(0,getN()-1,wA,wB,wC,wD,wE,wF);}
    int   randomize(int wA,int wB,int wC,int wD,int wE,int wF){return randomize(0,getN()-1,wA,wB,wC,wD,wE,wF);}
    int   randomize(long start, long end, unsigned seed, int wA,int wB,int wC,int wD,int wE,int wF){srand(seed);return randomize(start,end, wA,wB,wC,wD,wE,wF);}
    char  get(long n){return symseq::get(n,"ABCDEF");}     //get a value from x at location n
    void  put(long n, char symbol){symseq::put(n,symbol);}
    seqR1 dietoR1(void);         //map die face values to R^1
    seqC1 dietoC1(void);         //map die face values to R^1
    seqR1 dietoR1pam(void);      //map die face values to R^1 using PAM scheme (symmetric about zero)
    seqR3 dietoR3(void);         //map die face values to R^3
    seqR1 histogram(const long start, const long end, int display, FILE *fptr);//compute, display, and write histogram 
    seqR1 histogram(){return histogram(0,getN()-1,0,NULL);}//compute histogram
    seqR1 histogram(const long start, const long end){return histogram(start,end,0,NULL);}//compute histogram
    seqR1 histogram(int display,FILE *fptr){return histogram(0,getN()-1,1,fptr);}//print histogram to file
    seqR1 histogram(FILE *fptr){return histogram(0,getN()-1,0,fptr);}//print histogram to file
    void operator=(dieseq y);     //x=y
  };

extern int     die_domain(char c);//check if value is in the domain of die
extern double  die_dietoR1 (char c);
extern double  die_dietoR1pam(char c);
extern vectR3  die_dietoR3 (char c);
extern vectR6  die_dietoR6 (char c);
extern complex die_dietoC1c(char c);

