/*============================================================================
 * Daniel J. Greenhoe
 * header file for routines for DNA routines
 *============================================================================*/
class dnaseq: public symseq {
  public:
    dnaseq(long M) : symseq(M){};  //constructor initializing to '.'
    void   seed(unsigned seed){srand(seed);}    //
    void   randomize(void){symseq::randomize("ATCG");}    //
    void   randomize(unsigned seed){srand(seed); randomize();}
    int    randomize(long start, long end, int wA,int wT,int wC,int wG);
    int    randomize(unsigned seed,int wA,int wT,int wC,int wG){srand(seed);return randomize(0,getN()-1,wA,wT,wC,wG);}
    int    randomize(int wA,int wT,int wC,int wG){return randomize(0,getN()-1,wA,wT,wC,wG);}
    int    randomize(long start, long end, unsigned seed,int wA,int wT,int wC,int wG){srand(seed);return randomize(start,end,wA,wT,wC,wG);}
    char   get(long n){return symseq::get(n,"ATCG");}     //get a value from x at location n
    void   put(long n, char symbol){symseq::put(n,symbol);}
    void   put(dnaseq *y, const long n, const char symbol);
    void   put(dnaseq *y, long n){return put(y,n,'.');}
    void   put(const long start, const long end, char c);//put a value <c> at locations start to end
    seqR2  dnatoR2(void);         //map gsp face values to R^2
    dnaseq downsample(int factor);//downsample by a factor of <factor>
    seqR1  dnatoR1(void);         //map dna sequence to R^1
    seqC1  dnatoC1(void);         //map dna sequence to R^1
    seqR1  dnatoR1pam(void);      //map dna sequence to R^1 using PAM scheme (symmetric about zero)
    seqR1  dnatoR1bin(void);    //map dna sequence to R^1 using AT-CG binary scheme
    seqR4  dnatoR4(void);         //map dna sequence to R^4 
    double Rxx (const long m);
    int    Rxx (const seqR1 *Rxx, const int showcount);
    int    Rxxo(const seqR1 *rxx, const int showcount);

    seqR1  histogram(const long start, const long end, int display, FILE *fptr);//compute, display, and write histogram 
    seqR1  histogram(){return histogram(0,getN()-1,0,NULL);}//compute histogram
    seqR1  histogram(const long start, const long end){return histogram(start,end,0,NULL);}//compute histogram
    seqR1  histogram(int display,FILE *fptr){return histogram(0,getN()-1,display,fptr);}//print histogram to file
    seqR1  histogram(FILE *fptr){return histogram(0,getN()-1,0,fptr);}//print histogram to file
    void operator=(dnaseq y);     //x=y
  };

extern long    numsym_fasta_file(const char *filename);
extern int     read_fasta_file(const char *filename, char *description, dnaseq *x);
extern int     dna_domain (char symbol);
extern vectR2  dna_dnatoR2(char symbol);
extern double  dna_dnatoR1(char symbol);
extern complex dnatoC1c   (char symbol);
extern vectR4  dnatoR4c   (char symbol);

