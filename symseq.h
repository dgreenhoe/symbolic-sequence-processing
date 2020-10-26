/*============================================================================
 * Daniel J. Greenhoe
 * header file for routines for char sequence functions
 *============================================================================*/
class symseq {
  private:
    long N;
    char *x;
  public:
    symseq(const long M);  //constructor initializing to '.'
    symseq(const long M,const unsigned seed,const char *symbols); //constructor initializing using seed
    void clear(void);                 //fill sequence with the value 'A'
    char get(const long n) const;           //get a value from x at location n
    char get(const long n, const char *symbols) const;   //get a value from x at location n but exit if not in <symbols> string
    const void put(const long n, const char symbol);
    void put(const long start, const long end, const char symbol);
    const long getN(void) const {const long M=(const long)N; return M;}        //get N
    void downsample(int M, symseq *y);//downsample by a factor of <M> and write to <y>
    void  list(const long start, const long end, const char* str1, const char *str2, const int display, FILE *fptr);
    void  list(const long start, const long end, FILE *fptr){list(start,end, "","",1,fptr);}
    void  list(const long start, const long end, const int display, FILE *fptr){list(start,end, "","",1,fptr);}
    void  list(const long start, const long end, const char* str1, const char *str2, FILE *fptr){
          list(           start,            end,             str1,             str2, 1,    fptr);
          }
    void  list(long start, long end){list(start,end,"","",NULL);} //list contents of sequence
    void  list(void){                list(0,    N-1,"","",NULL);}    //list contents of sequence
    void  list(long start){          list(start,N-1,"","",NULL);} //list contents of sequence
    void shiftL(long n);              //shift symseq n elements to the left
    void shiftR(long n);              //shift symseq n elements to the right
    void prngseed(unsigned seed){srand(seed);}    //
    void randomize(const char *symbols);//randomize using a list of symbols
    void randomize(const unsigned seed,const char *symbols){prngseed(seed); randomize(symbols);}
    void operator=(symseq y);     //x=y
    void  operator>>=(long n){shiftR(n);} //shift symseq n elements to the right
    void  operator<<=(long n){shiftL(n);} //shift symseq n elements to the left
    int operator==(symseq y);     //test if x==y; 1 if yes, 0 if no
  };

extern long cmp(const symseq *x,const symseq *y, int showdiff, FILE *fptr);
void copy(const long start, const long end, const symseq *x, symseq *y);
void downsample(int M, symseq *x, symseq *y);

