/*============================================================================
 * Daniel J. Greenhoe
 *============================================================================*/
/*-------------------------------------------------------------------------
 * ordered pair (x,y)
 *-------------------------------------------------------------------------*/
class cfour {
  private:
    complex z[4];
  public:
    cfour(void);
    cfour(double ur, double ui);
    cfour   get(void)   const ;
    complex get (int n) const { return z[n]; }
    complex get1(void)  const { return z[0]; }
    complex get2(void)  const { return z[1]; }
    complex get3(void)  const { return z[2]; }
    complex get4(void)  const { return z[3]; }
    void puti(int i,complex uv){z[i]=uv;}
    void put(cfour uv);
    void put(double ur,double ui);
    void clear(void){put(0,0);}
    void   list(const char *str1, const char *str2) const;
    void   list(void)                               const { list("" , ""  ); }
    void   list(const char *str)                    const { list(str, "\n"); }
    void   listn(void)                              const { list("" , "\n"); }
  };

/*-------------------------------------------------------------------------
 *        | z1 |
 * vector | z2 | over C^4
 *        | z3 |
 *        | z4 |
 * vectors on C^4 are 4-tuples
 * (and hence inherit all the properties of class cfour)
 * but also have additional linear space (vector space) properties
 *-------------------------------------------------------------------------*/
class vectC4: public cfour {
  public:
    vectC4(double ur, double ui) : cfour(ur,ui){};
    vectC4(void) : cfour(){}; 
    vectC4  get(void)  const;
    complex get(int i) const { return cfour::get(i); }
    vectC4 conj(void);
    double mag(void);
    double norm(void){return mag();}
    double r(void)   {return mag();}
    void  operator+=(vectC4 q); 
    void  operator-=(vectC4 q); 
    void  operator*=(double a);
  };

class seqC4 {
  private:
    long N;
    vectC4 *x;
  public:
    seqC4(long M);               //constructor
    seqC4(long M, double u);     //constructor
    void fill(double u,double v);//fill seqC4 with the value (u,v)
    void clear(void){fill(0,0);}    //fill seqR1 with the value 0
    int  put(long n, vectC4 xyz);       //put a value <u> at location n in seq.
    int  put(long n, double u1,double u2,double u3,double u4);
    vectC4  get (long n) const { return x[n].get();  }  //get a value from x at location n
    complex get1(long n) const { return x[n].get1(); }  //get a value from x1 at location n
    complex get2(long n) const { return x[n].get2(); }  //get a value from x2 at location n
    complex get3(long n) const { return x[n].get3(); }  //get a value from x3 at location n
    complex get4(long n) const { return x[n].get4(); }  //get a value from x4 at location n
    long    getN(void)   const { return N; }      //get N
    void    list(const long start, const long end) const;//list contents of sequence
    void    list(void)   const {list(0,N-1);}   //list contents of sequence
    void    list1(const long start, const long end, const char *str1, const char *str2) const;//
    void    list1(const long start, const long end) const {list1(start,end,"","");}//
    void    list1(void)  const {list1(0,N-1);}  //list contents of seq. using 1 space each
    void    test(void);
  };

/*=====================================
 * operator overloading
 *=====================================*/
vectC4 operator-(vectC4 p);            // -p
vectC4 operator+(vectC4 p, vectC4 q);  // p+q
vectC4 operator-(vectC4 p, vectC4 q);  // p-q
vectC4 operator*(const double a, const vectC4 y);
vectC4 operator*(const complex z, const vectR4 x);
double operator^(vectC4 p,vectC4 q);   // "dot product" of p and q
seqC4  operator*(seqC4 xx,seqR1 y);

/*=====================================
 * functions
 *=====================================*/
extern int mag(seqC4 *xC4, seqR1 *ymag);




