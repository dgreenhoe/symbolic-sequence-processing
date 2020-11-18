/*============================================================================
 * Daniel J. Greenhoe
 *============================================================================*/
/*-------------------------------------------------------------------------
 * ordered pair (x,y)
 *-------------------------------------------------------------------------*/
class csix {
  private:
    complex z[6];
  public:
    csix(void);
    csix(double ur, double ui);
    csix    get(void)   const;
    complex get (int n) const { return z[n]; }
    complex get1(void)  const { return z[0]; }
    complex get2(void)  const { return z[1]; }
    complex get3(void)  const { return z[2]; }
    complex get4(void)  const { return z[3]; }
    complex get5(void)  const { return z[4]; }
    complex get6(void)  const { return z[5]; }
    void puti(int i,complex uv){z[i]=uv;}
    void put(csix uv);
    void put(double ur,double ui);
    void clear(void){put(0,0);}
    void list(const char *str1, const char *str2) const;
    void list(void)            const {list("","");}
    void list(const char *str) const {list(str,"\n");}
    void listn(void)           const {list("","\n");}
  };

/*-------------------------------------------------------------------------
 *        | z1 |
 *        | z2 |
 * vector | z3 | over C^6
 *        | z4 |
 *        | z5 |
 *        | z6 |
 * vectors on C^6 are 6-tuples
 * (and hence inherit all the properties of class csix)
 * but also have additional linear space (vector space) properties
 *-------------------------------------------------------------------------*/
class vectC6: public csix {
  public:
    vectC6(double ur, double ui) : csix(ur,ui){};
    vectC6(void) : csix(){}; 
    vectC6 get(void)   const;
    complex get(int i) const { return csix::get(i); }
    vectC6 conj(void)  const;
    double mag(void)   const;
    double norm(void)  const { return mag(); }
    double r(void)     const { return mag(); }
    void  operator += (vectC6 q); //{put(get1()+q.get1(),get2()+q.get2(),get3()+q.get3(),get4()+q.get4(),get5()+q.get5(),get6()+q.get6());} //p=p+q
    void  operator -= (vectC6 q); //{put(get1()-q.get1(),get2()-q.get2(),get3()-q.get3(),get4()-q.get4(),get5()-q.get5(),get6()-q.get6());} //p=p-q
    void  operator *= (double a);
  };

class seqC6 {
  private:
    long N;
    vectC6 *x;
  public:
    seqC6(long M);               //constructor
    seqC6(long M, double u);     //constructor
    void fill(double u,double v);//fill seqC6 with the value (u,v)
    void clear(void){fill(0,0);}    //fill seqR1 with the value 0
    int  put(long n, vectC6 xyz);       //put a value <u> at location n in seq.
    int  put(long n, double u1,double u2,double u3,double u4,double u5,double u6);
    vectC6 get  (long n){return x[n].get();}   //get a value from x at location n
    complex get1(long n){return x[n].get1();}  //get a value from x1 at location n
    complex get2(long n){return x[n].get2();}  //get a value from x2 at location n
    complex get3(long n){return x[n].get3();}  //get a value from x3 at location n
    complex get4(long n){return x[n].get4();}  //get a value from x4 at location n
    complex get5(long n){return x[n].get5();}  //get a value from x5 at location n
    complex get6(long n){return x[n].get6();}  //get a value from x6 at location n
    long getN(void){return N;}      //get N
    void list(const long start, const long end, const char *str1, const char *str2, FILE *ptr) const;
    void list(const long start, const long end, const char *str1, const char *str2, const int display, FILE *fptr) const
         {
           if(display)   list(start,end,str1,str2,stdout);
           if(1)         list(start,end,str1,str2,fptr);
         }
    void list(const char* str1, const char *str2, const int display, FILE *fptr) const
         {
           if(display)   list(0,N-1,str1,str2,stdout);
           if(1)         list(0,N-1,str1,str2,fptr);
         }
    void list(const long start, const long end)                    const {list(start,end,"","",stdout);}
    void list(void)                                                const {list(0,N-1,"","",stdout);}
    void list1(const long start, const long end, const char *str1, const char *str2)const ;//
    void list1(const long start, const long end)                   const {list1(start,end,"","");}//
    void list1(void)                                               const {list1(0,N-1);}  //list contents of seq. using 1 space each
    void test(void);
  };

/*=====================================
 * operator overloading
 *=====================================*/
vectC6 operator-(vectC6 p);            // -p
vectC6 operator+(vectC6 p, vectC6 q);  // p+q
vectC6 operator-(vectC6 p, vectC6 q);  // p-q
vectC6 operator*(const double a, const vectC6 y);
vectC6 operator*(const complex z, const vectR6 x);
//vectC6 operator*(const complex z, const vectR6 x){vectC6 y(
//z*(x.get1()),
//z*(x.get2()),
//z*(x.get3()),
//z*(x.get4()),
//z*(x.get5()),
//z*(x.get6())); 
//return y;}
double operator^(vectC6 p,vectC6 q);   // "dot product" of p and q
seqC6  operator*(seqC6 xx,seqR1 y);

/*=====================================
 * functions
 *=====================================*/
//extern double pqtheta(vectR6 p, vectR6 q); //return radians between vectors induced by p and q in R^6
extern int mag(seqC6 *xC6, seqR1 *ymag);




