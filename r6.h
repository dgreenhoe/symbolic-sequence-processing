//=============================================================================
//! Daniel J. Greenhoe
//=============================================================================
#include <vector>
//#include <Eigen/Dense>  // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=89325

//typedef Eigen::Matrix< double, 6, 1 > Vector6d;

//-----------------------------------------------------------------------------
//! \brief Ordered pair (x,y)
//-----------------------------------------------------------------------------
class osix 
{
  private:
    std::vector<double> x = {0, 0, 0, 0, 0, 0};
  public:
    osix(void);
    osix(double u0, double u1, double u2, double u3, double u4, double u5);
    osix(double u);
    const double* getdata(void) const { return x.data(); }
    double* getdataa(void){ return x.data(); }
    osix   get(void) const ;
    double get(int n) const {return x.at(n)  ;}
    double get1(void) const {return x.front();}; //get component x1
    double get2(void) const {return x.at(1)  ;}; //get component x2
    double get3(void) const {return x.at(2)  ;}; //get component x3
    double get4(void) const {return x.at(3)  ;}; //get component x4
    double get5(void) const {return x.at(4)  ;}; //get component x5
    double get6(void) const {return x.back() ;}; //get component x6
    double max(void)  const;
    double min(void)  const;
    void   put(double u0, double u1, double u2, double u3, double u4, double u5)
    {
      x.front() = u0;
      x.at(1)   = u1;
      x.at(2)   = u2;
      x.at(3)   = u3;
      x.at(4)   = u4;
      x.back()  = u5;
    }
    void put(int n,double u){ x.at(n) = u; }
    void put(double u);
    void put(osix u);
    void clear(void){ put(0); }
    void list(const char *str1, const char *str2);
    void list(void){list("","");}
    void list(const char *str){list(str,"\n");}
    void listn(void){list("","\n");}
};

//-----------------------------------------------------------------------------
//!        | x1 |
//!        | x2 |
//! vector | x3 | over R^6
//!        | x4 |
//!        | x5 |
//!        | x6 |
//! vectors on R^6 are 6-tuples
//! (and hence inherit all the properties of class osix)
//! but also have additional linear space (vector space) properties
//-----------------------------------------------------------------------------
class vectR6: public osix 
{
  public:
    vectR6(double u1,double u2,double u3,double u4,double u5,double u6) : osix(u1,u2,u3,u4,u5,u6){};
    vectR6(double u) : osix(u){};
    vectR6(void) : osix(){}; 
    vectR6 get(void) const;
    const double get(int i) const {return osix::get(i);}
    const double mag(void) const;
    double norm(void)const {return mag();}
    double r(void)   const {return mag();}
    vectR6 mpy       ( const double a );
    void   operator+=( const vectR6 q );
    void   operator-=( const vectR6 q );
    void   operator*=( const double a );
};

class seqR6 
{
  private:
    long N;
    vectR6 *x;
  public:
    seqR6(const long M);
    seqR6(const long M, const double u);
    void   fill(const double u);
    void   clear(void);
    int    put( const long n, const vectR6 xyz);
    int    put( const long n, double u1,double u2,double u3,double u4,double u5,double u6);
    vectR6 get (const long n) const { return x[n].get() ; }  //! \brief get a value from x at location n
    double get1(const long n) const { return x[n].get1(); }  //! \brief get a value from x1 at location n
    double get2(const long n) const { return x[n].get2(); }  //! \brief get a value from x2 at location n
    double get3(const long n) const { return x[n].get3(); }  //! \brief get a value from x3 at location n
    double get4(const long n) const { return x[n].get4(); }  //! \brief get a value from x4 at location n
    double get5(const long n) const { return x[n].get5(); }  //! \brief get a value from x5 at location n
    double get6(const long n) const { return x[n].get6(); }  //! \brief get a value from x6 at location n
    long   getN(void)         const { return N;           }  //! \brief get N
    void list(const long start, const long end, const char *str1, const char *str2, FILE *ptr);
    void list(const long start, const long end, const char *str1, const char *str2, int display, FILE *fptr){
         if(display) list(start,end,str1,str2,stdout);
         if(1)       list(start,end,str1,str2,fptr);
         }
    void list(const char* str1, const char *str2,int display,FILE *fptr){
         if(display) list(0,N-1,str1,str2,stdout);
         if(1)       list(0,N-1,str1,str2,fptr);
         }
    void list(const long start, const long end){list(start,end,"","",stdout);}
    void list(void){list(0,N-1,"","",stdout);}
    void list1(const long start, const long end, const char *str1, const char *str2, FILE *ptr);
    void list1(const long start, const long end, const char *str1, const char *str2, int display, FILE *fptr){
         if(display) list1(start,end,str1,str2,stdout);
         if(1)       list1(start,end,str1,str2,fptr);
         }
    void list1(const char* str1, const char *str2,int display,FILE *fptr){
         if(display) list1(0,N-1,str1,str2,stdout);
         if(1)       list1(0,N-1,str1,str2,fptr);
         }
    void list1(const long start, const long end){list1(start,end,"","",stdout);}
    void list1(void){list1(0,N-1,"","",stdout);}
    void test(void);
};

/*=====================================
 * operator overloading
 *=====================================*/
vectR6 operator-( const vectR6 p);            // -p
vectR6 operator+( const vectR6 p, const vectR6 q);  // p+q
vectR6 operator-( const vectR6 p, const vectR6 q);  // p-q
vectR6 operator*( const double a, const vectR6 y);
inline vectR6 operator*( const vectR6 y, const double a){return a*y;}
double operator^( const vectR6 p, const vectR6 q);   // "dot product" of p and q

/*=====================================
 * functions
 *=====================================*/
double pqtheta(const vectR6 p, const vectR6 q); //return radians between vectors induced by p and q in R^6




