//=============================================================================
//! \brief Daniel J. Greenhoe
//=============================================================================
#include <vector>
//-----------------------------------------------------------------------------
//! \brief ordered triple (x,y,y)
//-----------------------------------------------------------------------------
class otriple 
{
  private:
    std::vector<double> xyz = { 0, 0, 0 };
  public:
    otriple(double u, double v, double w){xyz.front()=u; xyz.at(1)=v; xyz.back()=w;}
    otriple(double u){ xyz.front()=u; xyz.at(1)=u; xyz.back()=u; }
    otriple(void)    { xyz.front()=0; xyz.at(1)=0; xyz.back()=0; }
    void   put(double u,double v,double w){xyz.front()=u; xyz.at(1)=v; xyz.back()=w; }
    //void   put(int u, int v, int w){ put((double)u, (double)v, (double)w); }
    void   clear(void) { xyz.front()=0; xyz.at(1)=0; xyz.back()=0; }
    const double *getdata(void) const { return xyz.data(); }
    double *getdataa(void) { return xyz.data(); }
    void   list(){ printf("(x=%12.8lf, y=%12.8lf, z=%12.8lf)\n", getx(), gety(), getz() ); }
    double getx(void) const {return xyz.front(); }; //get component x
    double gety(void) const {return xyz.at(1)  ; }; //get component y
    double getz(void) const {return xyz.back() ; }; //get component z
    double max(void) const;
    double min(void) const;
};

//-------------------------------------------------------------------------
//        | x |
// vector | y | over R^3
//        | z |
// vectors on R^3 are ordered pairs 
// (and hence inherit all the properties of class opair)
// but also have additional linear space (vector space) properties
//-------------------------------------------------------------------------*/
class vectR3: public otriple 
{
  public:
    vectR3(double u, double v, double w) : otriple(u,v,w){};
    vectR3(double u) : otriple(u){};
    vectR3(void) : otriple(){};
    void   put(vectR3 abc){ otriple::put(abc.getx(), abc.gety(), abc.getz()); }
    void   put(   int u,    int v,    int w){ otriple::put((double)u, (double)v, (double)w); }
    void   put(   int u,    int v, double w){ otriple::put((double)u, (double)v, (double)w); }
    void   put(   int u, double v,    int w){ otriple::put((double)u, (double)v, (double)w); }
    void   put(   int u, double v, double w){ otriple::put((double)u, (double)v, (double)w); }
    void   put(double u,    int v,    int w){ otriple::put((double)u, (double)v, (double)w); }
    void   put(double u,    int v, double w){ otriple::put((double)u, (double)v, (double)w); }
    void   put(double u, double v,    int w){ otriple::put((double)u, (double)v, (double)w); }
    void   put(double u, double v, double w){ otriple::put((double)u, (double)v, (double)w); }
    double mag(void)  const;
    double norm(void) const { return mag(); }
    void polartoxyz(const double r, const double theta, const double phi);
    void  operator+=(const vectR3 q);
    void  operator-=(const vectR3 q);
};

class seqR3 {
  private:
    long N;
    double *x,*y,*z;
    vectR3 *seqr3;
  public:
    seqR3(long M);               //constructor
    seqR3(long M, double u);     //constructor
    void clear(void);               //fill seqR1 with the value 0
    void fill(double u);            //fill seqR1 with the value <u>
    int  put(long n, vectR3 xyz);  //put a value <u> at location n in seq.
    int  put(long n, double u, double v, double w);
    void randomize(unsigned seed);   //
    void randomize(unsigned seed, int min, int max); //
    vectR3 get(long n);              //get a value from x at location n
    double getx(long n);             //get a value from x at location n
    double gety(long n);             //get a value from x at location n
    double getz(long n);             //get a value from x at location n
    long getN(void){return N;}      //get N
    void list(const long start, const long end, const char *str1, const char *str2, FILE *ptr);
    void list(const long start, const long end, const char *str1, const char *str2, int display, FILE *fptr){
         if(display) list(start,end,str1,str2,stdout);
         if(1)       list(start,end,str1,str2,fptr);
         }
    void list(const long start, const long end){list(start,end,"","",stdout);}
    void list(void){list(0,N-1,"","",stdout);}
    void list1(void);              //list contents of seq. using 1 space each
    void list1(long start, long end);//
    void test(void);
  };

//=====================================
// operator overloading
//=====================================
vectR3  operator-(vectR3 p);               // -p
vectR3  operator+(vectR3 p, vectR3 q);  // p+q
vectR3  operator-(vectR3 p, vectR3 q);  // p-q
vectR3  operator&(vectR3 p,double phi);    // <p> rotated counter-clockwise by <phi>
inline vectR3  operator*(vectR3 p,double a){vectR3 pa(p.getx()*a,p.gety()*a,p.getz()*a); return pa;}
inline double  operator^(vectR3 p,vectR3 q){return p.getx()*q.getx() + p.gety()*q.gety() + p.getz()*q.getz();}   // "dot product" of p and q
//double  operator*(vectR3 p,vectR3 q){return p.getx()*q.getx() + p.gety()*q.gety() + p.getz()*q.getz();}   // "dot product" of p and q

//=====================================
// functions
//=====================================
extern double pqtheta(const vectR3 p, const vectR3 q); //return radians between vectors induced by p and q in R^3




