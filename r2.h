//============================================================================
// Daniel J. Greenhoe
//============================================================================
#include <vector>
//-------------------------------------------------------------------------
// ordered pair (x,y) on R^2
//-------------------------------------------------------------------------
class opair 
{
//  private:
  //double x,y;
  public:
    std::vector<double> pairxy = { 0, 0 };
    opair(const double u, const double v)        { pairxy.front() = u; pairxy.back() = v;                }
    opair(const double u)                        { pairxy.front() = u; pairxy.back() = u;                }
    opair(void)                                  { pairxy.front() = 0; pairxy.back() = 0;                }
    void   clear(void)                           { pairxy.front() = 0; pairxy.back() = 0;                }
    void   put(  const double u, const double v) { pairxy.front() = u; pairxy.back() = v;                }
    opair  get(  void)                                         const { opair p(getx(),gety()); return p; }
    double getx( void)                                         const { return pairxy.front();            }
    double gety( void)                                         const { return pairxy.back();             }
    void   list(const char* str1, const char *str2, FILE *ptr) const;
    void   list(const char* str1, const char *str2)            const { list(str1,str2,NULL);             }
    void   list(FILE *fptr)                                    const { list("","",fptr);                 }
    void   list(void)                                          const { list("","",NULL);                 }
};

//-----------------------------------------------------------------------------
// vector | x | over R^2
//        | y |
// vectors on R^2 are ordered pairs
// (and hence inherit all the properties of class opair)
// but also have additional linear space (vector space) properties
//-----------------------------------------------------------------------------
class vectR2: public opair
{
  public:
    vectR2(const double u, const double v) : opair(u,v){};
    vectR2(const double u) : opair(u){};
    vectR2(void) : opair(){};
    double theta(void)                                           const;
    inline double mag(void)                                      const { return sqrt(getx()*getx()+gety()*gety()); }
    inline double norm(void)                                     const { return mag();                             }
    inline void   polartoxy( const double r, const double theta)       { put(r*cos(theta) , r*sin(theta));         }
    inline void   add(       const double u, const double v    )       { put(getx() + u   , gety() + v );          }
    void   operator+=(const vectR2 q)                                  { put(getx()+q.getx(), gety()+q.gety());    }
    void   operator-=(const vectR2 q)                                  { put(getx()-q.getx(), gety()-q.gety());    }
    void   operator&=(const double phi);
    void   operator*=(const double a)                                  { put(a*getx()     , a*gety()   );          }
    vectR2 operator* (const double a)                            const { vectR2 u(a*getx(), a*gety()   ); return u;}
};

vectR2  operator-( const vectR2 p                   );
vectR2  operator+( const vectR2 p, const vectR2 q   );
vectR2  operator-( const vectR2 p, const vectR2 q   );
vectR2  operator&( const vectR2 p, const double phi );
inline double  operator^( const vectR2 p, const vectR2 q ) {return p.getx()*q.getx() + p.gety()*q.gety(); }

//-------------------------------------------------------------------------
// class of sequences over R^2
//-------------------------------------------------------------------------
class seqR2 
{
  private:
    long N;
    vectR2 *xy;
  public:
    seqR2(long M);
    seqR2(long M, double u);
    void   clear(void);
    void   fill(double u);
    void   inc( double x0, double y0,double dx, double dy);
    int    put( long n,    vectR2 xy);
    int    put( long n,    double u, double v);
    vectR2 get(  const long n)                    const ;
    double getx( const long n)                    const ;
    double gety( const long n)                    const ;
    long   getN(void)                             const {return N;}
    double norm( const long n)                    const ;
    vectR2 max(const int verbose)                 const ;
    vectR2 max(void)                              const {return max(0);}
    void   list(const long start, const long end) const {list(start,end,"","",stdout);}
    void   list(void)                             const {list(0,N-1,"","",stdout);}
    void   list1(void)                            const ;
    void   list1(long start, long end)            const ;
    void   test(void)                             const ;
    void list(const long start, const long end, const char *str1, const char *str2, FILE *ptr) const;
    void list(const long start, const long end, const char *str1, const char *str2, int display, FILE *fptr) const {
         if(display) list(start,end,str1,str2,stdout);
         if(1)       list(start,end,str1,str2,fptr);
         }
    void list(const char* str1, const char *str2,int display,FILE *fptr) const {
         if(display) list(0,N-1,str1,str2,stdout);
         if(1)       list(0,N-1,str1,str2,fptr);
         }
};

//=====================================
// functions
//=====================================
extern double pqtheta(const vectR2 p, const vectR2 q);
inline double chordlength(const vectR2 p, const vectR2 q)
{
  const vectR2  pqd = p - q; // differnce of p and q
  return pqd.norm();         // "length" of difference
}