/*============================================================================
 * Daniel J. Greenhoe
 *============================================================================*/
/*-------------------------------------------------------------------------
 * ordered pair (x,y) on R^2
 *-------------------------------------------------------------------------*/
class opair {
  private:
    double x,y;
  public:
    opair(const double u, const double v)            { x = u; y = v;                     }
    opair(const double u)                            { x = u; y = u;                     }
    opair(void)                                      { x = 0; y = 0;                     }
    void   clear(void)                               { x = 0; y = 0;                     }
    void   put(  const double u, const double v)     { x = u; y = v;                     }
    opair  get(  void)                         const { opair p(getx(),gety()); return p; }
    double getx( void)                         const { return x;                         }
    double gety( void)                         const { return y;                         }
    void  list(const char* str1, const char *str2, FILE *ptr);
    void  list(const char* str1, const char *str2){list(str1,str2,NULL);}
    void  list(FILE *fptr){list("","",fptr);} //list contents of sequence
    void  list(void){list("","",NULL);}    //list contents of sequence
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
    void   operator+=(const vectR2 q)                                  { put(getx()+q.getx(), gety()+q.gety());    }
    void   operator-=(const vectR2 q)                                  { put(getx()-q.getx(), gety()-q.gety());    }
    void   operator&=(const double phi);
    inline double mag(void)                                      const { return sqrt(getx()*getx()+gety()*gety()); }
    inline double norm(void)                                     const { return mag();                             }
    inline void   polartoxy( const double r, const double theta)       { put(r*cos(theta) , r*sin(theta));         }
    inline void   add(       const double u, const double v    )       { put(getx() + u   , gety() + v );          }
    inline void   operator*=(const double a)                           { put(a*getx()     , a*gety()   );          }
    inline vectR2 operator* (const double a)                     const { vectR2 u(a*getx(), a*gety()   ); return u;}
};

vectR2  operator-( const vectR2 p                   );
vectR2  operator+( const vectR2 p, const vectR2 q   );
vectR2  operator-( const vectR2 p, const vectR2 q   );
vectR2  operator&( const vectR2 p, const double phi );
inline double  operator^(const vectR2 p, const vectR2 q){return p.getx()*q.getx() + p.gety()*q.gety();}

/*-------------------------------------------------------------------------
 * class of sequences over R^2
 *-------------------------------------------------------------------------*/
class seqR2 {
  private:
    long N;
    vectR2 *xy;
  public:
    seqR2(long M);             //constructor
    seqR2(long M, double u);   //constructor
    void   clear(void);              //fill seqR1 with the value 0
    void   fill(double u);           //fill seqR1 with the value <u>
    void   inc(double x0, double y0,double dx, double dy);
    int    put(long n, vectR2 xy);   //put a value <u> at location n in seq.
    int    put(long n, double u, double v);
    vectR2 get(const long n)  const ;              //get a value from x at location n
    double getx(const long n) const ;             //get a value from x at location n
    double gety(const long n) const ;             //get a value from x at location n
    long   getN(void)         const {return N;}     //get N
    double norm(const long n) const ;             //norm of element x_n
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
    void   list1(void);               //list contents of seq. using 1 space each
    void   list1(long start, long end);//
    void   test(void);
    vectR2 max(const int verbose) const;
    vectR2 max(void) const {return max(0);}   //max mode=0=no print        
  };

/*=====================================
 * functions
 *=====================================*/
extern double pqtheta(const vectR2 p, const vectR2 q); //return radians between vectors induced by p and q in R^2
inline double chordlength(vectR2 p, vectR2 q){
  vectR2  pqd=p-q; // differnce of p and q
  return pqd.norm(); // "length" of difference
  }