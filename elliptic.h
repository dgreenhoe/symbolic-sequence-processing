/*============================================================================
 * Daniel J. Greenhoe
 *============================================================================*/
/*=====================================
 * defines
 *=====================================*/
/*=====================================
 * classes
 *=====================================*/
//-----------------------------------------------------------------------------
//! \brief Ellipse with parameters (a,b), phi, and (xo,yo) where (x,y) is all the points in R^2
//! defined by 
//!   x^2    y^2   
//!   ---  + ---  = 1  
//!   a^2    b^2
//! rotated in R^2 by <phi> radians and offset by (xo+yo)
//------------------------------------------------------------------------------
class ellipsec{
  private:
    double a,b,phi,xo,yo;
  public:
    ellipsec( const double aa, const double bb ){ a=aa; b=bb; phi=0; xo=0; yo=0; }
    ellipsec( const double aa, const double bb, const double pp, const double xx, const double yy)
    {
      a   = aa;
      b   = bb;
      phi = pp;
      xo  = xx;
      yo  = yy;
    }
    ellipsec(void)
    { 
      a   = 1; 
      b   = 1; 
      phi = 0; 
      xo  = 0; 
      yo  = 0; 
    }
    void set( const double aa, const double bb, const double pp, const double xx, const double yy )
    { 
      a   = aa; 
      b   = bb; 
      phi = pp; 
      xo  = xx; 
      yo  = yy;
    }
    void   set(   const double u,  const double v){a=u; b=v;}
    void   set(   const double aa, const double bb, const double pp){a=aa; b=bb; phi=pp;} //set parameters
    void   seta(  const double aa){ a   = aa; }
    void   setb(  const double bb){ b   = bb; }
    void   setphi(const double pp){ phi = pp; }
    double geta(void){   return a  ; }
    double getb(void){   return b  ; }
    double getphi(void){ return phi; }
    double estimate(void);
    double x(const double theta); //get x component of a point at angle theta on ellipse(a,b)
    double y(const double theta); //get y component of a point at angle theta on ellipse(a,b)
    vectR2  xy(const double theta);
    double pathlength(const double ta, const double tb, const long N);
    double perimeter(const long N){ return 4*pathlength( 0, M_PI/2, N ); } 
    int    findt_dfroms(   const double s, const double targetd, const int direction, const long N, double *t, double *errord);
    int    setab_givenxyb( const double x, const double y, const double bb);
    int    setab_givenxya( const double x, const double y, const double aa);
    int    tgivenxy(       const double x, const double y, double *t);
    int    setab_givenxyb( const vectR2 p, const double bb );
    int    setab_givenxya( const vectR2 p, const double aa );
    int    tgivenxy(       const vectR2 p, double *t       );
    vectR2 normalize(      const vectR2 p                  );
  };

//typedef vectR2 (*fncttype)(double,double,double);
typedef vectR2 (*fncttype)(double);


extern double metric_balloon(vectR2 p, vectR2 q);


