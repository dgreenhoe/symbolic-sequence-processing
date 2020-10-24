/*============================================================================
 * Daniel J. Greenhoe
 *============================================================================*/
/*=====================================
 * defines
 *=====================================*/
/*=====================================
 * classes
 *=====================================*/
/*-------------------------------------------------------------------------
 * ellipse with parameters (a,b), phi, and (xo,yo) where (x,y) is all the points in R^2
 * defined by 
 *   x^2    y^2   
 *   ---  + ---  = 1  
 *   a^2    b^2
 * rotated in R^2 by <phi> radians and offset by (xo+yo)
 *-------------------------------------------------------------------------*/
class ellipsec{
  private:
    double a,b,phi,xo,yo;
  public:
    ellipsec(double aa, double bb){a=aa;b=bb;phi=0;xo=0;yo=0;}           //constructor using 2 long float arguments
    ellipsec(double aa, double bb, double pp, double xx, double yy){a=aa;b=bb;phi=pp;xo=xx;yo=yy;} //constructor
    ellipsec(void){a=1; b=1;;phi=0;xo=0;yo=0;}                        //constructor using no arguments (set to 0,0)
    void   set(double aa, double bb, double pp, double xx, double yy){a=aa;b=bb;phi=pp;xo=xx;yo=yy;} //constructor
    void   set(double u,double v){a=u; b=v;}         //set (a,b)=(u,v)
    void   set(double aa,double bb, double pp){a=aa; b=bb; phi=pp;} //set parameters
    void   seta(double aa){a=aa;}                    //a = seta(aa)
    void   setb(double bb){b=bb;}                    //b = setb(bb)
    void   setphi(double pp){phi=pp;}                //phi = setphi(phi)
    double geta(void){return a;}                     //get a parameter of ellipse(a,b)
    double getb(void){return b;}                     //get b parameter of ellipse(a,b)
    double getphi(void){return phi;}                     //get b parameter of ellipse(a,b)
    double x(double t); //{return a*cos(phi)*cos(t)-b*sin(phi)*sin(t)+xo;}  //get x component of a point at angle theta on ellipse(a,b)
    double y(double t); //{return a*sin(phi)*cos(t)+b*cos(phi)*sin(t)+yo;}  //get y component of a point at angle theta on ellipse(a,b)
    vectR2  xy(double t);
    double pathlength(double ta, double tb, long N);
    double perimeter(long N){return 4*pathlength(0,PI/2, N);} 
    int    findt_dfroms(double s, double targetd, int direction, long N, double *t, double *errord); //find parameter t that is a distance <targetd> from parameter <s>
    int    setab_givenxyb(double x, double y, double bb);
    int    setab_givenxya(double x, double y, double aa);
    int    setab_givenxyb(vectR2 p, double bb);
    int    setab_givenxya(vectR2 p, double aa);
    int    tgivenxy(vectR2 p, double *t);             //get parameter t givein point p=(x,y)
    int    tgivenxy(double x, double y, double *t);  //get parameter t givein (x,y)
    vectR2  normalize(vectR2 p);
    double estimate(void);
  };

//typedef vectR2 (*fncttype)(double,double,double);
typedef vectR2 (*fncttype)(double);


extern double metric_balloon(vectR2 p, vectR2 q);


