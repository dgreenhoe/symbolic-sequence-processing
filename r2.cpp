//=============================================================================
// Daniel J. Greenhoe
// normed linear space R^2
//=============================================================================
//=====================================
// headers
//=====================================
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <complex>
#include <Eigen/Dense>
#include "r1.h"
#include "r2.h"
#include "r2op.h"

//-----------------------------------------------------------------------------
//! List value of opair
//-----------------------------------------------------------------------------
void opair::list(const char *str1, const char *str2, FILE *ptr) const
{
  if( strlen(str1)>0 ) fprintf(ptr,"%s",str1);
  fprintf(ptr,"(%9.6lf,%9.6lf)",getx(),gety());
  if( strlen(str2)>0 ) fprintf(ptr,"%s",str2);
}

//-----------------------------------------------------------------------------
//! \brief polar rotation coordinate <theta> of opair point (x,y)
//! \returns Return value is in the half open interval [0:2pi), 
//!          -1 on error
//-----------------------------------------------------------------------------
double vectR2::theta(void) const
{
  const double x = getx();
  const double y = gety();
  if(x==0){
    if(     y == 0 ) return -1           ;
    else if(y > 0  ) return  1 * M_PI / 2;
    else             return  3 * M_PI / 2;
  }
  if(y==0){
    if(     x>0    ) return 0;
    else             return M_PI;
  }
  if( x>0 && y>0 ) return        atan( y/x); // 1st quadrant
  if( x>0 && y<0 ) return 2*M_PI-atan(-y/x); // 4th quadrant
  if( x<0 && y>0 ) return   M_PI-atan(-y/x); // 2nd quadrant
  if( x<0 && y<0 ) return   M_PI+atan( y/x); // 3rd quadrant
  else             return -1;
}

//-----------------------------------------------------------------------------
//! \brief Operator to rotate vector R2 [x,y] counter-clockwise by <phi> radians
//-----------------------------------------------------------------------------
void vectR2::operator&=(const double phi)
{
  vectR2  p(getx(),gety());
  p = p & phi;
  put(p.getx(),p.gety());
}

//=====================================
// seqR2
//=====================================
//-----------------------------------------------------------------------------
//! \brief Constructor initializing seqR1 to 0s
//-----------------------------------------------------------------------------
seqR2::seqR2(const long M)
{
  N  = M;
  xy = new vectR2[N];
  clear();
}

//-----------------------------------------------------------------------------
// constructor initializing seqR1 to <u>
//-----------------------------------------------------------------------------
seqR2::seqR2(long M, double u)
{
  long n;
  N=M;
  xy = new vectR2[N];
  for(n=0; n<N; n++)
  {
    xy[n].put(u,u);
  }
}

//-----------------------------------------------------------------------------
//! \brief Fill the seqR1 with a value 0
//-----------------------------------------------------------------------------
void seqR2::clear(void)
{
  for( long n=0; n<N; n++ ) xy[n].clear();
}

//-----------------------------------------------------------------------------
//! \brief Fill the seqR1 with a value <u>
//-----------------------------------------------------------------------------
void seqR2::fill(double u)
{
  long n;
  for(n=0; n<N; n++) xy[n].put(u,u);
}

//-----------------------------------------------------------------------------
//! \brief Fill the seqR1 with (x_0, x_1, x_2, ...)
//!        where x_n = x_{n-1} + (dx,dy)
//-----------------------------------------------------------------------------
void seqR2::inc(double x0, double y0,double dx, double dy)
{
  for( long n=0; n<N; n++ )
  {
    xy[n].put(x0,y0);
    x0 += dx;
    y0 += dy;
  }
}

//-----------------------------------------------------------------------------
//! \brief Put an order pair (u,v) into the seqR1 x at location n
//-----------------------------------------------------------------------------
int seqR2::put(long n, double u, double v)
{
  int retval=0;
  if(n>=N){
    fprintf(stderr,"n=%ld larger than seqR1 size N=%ld\n",n,N);
    retval=-1;
    }
  else xy[n].put(u,v);
  return retval;
}

//-----------------------------------------------------------------------------
//! \brief Put an R2 vector xyz into the seqR1 x at location n
//-----------------------------------------------------------------------------
int seqR2::put(long n, vectR2 xya)
{
  int retval=0;
  if(n>=N){
    fprintf(stderr,"n=%ld larger than seqR1 size N=%ld\n",n,N);
    retval=-1;
    }
  else xy[n]=xya;
  return retval;
}

//-----------------------------------------------------------------------------
//! Get a single value from the seqR1 x at location n
//-----------------------------------------------------------------------------
vectR2 seqR2::get(const long n) const
{
  vectR2 xya(0,0);
  if(n<N)xya=xy[n];
  else   fprintf(stderr,"n=%ld larger than seqR1 size N=%ld\n",n,N);
  return xya;
}

//-----------------------------------------------------------------------------
//! \brief Get a single value from the seqR1 x element x at location n
//-----------------------------------------------------------------------------
double seqR2::getx(const long n) const
{
  double u=0;
  if(n<N)u=xy[n].getx();
  else   fprintf(stderr,"n=%ld larger than x seqR1 size N=%ld\n",n,N);
  return u;
}

//-----------------------------------------------------------------------------
//! \brief Get a single value from the seqR1 x element y at location n
//-----------------------------------------------------------------------------
double seqR2::gety(const long n) const
{
  double u=0;
  if(n<N)u=xy[n].gety();
  else   fprintf(stderr,"n=%ld larger than y seqR1 size N=%ld\n",n,N);
  return u;
}

//-----------------------------------------------------------------------------
//! \brief List contents of sequence
//-----------------------------------------------------------------------------
void seqR2::list(const long start, const long end, const char *str1, const char *str2, FILE *ptr) const
{
  long n,m;
  vectR2 x;
  if(strlen(str1)>0)fprintf(ptr,"%s",str1);
  for(n=start,m=1; n<=end; n++,m++){
    x=get(n);
    fprintf(ptr,"  ");
    x.list(ptr);
    if(m%3==0)fprintf(ptr,"\n");
    }
  if(strlen(str2)>0)fprintf(ptr,"%s",str2);
}

//-----------------------------------------------------------------------------
//! \brief List contents of seqC1 using 1 digit per element
//-----------------------------------------------------------------------------
void seqR2::list1(void) const
{
  long n,m;
  for(n=0,m=1; n<N; n++,m++){
    printf("(%2.0lf,%2.0lf)   ",getx(n),gety(n));
    if(m%5==0)printf("\n");
    }
}

//-----------------------------------------------------------------------------
//! \brief List contents of seqC1 using 1 digit per element
//-----------------------------------------------------------------------------
void seqR2::list1(long start, long end) const
{
  long n,m;
  for(n=start,m=1; n<=end; n++,m++){
    printf("(%2.0lf,%2.0lf)   ",getx(n),gety(n));
    if(m%50==0)printf("\n");
    else if(m%10==0)printf(" ");
    }
}

//-----------------------------------------------------------------------------
//! \brief Return the largest pair of values in the seqR1 as measured by norm()
//-----------------------------------------------------------------------------
double seqR2::norm(const long n) const
{
  vectR2 xya=get(n);
  return xya.norm();
}

//-----------------------------------------------------------------------------
//! \brief Return the largest pair of values in the seqR1 as measured by norm()
//-----------------------------------------------------------------------------
vectR2 seqR2::max(const int verbose) const
{
  long n;
  double maxnorm=0;
  long   maxn=0;
  vectR2  maxpair;
  for(n=0; n<N; n++)if(norm(n)>maxnorm){maxnorm=norm(n); maxn=n;}
  maxpair=get(maxn);
  if(verbose)
  {
    for(n=0; n<N; n++)
    {
      if(norm(n)>=(maxnorm*0.999))
        printf("max=(%lf,%lf) at n=%ld\n",maxpair.getx(),maxpair.gety(),n);
    }
  }
  return maxpair;
}

//=====================================
// external operations
//=====================================
//-----------------------------------------------------------------------------
//! \brief Operator: return p+q
//-----------------------------------------------------------------------------
vectR2 operator+(const vectR2 p, const vectR2 q)
{
  const Eigen::Map< const Eigen::Vector2d > a( p.getdata() );
  const Eigen::Map< const Eigen::Vector2d > b( q.getdata() );
  const Eigen::Vector2d c = a + b;
  const vectR2 r( c(0), c(1) );
  return r;
}

//-----------------------------------------------------------------------------
//! \brief Operator: return p-q
//-----------------------------------------------------------------------------
vectR2 operator-(const vectR2 p, const vectR2 q)
{
  const Eigen::Map< const Eigen::Vector2d > a( p.getdata() );
  const Eigen::Map< const Eigen::Vector2d > b( q.getdata() );
  const Eigen::Vector2d c = a - b;
  const vectR2 r( c(0), c(1) );
  return r;
}

//-----------------------------------------------------------------------------
//! \brief Operator: return -p
//-----------------------------------------------------------------------------
vectR2 operator-(const vectR2 p)
{
  const Eigen::Map< const Eigen::Vector2d > a( p.getdata() );
  const Eigen::Vector2d b = -a;
  const vectR2 c( b(0), b(1) );
  return c;
}

//-----------------------------------------------------------------------------
//!  \brief operator: return <p> rotated counter-clockwise by <phi> radians
//!  \details
//!    https://stackoverflow.com/questions/17036818/
//!    https://eigen.tuxfamily.org/dox/group__TopicStorageOrders.html
//!    https://stackoverflow.com/questions/28722899/
//-----------------------------------------------------------------------------
vectR2 operator&(const vectR2 p, const double phi)
{
  double cosphi, sinphi;
  if(phi==0){ cosphi = 1;        sinphi = 0;        }
  else      { cosphi = cos(phi); sinphi = sin(phi); }
  const std::vector<double> rr = { cosphi, -sinphi, 
                                   sinphi,  cosphi  };
  const Eigen::Matrix<double, 2, 2, Eigen::RowMajor> R(rr.data()); 
//printf("     _                    _\n");
//printf("R = | %9.6lf %9.6lf  |\n", R(0,0), R(0,1) );
//printf("    |_%9.6lf %9.6lf _|\n", R(1,0), R(1,1) );
//const Eigen::Map< const Eigen::Vector2d > x( p.pairxy.data() );
  const Eigen::Map< const Eigen::Vector2d > x( p.getdata() );
  const Eigen::Vector2d y = R * x;
  const vectR2 q( y(0), y(1) );
  return q;
}

//-----------------------------------------------------------------------------
//! \brief Return the angle theta in radians between the two vectors induced by
//!        the points <p> and <q> in the plane R^2
//! \returns On SUCCESS return theta in the closed interval [0:M_PI];
//!          On ERROR   exit with value EXIT_FAILURE
//-----------------------------------------------------------------------------
double pqtheta(const vectR2 p, const vectR2 q)
{
  const double rp=p.mag(), rq=q.mag();
  double y;
  if(rp==0) return -1;
  if(rq==0) return -2;
  y = (p^q)/(rp*rq);
  if( y > +1.0 )
  {
    if(y-1 < 1e-13 ) y = 1.0;
    else{
      fprintf(stderr,"\nERROR using pqtheta(vectR2 p, vectR2 q): (p^q)/(rp*rq)=%.12lf>+1\n",y);
      exit(EXIT_FAILURE);
      }
  }
  if(y<(-1.0))
  {
    if(y+1> -1e-13 ) y = -1.0;
    else{
      fprintf(stderr,"\nERROR using pqtheta(vectR2 p, vectR2 q): (p^q)/(rp*rq)=%.12lf<-1\n",y);
      exit(EXIT_FAILURE);
      }
  }
  const double theta = acos(y);
  return theta;
}

