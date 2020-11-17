//============================================================================
// Daniel J. Greenhoe
// normed linear space R^2
//============================================================================
//=====================================
// headers
//=====================================
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <Eigen/Dense>
#include "r1.h"
#include "r3.h"


//=====================================
// VectR3
//=====================================*/
//-----------------------------------------------------------------------------
//! \brief Calculate magnitude of vector
//-----------------------------------------------------------------------------
double vectR3::mag(void) const 
{
  const double x = getx();
  const double y = gety();
  const double z = getz();
  const double magSq  = x*x + y*y + z*z;
  const double magVal = sqrt( magSq );
  return magVal;
}

//-----------------------------------------------------------------------------
//! \brief Convert polar (r, theta, phi) coordinate to rectangular (x,y,z) coordinate
//-----------------------------------------------------------------------------
void vectR3::polartoxyz(const double r, const double theta, const double phi)
{
  const double x = r*cos(phi)*cos(theta);
  const double y = r*cos(phi)*sin(theta);
  const double z = r*sin(phi);
  otriple::put( x, y, z );
}

//-----------------------------------------------------------------------------
//! \brief Add vector q to vectR3 vector
//-----------------------------------------------------------------------------
void vectR3::operator+=(const vectR3 q)
{
  Eigen::Map< Eigen::Vector3d > a(   getdataa() );
  const Eigen::Map< const Eigen::Vector3d > b( q.getdata() );
  a += b;
}

//-----------------------------------------------------------------------------
//! \brief Subtract vector q from vectR3 vector
//-----------------------------------------------------------------------------
void vectR3::operator-=(const vectR3 q)
{
  Eigen::Map< Eigen::Vector3d > a(   getdataa() );
  const Eigen::Map< const Eigen::Vector3d > b( q.getdata() );
  a -= b;
}

//=====================================
// seqR3
//=====================================*/
//-----------------------------------------------------------------------------
//! \brief Constructor initializing seqR3 to 0
//-----------------------------------------------------------------------------
seqR3::seqR3(long M)
{
  N=M;
  seqr3 = new vectR3[N];
  clear();
}

//-----------------------------------------------------------------------------
//! \brief Constructor initializing seqR3 to <u>
//-----------------------------------------------------------------------------
seqR3::seqR3(long M, double u)
{
  N = M;
  seqr3 = new vectR3[N];
  fill( u );
}

//-----------------------------------------------------------------------------
//! \brief Fill the seqR3 with a value 0
//-----------------------------------------------------------------------------
void seqR3::clear(void)
{
  for( long n=0; n<N; n++) put( n, 0 );
}

//-----------------------------------------------------------------------------
//! \brief Fill the seqR3 with a value <u>
//-----------------------------------------------------------------------------
void seqR3::fill(double u)
{
  long n;
  for(n=0; n<N; n++) put( n, u );
}

//-----------------------------------------------------------------------------
//! \brief Put a single value <u> into the seqR3 x at location n
//-----------------------------------------------------------------------------
int seqR3::put(long n, double u, double v, double w)
{
  if(n<N){
    seqr3[n].otriple::put( u, v, w );
    return 0;
    }
  else{   
    fprintf(stderr,"n=%ld larger than seqR3 size N=%ld\n",n,N);
    return -1;
    }
}

//-----------------------------------------------------------------------------
//! \brief Put a single value <u> into the seqR3 x at location n
//-----------------------------------------------------------------------------
int seqR3::put(long n, vectR3 abc)
{
  if(n<N)
  {
    (seqr3[n]).otriple::put( abc.getx(), abc.gety(), abc.gety() );
    return 0;
  }
  else
  {   
    fprintf(stderr,"n=%ld larger than seqR3 size N=%ld\n",n,N);
    return -1;
  }
}

//-----------------------------------------------------------------------------
//! \brief Get a single value from the sequence x at location n
//-----------------------------------------------------------------------------
vectR3 seqR3::get(long n)
{
  return seqr3[n];
}

//-----------------------------------------------------------------------------
//! \brief Get the x element from the sequence at location n
//-----------------------------------------------------------------------------
double seqR3::getx(long n)
{
  double u=0;
  if(n<N)u = seqr3[n].getx();
  else   fprintf(stderr,"n=%ld larger than x seqR3 size N=%ld\n",n,N);
  return u;
}

//-----------------------------------------------------------------------------
//! \brief Get the y element from the sequence at location n
//-----------------------------------------------------------------------------
double seqR3::gety(long n){
  double u=0;
  if(n<N)u = seqr3[n].gety();
  else   fprintf(stderr,"n=%ld larger than y seqR3 size N=%ld\n",n,N);
  return u;
  }

//-----------------------------------------------------------------------------
//! \brief Get the z element from the sequence at location n
//-----------------------------------------------------------------------------
double seqR3::getz(long n){
  double u=0;
  if(n<N)u = seqr3[n].getz();
  else   fprintf(stderr,"n=%ld larger than z seqR3 size N=%ld\n",n,N);
  return u;
  }

//-----------------------------------------------------------------------------
//! \brief list contents of sequence
//-----------------------------------------------------------------------------
void seqR3::list(const long start, const long end, const char *str1, const char *str2, FILE *ptr){
  long n,m;
  if(ptr!=NULL){
    if(strlen(str1)>0) fprintf(ptr,"%s",str1);
    for(n=start,m=1; n<=end; n++,m++){
      fprintf(ptr,"(%6.3lf,%6.3lf,%6.3lf) ", seqr3[n].getx(), seqr3[n].gety(), seqr3[n].getz() );
      if(m%3==0)fprintf(ptr,"\n");
      }
    if(strlen(str2)>0)fprintf(ptr,"%s",str2);
    }
  }

//-----------------------------------------------------------------------------
//! \brief list contents of seqR3 using 1 digit per element
//-----------------------------------------------------------------------------
void seqR3::list1(void){
  long n,m;
  for(n=0,m=1; n<N; n++,m++){
    printf("(%2.0lf,%2.0lf,%2.0lf)   ", seqr3[n].getx(), seqr3[n].gety(), seqr3[n].getz() );
    if(m%5==0)printf("\n");
    }
  }
void seqR3::list1(long start, long end){
  long n,m;
  for(n=start,m=1; n<=end; n++,m++){
    printf("(%2.0lf,%2.0lf,%2.0lf)   ", seqr3[n].getx(), seqr3[n].gety(), seqr3[n].getz() );
    if(m%50==0)printf("\n");
    else if(m%10==0)printf(" ");
    }
  }


//=====================================
// external operations
//=====================================
//-----------------------------------------------------------------------------
//! \brief operator: return p+q
//-----------------------------------------------------------------------------
vectR3 operator+(vectR3 p, vectR3 q)
{
  const Eigen::Map< const Eigen::Vector3d > a( p.getdataa() );
  const Eigen::Map< const Eigen::Vector3d > b( q.getdataa() );
  const Eigen::Vector3d c = a + b;
  const vectR3 r( c(0), c(1), c(2) );
  return r;
}

//-----------------------------------------------------------------------------
//! \brief operator: return p-q
//-----------------------------------------------------------------------------
vectR3 operator-(vectR3 p, vectR3 q){
  const Eigen::Map< const Eigen::Vector3d > a( p.getdataa() );
  const Eigen::Map< const Eigen::Vector3d > b( q.getdataa() );
  const Eigen::Vector3d c = a - b;
  const vectR3 r( c(0), c(1), c(2) );
  return r;
  }

//-----------------------------------------------------------------------------
//! \brief  operator: return -p
//-----------------------------------------------------------------------------
vectR3 operator-(vectR3 p)
{
  const Eigen::Map< const Eigen::Vector3d > a( p.getdataa() );
  const Eigen::Vector3d b = -a;
  const vectR3 q( b(0), b(1), b(2) );
  return q;
}

//-----------------------------------------------------------------------------
//! \brief return the angle theta in radians between the two vectors induced by 
//! the points <p> and <q> in the space R^3.
//! \returns On SUCCESS return theta in the closed interval [0:PI];
//!          On ERROR   return negative value or exit with value EXIT_FAILURE
//-----------------------------------------------------------------------------
double pqtheta(const vectR3 p, const vectR3 q)
{
  const double rp=p.mag();
  const double rq=q.mag();
  if(rp==0) return -1;
  if(rq==0) return -2;
  const double y = (p^q)/(rp*rq);
  if(y>+1)  {fprintf(stderr,"\nERROR using pqtheta(vectR3 p, vectR3 q): (p^q)/(rp*rq)=%lf>+1\n",y); exit(EXIT_FAILURE);}
  if(y<-1)  {fprintf(stderr,"\nERROR using pqtheta(vectR3 p, vectR3 q): (p^q)/(rp*rq)=%lf<-1\n",y); exit(EXIT_FAILURE);}
  const double theta = acos(y);
  return theta;
}

//-----------------------------------------------------------------------------
//! \brief Return the minimum element of the tupple
//-----------------------------------------------------------------------------
double otriple::min(void) const
{
  const Eigen::Map< const Eigen::Vector3d > abc( getdata() );
  int row, col;
  const double min = abc.minCoeff( &row, &col );
  return min;
}

//-----------------------------------------------------------------------------
//! \brief Return the maximum element of the tupple
//-----------------------------------------------------------------------------
double otriple::max(void) const
{
  const Eigen::Map< const Eigen::Vector3d > abc( getdata() );
  int row, col;
  const double max = abc.maxCoeff( &row, &col );
  return max;
}
