//=============================================================================
//! \brief Daniel J. Greenhoe
//! \brief normed linear space R^2
//=============================================================================
//=====================================
//! \brief headers
//=====================================
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include <Eigen/Dense>  // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=89325
#include "r1.h"
#include "r6.h"
typedef Eigen::Matrix< double, 6, 1 > Vector6d;
//=====================================
//! \brief osix
//=====================================
//-----------------------------------------------------------------------------
//! \brief osix constructors
//-----------------------------------------------------------------------------
osix::osix(void)
{
  x.at(0) = 0.0;
  x.at(1) = 0.0;
  x.at(2) = 0.0;
  x.at(3) = 0.0;
  x.at(4) = 0.0;
  x.at(5) = 0.0;
}

osix::osix(const double u0, const double u1, const double u2, const double u3, const double u4, const double u5)
{
  x.at(0) = u0;
  x.at(1) = u1;
  x.at(2) = u2;
  x.at(3) = u3;
  x.at(4) = u4;
  x.at(5) = u5;
}
   
osix::osix(const double u)
{
  x.at(0) = u;
  x.at(1) = u;
  x.at(2) = u;
  x.at(3) = u;
  x.at(4) = u;
  x.at(5) = u;
}

//-----------------------------------------------------------------------------
//! \brief osix put member functions
//-----------------------------------------------------------------------------
void osix::put(const double u)
{
  x.at(0) = u;
  x.at(1) = u;
  x.at(2) = u;
  x.at(3) = u;
  x.at(4) = u;
  x.at(5) = u;
}

void osix::put(const osix u)
{
  x.at(0) = u.get(0);
  x.at(1) = u.get(1);
  x.at(2) = u.get(2);
  x.at(3) = u.get(3);
  x.at(4) = u.get(4);
  x.at(5) = u.get(5);
} 

//-----------------------------------------------------------------------------
//! \brief return the 6-tuple value
//-----------------------------------------------------------------------------
osix osix::get(void) const
{
  osix u;
  u.put( get(0) );
  u.put( get(1) );
  u.put( get(2) );
  u.put( get(3) );
  u.put( get(4) );
  u.put( get(5) );
  return u;
}

//-----------------------------------------------------------------------------
//! \brief return the minimum element of the 6 tupple
//-----------------------------------------------------------------------------
double osix::min(void) const
{
  const Eigen::Map< const Vector6d > a( getdata() );
  const double minVal = a.minCoeff();
  return minVal;
}

//-----------------------------------------------------------------------------
//! \brief return the maximum element of the 6 tupple
//-----------------------------------------------------------------------------
double osix::max(void) const
{
  const Eigen::Map< const Vector6d > a( getdata() );
  const double maxVal = a.maxCoeff();
  return maxVal;
}

//-----------------------------------------------------------------------------
//! \brief print the tuple
//-----------------------------------------------------------------------------
void osix::list(const char *str1, const char *str2)
{
  int i;
  if(strlen(str1)!=0) printf("%s",str1);
  putchar('(');
  for( i=0; i<5; i++ ) printf("%9.6lf,", get(i));
  printf("%9.6lf)", get(5) );
  if(strlen(str2)!=0) printf("%s",str2);
}

//=====================================
//! \brief vectR6 functions
//=====================================
//-----------------------------------------------------------------------------
//! \brief return the 6-tuple value
//-----------------------------------------------------------------------------
vectR6 vectR6::get(void) const 
{
  vectR6 u;
  int i;
  for( i=0; i<6; i++ ) u.put( i, get(i) );
  return u;
}

//-----------------------------------------------------------------------------
//! \brief magnitude
//-----------------------------------------------------------------------------
const double vectR6::mag(void) const
{
  const Eigen::Map< const Vector6d > a( getdata() );
  return a.norm();
}

//-----------------------------------------------------------------------------
//! \brief Multiply the vector by a scalar a
//-----------------------------------------------------------------------------
vectR6 vectR6::mpy(const double a)
{
  vectR6 w;
  const Eigen::Map< const Vector6d > vv( getdata() );
  Eigen::Map< Vector6d > ww( w.getdataa() );
  ww = a * vv;
  return w;
}

//-----------------------------------------------------------------------------
//! \brief operator: +=
//-----------------------------------------------------------------------------
void vectR6::operator+=(const vectR6 q)
{
  Eigen::Map< Vector6d > pp( getdataa() );
  const Eigen::Map< const Vector6d > qq( q.getdata() );
  pp = pp + qq;
}

//-----------------------------------------------------------------------------
//! \brief operator: -=
//-----------------------------------------------------------------------------
void vectR6::operator-=(const vectR6 q)
{
  Eigen::Map< Vector6d > pp( getdataa() );
  const Eigen::Map< const Vector6d > qq( q.getdata() );
  pp = pp - qq;
}

//-----------------------------------------------------------------------------
//! \brief operator: -=
//-----------------------------------------------------------------------------
void vectR6::operator*=(const double a)
{
  vectR6  p=get();
  p = a*p;
  put(p);
}


//-----------------------------------------------------------------------------
//! \brief operator: a*x
//-----------------------------------------------------------------------------
vectR6 operator*(const double a, const vectR6 x)
{
  vectR6 y;
  const Eigen::Map< const Vector6d > xx( x.getdata() );
  Eigen::Map< Vector6d > yy( y.getdataa() );
  yy = a * xx;
  return y;
}

//-----------------------------------------------------------------------------
//! \brief operator: dot product of p and q
//-----------------------------------------------------------------------------
double operator^(const vectR6 p, const vectR6 q)
{
  const  Eigen::Map< const Vector6d > pp( p.getdata() );
  const  Eigen::Map< const Vector6d > qq( q.getdata() );
  double innerProduct = pp.adjoint() * qq;
  return innerProduct;
}

//=====================================
//! \brief seqR6
//=====================================
//-----------------------------------------------------------------------------
//! \brief constructor initializing seqR1 to 0
//-----------------------------------------------------------------------------
seqR6::seqR6(const long M)
{
  N = M;
  x = (vectR6 *)malloc( N * sizeof(vectR6) );
  clear();
}

//-----------------------------------------------------------------------------
//! \brief constructor initializing seqR1 to <u>
//-----------------------------------------------------------------------------
seqR6::seqR6(const long M, const double u)
{
  N=M;
  x = (vectR6 *)malloc(N*sizeof(vectR6));
  fill( u );
}

//-----------------------------------------------------------------------------
//! \brief fill the seqR1 with a value <u>
//-----------------------------------------------------------------------------
void seqR6::fill(const double u)
{
  long n;
  for(n=0; n<N; n++)x[n].put(u);
}

//-----------------------------------------------------------------------------
//! \brief fill the seqR1 with a values 0
//-----------------------------------------------------------------------------
void seqR6::clear(void)
{
  fill(0.0);
}

//-----------------------------------------------------------------------------
//! \brief put a single value u=(u1,u2,u3,u4,u5,u6) into the seqR1 x at location n
//-----------------------------------------------------------------------------
int seqR6::put(const long n, const double u1, const double u2, const double u3, const double u4, const double u5, const double u6)
{
  int retval=0;
  if(n<N)x[n].put(u1,u2,u3,u4,u5,u6);
  else{
    fprintf(stderr,"n=%ld larger than seqR1 size N=%ld\n",n,N);
    retval=-1;
    }
  return retval;
}

int seqR6::put(const long n, const vectR6 u)
{
  int retval=0;
  if(n<N)x[n].put(u);
  else{   
    fprintf(stderr,"n=%ld larger than seqR6 size N=%ld\n",n,N);
    retval=-1;
    }
  return retval;
}

//-----------------------------------------------------------------------------
//! \brief list contents of sequence
//-----------------------------------------------------------------------------
void seqR6::list(const long start, const long end, const char* str1, const char *str2, FILE *ptr)
{
  long n,m;
  vectR6 p;
  if(ptr!=NULL){
    if(strlen(str1)>0)fprintf(ptr,"%s",str1);
    for(n=start,m=1; n<=end; n++,m++){
      p=x[n];
      fprintf(ptr,"(%5.2lf,%5.2lf,%5.2lf,%5.2lf,%5.2lf,%5.2lf) ",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6());
      if(m%2==0)fprintf(ptr,"\n");
      }
    if(strlen(str2)>0)fprintf(ptr,"%s",str2);
    }
}

//-----------------------------------------------------------------------------
//! \brief list contents of seqR1 using 1 digit per element
//-----------------------------------------------------------------------------
void seqR6::list1(const long start, const long end, const char *str1, const char *str2,FILE *ptr)
{
  long n,m;
  vectR6 p;
  if(ptr!=NULL){
    if(strlen(str1)>0)fprintf(ptr,"%s",str1);
    for(n=start,m=1; n<=end; n++,m++){
      p=x[n];
      fprintf(ptr," %1.0lf%1.0lf%1.0lf%1.0lf%1.0lf%1.0lf",p.get1(),p.get2(),p.get3(),p.get4(),p.get5(),p.get6());
      if(m%10==0)fprintf(ptr,"\n");
      }
    if(strlen(str2)>0)fprintf(ptr,"%s",str2);
    }
}

//=====================================
// External operations
//=====================================
//-----------------------------------------------------------------------------
//! \brief operator: return p+q
//-----------------------------------------------------------------------------
vectR6 operator+(const vectR6 p, const vectR6 q)
{
  vectR6 y;
  Eigen::Map< Vector6d > yy( y.getdataa() );
  const Eigen::Map< const Vector6d > pp( p.getdata() );
  const Eigen::Map< const Vector6d > qq( q.getdata() );
  yy = pp + qq;
  return y;
}

//-----------------------------------------------------------------------------
//! \brief operator: return p-q
//-----------------------------------------------------------------------------
vectR6 operator-(const vectR6 p, const vectR6 q)
{
  vectR6 y;
  Eigen::Map< Vector6d > yy( y.getdataa() );
  const Eigen::Map< const Vector6d > pp( p.getdata() );
  const Eigen::Map< const Vector6d > qq( q.getdata() );
  yy = pp - qq;
  return y;
}

//-----------------------------------------------------------------------------
//! \brief operator: return -p
//-----------------------------------------------------------------------------
vectR6 operator-(const vectR6 p)
{
  vectR6 q;
  Eigen::Map< Vector6d > qq( q.getdataa() );
  const Eigen::Map< const Vector6d > pp( p.getdata() );
  qq = - pp;
  return q;
}

//-----------------------------------------------------------------------------
//! \brief return the angle theta in radians between the two vectors induced by 
//!        the points <p> and <q> in the space R^6.
//! \returns On SUCCESS return theta in the closed interval [0:PI];
//!          On ERROR   return negative value or exit with value EXIT_FAILURE
//-----------------------------------------------------------------------------
double pqtheta(const vectR6 p, const vectR6 q)
{
  const double rp = p.r();
  const double rq = q.r();
  double y,theta;
  if(rp==0) 
  {
    fprintf(stderr,"\nERROR using pqtheta(vectR6 p, vectR6 q): |p| = %lf\n", rp);
    exit(EXIT_FAILURE);
  }
  if(rq==0)
  {
    fprintf(stderr,"\nERROR using pqtheta(vectR6 p, vectR6 q): |q| = %lf\n", rq);
    exit(EXIT_FAILURE);
  }
  y = (p^q)/(rp*rq);
  if(y>+1)  
  {
    fprintf(stderr,"\nERROR using pqtheta(vectR6 p, vectR6 q): (p^q)/(rp*rq)=%lf>+1\n",y); 
    exit(EXIT_FAILURE);
  }
  if(y<-1)  
  {
    fprintf(stderr,"\nERROR using pqtheta(vectR6 p, vectR6 q): (p^q)/(rp*rq)=%lf<-1\n",y); 
    exit(EXIT_FAILURE);
  }
  theta = acos(y);
  return theta;
}
