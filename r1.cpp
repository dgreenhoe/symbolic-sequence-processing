//=============================================================================
// Daniel J. Greenhoe
// routines for seqR1 routines
//=============================================================================
//=====================================
// headers
//=====================================
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "r1.h"

//-----------------------------------------------------------------------------
//! \brief Constructor initializing sequence to 0
//-----------------------------------------------------------------------------
seqR1::seqR1(long M){
  N = M;
  x = new double[N];
  clear();
}

//-----------------------------------------------------------------------------
//! \brief Constructor initializing sequence to u
//-----------------------------------------------------------------------------
seqR1::seqR1(long M, double u){
  N = M;
  x = new double[N];
  fill( u );
}

//-----------------------------------------------------------------------------
//! \brief Constructor initializing sequence to an array list
//! \param[in] list: pointer to array of initialization values
//-----------------------------------------------------------------------------
seqR1::seqR1(long M, double *list)
{
  N = M;
  x = new double[N];
  for(long n=0; n<N; n++) x[n] = list[n];
}

//-----------------------------------------------------------------------------
//! \brief fill the sequence with a value 0
//-----------------------------------------------------------------------------
void seqR1::clear(void)
{
  fill( 0.0 );
}

//-----------------------------------------------------------------------------
//! \brief fill the sequence with a value u
//! \param[in] u : Value to fill sequence with
//-----------------------------------------------------------------------------
void seqR1::fill(double u)
{
  long n;
  for(n=0; n<N; n++)x[n]=u;
}

//-----------------------------------------------------------------------------
//! \brief fill the sequence with Hanning filter coefficients
//-----------------------------------------------------------------------------
void seqR1::hanning(void)
{
  long n;
  for(n=0; n<N; n++) x[n]=0.5*(1-cos(2*M_PI*n/(N-1)));
}

//-----------------------------------------------------------------------------
//! \brief Convert low pass filter to high pass filter
//-----------------------------------------------------------------------------
void seqR1::lptohp(void)
{
  long n;
  double s=1;
  for(n=0; n<N; n++,s*=-1) x[n]*=s;
}

//-----------------------------------------------------------------------------
//! \brief put a single value u into the seqR1 x at location n
//! \param[in] n : Index value
//! \param[in] u : Value to be put in sequence at index n
//-----------------------------------------------------------------------------
const int seqR1::put(const long n, const double u)
{
  if(n<N){
    x[n]=u;
    return 0;
  }
  else return -1;
}

//-----------------------------------------------------------------------------
//! \brief Put the value u into the sequence x from location start to end
//! \param[in] u     : Value to be put in sequence
//! \param[in] start : Starting index value
//! \param[in] end   : Ending index value
//-----------------------------------------------------------------------------
int seqR1::put(const long start, const long end, const double u)
{
  long n;
  int retval=0;
  if(start<0||end>=N||start>end){
    fprintf(stderr,"ERROR using put(%ld,%ld,%lf)\n",start,end,u);
    retval=-1;
  }
  else for(n=start;n<=end;n++) x[n]=u;
  return retval;
}

//-----------------------------------------------------------------------------
//! \brief add a single value u into the seqR1 x at location n
//! \param[in] n: Index of sequence
//! \param[in] u: Add value u to element n of sequence
//-----------------------------------------------------------------------------
int seqR1::add(const long n, const double u){
  if(n<N){
    x[n]+=u;
    return 0;
  }
  else return -1;
}

//-----------------------------------------------------------------------------
//! \brief add a single value u to every element in the sequence x
//! \param[in] u: Add value u to each element of sequence
//-----------------------------------------------------------------------------
const int seqR1::add(const double u){
  long n;
  for(n=0; n<N; n++) x[n] += u;
  return 0;
}

//-----------------------------------------------------------------------------
//! \brief get a single value from the seqR1 x at location n
//! \param[in] n: sequence index value
//-----------------------------------------------------------------------------
double seqR1::get(long n) const
{
  if(n>=N){
    fprintf(stderr,"\nERROR using seqR1::get(n): index n=%ld larger than sequence length N=%ld\n",n,N);
    exit(EXIT_FAILURE);
  }
  return x[n];
}

//-----------------------------------------------------------------------------
//! \brief shift seqR1 n elements to the right inserting zeros on the left
//! \details example: if x = [ a b c d e f ] (N=6), then shiftR(2) results in
//!             x = [ 0 0 a b c d ] (N=6).
//! \param[in] n: number of elements to shift by
//-----------------------------------------------------------------------------
void seqR1::shiftR(long n)
{
  long m;
  for(m=N-1;m-n>=0;m--)x[m]=x[m-n];
  for(m=0;m<n;m++)     x[m]=0;
}

//-----------------------------------------------------------------------------
//! \brief shift seqR1 n elements to the left inserting zeros on the right
//! \details example: if x = [ a b c d e f ] (N=6), then shiftL(2) results in
//!                      x = [ c d e f 0 0 ] (N=6).
//-----------------------------------------------------------------------------
void seqR1::shiftL(long n){
  long m;
  for(m=0;m<N-n;m++) x[m]=x[m+n];
  for(m=N-2;m<N;m++) x[m]=0;
}

//-----------------------------------------------------------------------------
//! \brief return the maximum value of the seqR1 x
//! \param[in] start starting index of sequence
//! \param[in] end   ending index of sequence
//-----------------------------------------------------------------------------
double seqR1::max(const long start,const long end) const
{
  long n=start;
  double u,max;
  max=x[n];
  n++;
  for(;n<=end;n++){
    u=x[n];
    if(u>max) max=u;
  }
  return max;
}

//-----------------------------------------------------------------------------
//! \brief return the minimum value of the seqR1 x
//! \param[in] start starting index of sequence
//! \param[in] end   ending index of sequence
//-----------------------------------------------------------------------------
double seqR1::min(const long start,const long end) const
{
  long n=start;
  double u,min;
  min=x[start];
  n++;
  for(;n<=end;n++){
    u=x[n];
    if(u<min) min=u;
  }
  return min;
}

//-----------------------------------------------------------------------------
//! \brief return the first argument n yielding the maximum value x_n of the seqR1 x
//-----------------------------------------------------------------------------
long seqR1::argmax(void)
{
  long n,arg;
  double u,max;
  max=x[0];
  for(n=1; n<N; n++){
    u=x[n];
    if(u>max){
      max=u;
      arg=n;
    }
  }
  return arg;
}

//-----------------------------------------------------------------------------
//! \brief return the first argument n yielding the minimum value x_n of the seqR1 x
//-----------------------------------------------------------------------------
long seqR1::argmin(void)
{
  long n,arg;
  double u,min;
  min=x[0];
  for(n=1; n<N; n++){
    u=x[n];
    if(u<min){ 
      min=u;
      arg=n;
    }
  }
  return arg;
}

//-----------------------------------------------------------------------------
//! \brief return the lp norm of the seqR1 x
//! \details References: Giles (1987), page 9, 1.18 Example.
//!   John Robilliard Giles (1987)
//!   Introduction to the Analysis of Metric Spaces. 
//!   Number 3 in Australian Mathematical Society lecture series. 
//!   Cambridge University Press, Cambridge, ISBN 0521359287. 
//!   https://books.google.com/books?vid=ISBN0521359287&pg=PA9
//! \param[in] r: norm parameter
//-----------------------------------------------------------------------------
double seqR1::norm( double r )
{
  long n;
  double u,y=0;
  for(n=0;n<N;n++){
    u=fabs(x[n]);
    y+=pow(u,r);
  }
  return pow(y,1.0/r);
}

//-----------------------------------------------------------------------------
//! \brief return the Euclidean norm of the seqR1 x
//!        equivalent to norm(2), see above, but possibly more efficient
//! \reference: Giles (1987), page 3, 1.6 Examples.
//-----------------------------------------------------------------------------
double seqR1::norm(void)
{
  long n;
  double u,y=0;
  for(n=0;n<N;n++){
    u=x[n];
    y+=u*u;
  }
  return sqrt(y);
}

//-----------------------------------------------------------------------------
//! \brief return the supremum norm of the seqR1 x
//!        equivalent to norm(+infty), see above
//!        reference: Giles (1987), page 5, 1.9 Examples.
//-----------------------------------------------------------------------------
double seqR1::supnorm(void)
{
  double u,max;
  max=fabs(x[0]);
  for(long n=1; n<N; n++){
    u=fabs(x[n]);
    if(u>max) max=u;
  }
  return max;
}

//-----------------------------------------------------------------------------
//! \brief return the infinimum norm of the seqR1 x 
//!        equivalent to norm(-infty), see above
//-----------------------------------------------------------------------------
double seqR1::infnorm(void)
{
  double u,min;
  min=fabs(x[0]);
  for(long n=1; n<N; n++){
    u=fabs(x[n]);
    if(u<min) min=u;
  }
  return min;
}

//-----------------------------------------------------------------------------
//! \brief fill the seqR1 with pseudo-random numbers 
//!        uniformly distributed in the closed interval [0:1]
//!        using seed value <seed>
//! \param[in] seed: pseudo-random generator seed value
//-----------------------------------------------------------------------------
void seqR1::randomize( unsigned seed )
{
  double r;
  srand(seed);
  for(long n=0; n<N; n++){
    r=(double)rand();
    x[n]=r/((double)RAND_MAX+1.0);
  }
}

//-----------------------------------------------------------------------------
//! \brief fill the seqR1 with pseudo-random integers
//!        uniformly distributed in the closed interval [min:max]
//!        using seed value <seed>
//! \param[in] seed: pseudo-random generator seed value
//! \param[in] min:  minimum value
//! \param[in] max:  maximum value
//-----------------------------------------------------------------------------
void seqR1::randomize(unsigned seed, int min , int max )
{
  int r;
  unsigned u;
  int  order=abs(max-min)+1;
  srand(seed);
  for(long n=0; n<N; n++){
    r=rand();
    u = r%order+min;
    x[n]=(double)u;
  }
}

//-----------------------------------------------------------------------------
//! \brief round each element in seqR1 to nearest integer
//-----------------------------------------------------------------------------
void seqR1::increment(long n)
{
  x[n] += 1;
}

//-----------------------------------------------------------------------------
//! \brief Round each element in seqR1 to nearest integer
//-----------------------------------------------------------------------------
void seqR1::round(void)
{
  for(long n=0; n<N; n++) x[n] = (double)((long)(x[n]+0.5));
}

//-----------------------------------------------------------------------------
//! \brief Print contents of sequence to stream <*ptr> 
//!        first printing <*str1>, 
//!        followed by the sequence [x_n|n=start...end]
//!        followed by the <*str2>
//! \param[in]  start: sequence starting index value
//! \param[in]  end  : sequence ending index value
//! \param[in]  str1 : string to print before the sequence
//! \param[in]  str2 : string to print after the sequence
//! \param[out] ptr  : output string
//-----------------------------------------------------------------------------
void seqR1::list( const long  start, const long  end, const char *str1, const char *str2, FILE *ptr )
{
  long n,m;
  if(ptr!=NULL){
    if(strlen(str1)>0)fprintf(ptr,"%s",str1);
    n=start;
    fprintf(ptr,"%8ld:  %10.6lf",n,get(n));
    n++;
    for(m=1; n<=end; n++,m++){
        if(n%5==0)fprintf(ptr,"\n%8ld: ",n);
        fprintf(ptr," %10.6lf",get(n));
      }
    if(strlen(str2)>0)fprintf(ptr,"%s",str2);
  }
}

//-----------------------------------------------------------------------------
//! \brief List contents of sequence in a format convenient for LaTeX typesetting
//! \param[in] start: starting index of sequence
//! \param[in] end:   ending index of sequence
//-----------------------------------------------------------------------------
void seqR1::listx(long start, long end)
{
  long n,m;
  for(n=start,m=0; n<=end; n++,m++){
    printf("(%6ld, %10.5lf)%%%6ld\n",n,x[n],m);
  }
}


//-----------------------------------------------------------------------------
//! \brief List contents of seqR1 in long format
//! \param[in] start: starting index of sequence
//! \param[in] end:   ending index of sequence
//-----------------------------------------------------------------------------
void seqR1::listL(void)
{
  long n,m;
  for(n=0,m=1; n<N; n++,m++){
    printf("%12.9lf ",x[n]);
    if(m%5==0)printf("\n");
  }
}

//-----------------------------------------------------------------------------
//! \brief List contents of seqR1 using 1 digit per element
//! \param[in] start: starting index of sequence
//! \param[in] end:   ending index of sequence
//-----------------------------------------------------------------------------
void seqR1::list1(long start, long end)
{
  long n,m;
  for(n=start,m=1; n<=end; n++,m++){
    printf("%1.0lf",get(n));
    if(m%50==0)printf("\n");
    else if(m%10==0)printf(" ");
  }
}

//-----------------------------------------------------------------------------
//! \brief List contents of seqR1 using integer format
//! \param[in] start: starting index of sequence
//! \param[in] end:   ending index of sequence
//-----------------------------------------------------------------------------
void seqR1::listi(long start, long end)
{
  long n,m;
  for(n=start,m=1; n<=end; n++,m++){
    printf("%.0lf ",get(n));
    if(m%50==0)printf("\n");
    else if(m%10==0)printf(" ");
  }
}

//-----------------------------------------------------------------------------
//! \brief Normalize sequence
//-----------------------------------------------------------------------------
double seqR1::normalize(void)
{
  long n;
  double sum,scale;
  for(n=0,sum=0;n<N;n++) sum += fabs(x[n]);
  scale = 1.0/sum;
  for(n=0;n<N;n++) x[n] *= scale;
  return scale;  
}


//-----------------------------------------------------------------------------
//! \brief Return the largest pair of values in the seqR1 as measured by norm()
//! \returns Largest value 
//-----------------------------------------------------------------------------
double seqR1::max(const int mode) const
{
  long n;
  double maxnorm=0;
  long   maxn=0;
  double maxval,u;
  for(n=0; n<N; n++){
    u=fabs(get(n));
    if(u>maxnorm){maxnorm=u; maxn=n;}
  }
  maxval=get(maxn);
  if(mode){
    for(n=0; n<N; n++){
      u=fabs(get(n));
      if(u>=(maxnorm*0.999))printf("max=%lf at n=%ld\n",maxval,n);
    }
  }
  return maxval;
}

//-----------------------------------------------------------------------------
//! \brief List all values (n,x_n) such that x_n > <threshold>
//! \param[in]  threshold : threshold above which a value will be listed
//! \param[in]  start     : starting index of sequence
//! \param[in]  end       : ending index of sequence
//! \param[in]  str1      : String to print before listing
//! \param[in]  str2      : String to print after listing
//! \param[out] ptr       : Pointer to output file
//-----------------------------------------------------------------------------
long seqR1::gt(const double threshold,const long start, const long end, const char *str1, const char *str2, FILE *ptr){
  long n;
  double v;
  long count=0;
  char buf[256];
  
  if(strlen(str1)>0)if(ptr!=NULL){fprintf(ptr,"%s",str1);}
  for(n=start; n<=end; n++){
    v=get(n);
    if(v>threshold){
      count++;
      sprintf(buf," (%5ld,%9.6lf)",n,v);
      if(ptr!=NULL){fprintf(ptr,"%s",buf);}
      if(count%4==0&&count!=0)if(ptr!=NULL){fprintf(ptr,"\n");}
    }
  }
  if(count%4!=0)if(ptr!=NULL){fprintf(ptr,"\n");}
  sprintf(buf,"total number of values > %lf is %ld\n",threshold,count);
  if(ptr!=NULL)if(ptr!=NULL){fprintf(ptr,"%s",buf);}
  if(strlen(str2)>0)if(ptr!=NULL){fprintf(ptr,"%s",str2);}
  return count;
}

//-----------------------------------------------------------------------------
//! \brief List all values (n,x_n) such that x_n >= <threshold>
//! \param[in]  threshold : threshold above which a value will be listed
//! \param[in]  start     : starting index of sequence
//! \param[in]  end       : ending index of sequence
//! \param[in]  str1      : String to print before listing
//! \param[in]  str2      : String to print after listing
//! \param[out] ptr       : Pointer to output file
//-----------------------------------------------------------------------------
long seqR1::gte(const double threshold,const long start, const long end, const char *str1, const char *str2, FILE *ptr){
  long n;
  double v;
  long count=0;
  char buf[256];
  if(strlen(str1)>0){
    if(ptr!=NULL)fprintf(ptr,"%s",str1);
  }
  for(n=start; n<=end; n++){
    v=get(n);
    if(v>=threshold){
      count++;
      sprintf(buf," (%5ld,%9.6lf)",n,v);
      if(ptr!=NULL)fprintf(ptr,"%s",buf);
      if(count%4==0&&count!=0){if(ptr!=NULL)fprintf(ptr,"\n");}
    }
  }
  if(count%4!=0)if(ptr!=NULL)fprintf(ptr,"\n");
  sprintf(buf,"total number of values >= %lf is %ld\n",threshold,count);
  if(ptr!=NULL)fprintf(ptr,"%s",buf);
  if(strlen(str2)>0){
    if(ptr!=NULL)fprintf(ptr,"%s",str2);
  }
  return count;
}

//-----------------------------------------------------------------------------
// \brief List contents of seqR1 using 1 digit per element
//-----------------------------------------------------------------------------
//void seqR1::test(void)
//{
//  printf("seqR1::test N=%ld x[0,1,2] = %lf %lf %lf ...\n",N,x[0],x[1],x[2]);
//}

//-----------------------------------------------------------------------------
//! \brief operator seqR1 x = seqR1 y
//! \param[in] y: Set sequence x to sequence y
//-----------------------------------------------------------------------------
void seqR1::operator=(seqR1 y)
{
  const long M=y.getN();
  long n;
  if(N!=M){
    fprintf(stderr,"ERROR using seqR1 x = seqR1 y operation: lengths of x (%ld) and y (%ld) differ\n",N,M);
    exit(EXIT_FAILURE);
  }
  for(n=0;n<N;n++)x[n]=y.get(n);
}
 
//-----------------------------------------------------------------------------
//! \brief element in seqR1 x = p
//! \param[in] p: Set each element in sequence x to the value p
//-----------------------------------------------------------------------------
void seqR1::operator=(const double p)
{
  long n;
  for(n=0;n<N;n++)x[n]=p;
}
 
//-----------------------------------------------------------------------------
//! \brief operator: return x += y where y is another seqR1 of length M
//-----------------------------------------------------------------------------
void seqR1::operator+=(seqR1 y)
{
  long n,M=y.N,NM;
  if(N<=M)NM=N;
  else    NM=M;
  for(n=0;n<NM;n++)x[n]+=y.x[n];
}

//-----------------------------------------------------------------------------
//! \brief operator: return x -= y where y is another seqR1 of length M
//-----------------------------------------------------------------------------
void seqR1::operator-=(seqR1 y)
{
  long n,NM;
  const long M=y.getN();
  if(N<=M)NM=N;
  else    NM=M;
  for(n=0;n<NM;n++)x[n]-=y.x[n];
}

//-----------------------------------------------------------------------------
//! \brief operator: return x = alpha*x where alpha is a constant
//! \param[in] alpha: Multiply each element in sequence x by the value alpha
//-----------------------------------------------------------------------------
void seqR1::operator*=(double alpha)
{
  long n;
  for(n=0;n<N;n++)x[n]*=alpha;
}

