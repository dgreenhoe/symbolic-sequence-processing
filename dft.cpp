/*============================================================================
 * Daniel J. Greenhoe
 * DSP operations on R^1
 *============================================================================*/
/*=====================================
 * headers
 *=====================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "main.h"
#include "r1.h"
#include "r2.h"
#include "r4.h"
#include "r6.h"
#include "c1.h"
#include "c4.h"
#include "c6.h"
#include "dft.h"

/*-------------------------------------------------------------------------
 * <n>th element of unitary Discrete Fourier Transform (DFT) of seqR1 <x>
 *-------------------------------------------------------------------------*/
complex dftn(seqR1 *x, const long n){
  const long N=x->getN();
  const double scale=1.0/sqrt((double)N);
  long m;
  double theta,xm;
  complex yn(0,0);
  //printf("\n");
  if((n>=N)||(n<0))fprintf(stderr,"WARNING for y(n)=DFT(x,n): n (%ld) is outside domain of x [0,%ld]. Returning (0,0).\n",n,N-1);
  else{
yn.put(0,0);
    for(m=0;m<N;m++){
      theta = (-2.0*PI/(double)N)*(double)n*(double)m;
      xm=x->get(m);
      yn+=(xm*expi(theta));
      //printf("n=%2ld, m=%2ld, xm=%2.0lf theta=%lf (%3.0lf), expi(theta)=(%+6.3lf,%+6.3lf) yn=%5.3lf\n",n,m,xm,theta,theta/PI*180,yn.real(),yn.imag(), yn.mag());
      }
    yn *= scale;
    }
  return yn;
  }

/*-------------------------------------------------------------------------
 * <n>th element of unitary Discrete Fourier Transform (DFT) of seqC1 <x>
 *-------------------------------------------------------------------------*/
complex dftn(seqC1 *x, const long n){
  const long N=x->getN();
  const double scale=1.0/sqrt((double)N);
  long m;
  double theta;
  complex xm;
  complex yn(0,0);

  if((n>=N)||(n<0))fprintf(stderr,"WARNING for y(n)=DFT(x,n): n (%ld) is outside domain of x [0,%ld]. Returning (0,0).\n",n,N-1);
  else{
    for(m=0;m<N;m++){
      theta = (-2.0*PI/(double)N)*(double)n*(double)m;
      xm=x->get(m);
      yn+=(xm*expi(theta));
      }
    yn *= scale;
    }
  return yn;
  }

/*-------------------------------------------------------------------------
 * <n>th element of unitary Discrete Fourier Transform (DFT) of seqR4 <x>
 *-------------------------------------------------------------------------*/
vectC4 dftn(seqR4 *x, const long n){
  const long N=x->getN();
  const double scale=1.0/sqrt((double)N);
  long m;
  double theta;
  complex etheta;
  vectR4 xm;
  vectC4 yn(0,0);
  vectC4 ethetaxm;

  if((n>=N)||(n<0))fprintf(stderr,"WARNING for y(n)=DFT(x,n): n (%ld) is outside domain of x [0,%ld]. Returning (0,0).\n",n,N-1);
  else{
    for(m=0;m<N;m++){
      theta = (-2.0*PI/(double)N)*(double)n*(double)m;
      etheta= expi(theta);
      xm=x->get(m);
      ethetaxm = etheta*xm;
      yn+=ethetaxm;
      }
    yn *= scale;
    }
  return yn;
  }

/*-------------------------------------------------------------------------
 * <n>th element of unitary Discrete Fourier Transform (DFT) of seqR6 <x>
 *-------------------------------------------------------------------------*/
vectC6 dftn(seqR6 *x, const long n){
  const long N=x->getN();
  const double scale=1.0/sqrt((double)N);
  long m;
  double theta;
  complex etheta;
  vectR6 xm;
  vectC6 yn(0,0);
  vectC6 ethetaxm;

  if((n>=N)||(n<0))fprintf(stderr,"WARNING for y(n)=DFT(x,n): n (%ld) is outside domain of x [0,%ld]. Returning (0,0).\n",n,N-1);
  else{
    for(m=0;m<N;m++){
      theta = (-2.0*PI/(double)N)*(double)n*(double)m;
      etheta= expi(theta);
      xm=x->get(m);
      ethetaxm = etheta*xm;
      yn+=ethetaxm;
      }
    yn *= scale;
    }
  return yn;
  }

/*-------------------------------------------------------------------------
 * unitary Discrete Fourier Transform (DFT) from seqR1 to seqC1
 *                       N-1
 * note that X_{N-n} = a sum x_m exp(-i2pi(N-n)m/N)
 *                       m=0
 *                       N-1
 *                   = a sum x_m exp(-i2pi(-n)m/N) exp(-i2pi(N)m/N)
 *                       m=0
 *                     (  N-1                     )*
 *                   = (a sum x_m exp(-i2pi nm/N) )   for ((x_m)) real
 *                     (  m=0                     )
 *                   = (X_n)*                         for ((x_m)) real
 *-------------------------------------------------------------------------*/
int dft(seqR1 *x, seqC1 *y){
  const long Nx=x->getN();
  const long Ny=y->getN();
  const long N=(Nx<Ny)?Nx:Ny;
  int retval=0;
  long n=0;
  complex yn,ync;

  y->clear();
  if(Nx!=Ny){
    fprintf(stderr,"WARNING for y=DFT(x): lengths of x (%ld) and y (%ld) differ; will use length=%ld",Nx,Ny,N);
    retval=-1;
    }
  printf(" %10ld %10ld",n,N-n);
  for(n=0;n<=(N/2+1);n++){
    printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%10ld %10ld",n,N-n);
    yn =dftn(x,n);
    ync=yn.conj();
    y->put(n,yn);
    if(n!=0)y->put(N-n,ync);
    }
  return retval;
  }

/*-------------------------------------------------------------------------
 * unitary Discrete Fourier Transform (DFT) from seqC1 to seqC1
 *-------------------------------------------------------------------------*/
int dft(seqC1 *x, seqC1 *y){
  const long Nx=x->getN();
  const long Ny=y->getN();
  const long N=(Nx<Ny)?Nx:Ny;
  int retval=0;
  long n=0;
  complex yn;

  y->clear();
  if(Nx!=Ny){
    fprintf(stderr,"WARNING for y=DFT(x): lengths of x (%ld) and y (%ld) differ; will use length=%ld",Nx,Ny,N);
    retval=-1;
    }
  printf(" %10ld",n);
  for(n=0;n<N;n++){
    printf("\b\b\b\b\b\b\b\b\b\b%10ld",n);
    yn=dftn(x,n);
    y->put(n,yn);
    }
  return retval;
  }

/*-------------------------------------------------------------------------
 * unitary Discrete Fourier Transform (DFT) from seqR1 to seqC1
 *-------------------------------------------------------------------------*/
int dft(seqR4 *x, seqC4 *y)
{
  const long Nx=x->getN();
  const long Ny=y->getN();
  const long N=(Nx<Ny)?Nx:Ny;
  long n=0;
  vectC4 yn,ync;
  int retval=0;

  if(Nx!=Ny){
    fprintf(stderr,"WARNING for y=DFT(x): lengths of x (%ld) and y (%ld) differ; will use length=%ld",Nx,Ny,N);
    retval=-1;
    }
  y->clear();
  printf(" %10ld %10ld",n,N-n);
  for(n=0;n<=(N/2+1);n++){
    printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%10ld %10ld",n,N-n);
    yn=dftn(x,n);
    ync = yn.conj();
    y->put(n,yn);
    if(n!=0)y->put(N-n,ync);
    }
  return retval;
}

/*-------------------------------------------------------------------------
 * unitary Discrete Fourier Transform (DFT) from seqR1 to seqC1
 *-------------------------------------------------------------------------*/
int dft(seqR6 *x, seqC6 *y){
  const long Nx=x->getN();
  const long Ny=y->getN();
  const long N=(Nx<Ny)?Nx:Ny;
  long n=0;
  vectC6 yn;
  int retval=0;

  if(Nx!=Ny){
    fprintf(stderr,"WARNING for y=DFT(x): lengths of x (%ld) and y (%ld) differ; will use length=%ld",Nx,Ny,N);
    retval=-1;
    }
  y->clear();
  printf(" %10ld %10ld",n,N-n);
  for(n=0;n<=(N/2+1);n++){
    printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%10ld %10ld",n,N-n);
    yn=dftn(x,n);
    y->put(n,yn);
    if(n!=0)y->put(N-n,yn.conj());
    }
  return retval;
  }

