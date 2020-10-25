/*============================================================================
 * Daniel J. Greenhoe
 * routines for Real spin spinseqs
 *============================================================================*/
/*=====================================
 * headers
 *=====================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "main.h"
#include "symseq.h"
#include "r1.h"
#include "r2.h"
#include "r3.h"
#include "r4.h"
#include "r6.h"
#include "c1.h"
#include "euclid.h"
#include "larc.h"
#include "die.h"
#include "spinner.h"

/*-------------------------------------------------------------------------
 * display spinner metric table
 *-------------------------------------------------------------------------*/
int spinseq::metrictbl(void){
  char a,b;
  for(a='A';a<='F';a++){
    for(b='A';b<='F';b++)printf("d(%c,%c)=%.1lf ",a,b,spin_metric(a,b));
    printf("\n");
    }
  return 1;
  }

/*-------------------------------------------------------------------------
 * autocorrelation Rxx of a spinner seqR1 x with 2N offset
 *-------------------------------------------------------------------------*/
int spinseq::Rxxo(seqR1 *rxx, const int showcount){
  long N = getN();
  int rval;
  rval=Rxx(rxx,showcount);
  rxx->add(2*N);
  return rval;
  }

/*-------------------------------------------------------------------------
 * autocorrelation Rxx of a spinner seqR1 x
 *-------------------------------------------------------------------------*/
int spinseq::Rxx(seqR1 *rxx, const int showcount){
  long m;
  const long N=getN();
  int rval=0;
  double rxxm;
  if(showcount)fprintf(stderr,"  Calculate %ld auto-correlation values ... n=",2*N+1);
  for(m=-N;m<=N;m++){
    if(showcount)fprintf(stderr,"%8ld",m+N);
    rxxm=Rxx(m);
    if(rxxm>0)rval=-1;
    rxx->put(m+N,rxxm);
    if(showcount)fprintf(stderr,"\b\b\b\b\b\b\b\b");
    }
  if(showcount)fprintf(stderr,"%8ld .... done.\n",m+N);
  return rval;
  }

/*-------------------------------------------------------------------------
 * autocorrelation Rxx(m)
 *-------------------------------------------------------------------------*/
double spinseq::Rxx(const long m){
  const long mm=labs(m);
  const long N=getN();
  long n,nmm;
  double d,sum;
  char a,b;
  for(n=0,sum=0;n<(N+mm);n++){
    nmm=n-mm;
    a=(n  <0 || n  >=N)? 0.0 : get(n);
    b=(nmm<0 || nmm>=N)? 0.0 : get(nmm);
    d=(a==0 || b==0)?    1.0 : spin_metric(a,b);
    sum+=d;
    }
  return -sum;
  }

/*-------------------------------------------------------------------------
 * downsample sequence by a factor of <factor>
 *-------------------------------------------------------------------------*/
spinseq spinseq::downsample(int factor){
  const long N=getN();
  long n,m;
  long M;
  if(factor<1){
    fprintf(stderr,"ERROR using dieseq::downsample: factor=%d must be at least 1\n",factor);
    exit(EXIT_FAILURE);
    }
  M=N/factor;
  spinseq newseq(M);
  for(n=0,m=0; m<M; n+=factor,m++)newseq.put(m,get(n));
  return newseq;
  }


/*-------------------------------------------------------------------------
 * map spin face values to R^1
 * A-->1  B-->2  C-->3  D-->4  E-->5  F-->6
 *-------------------------------------------------------------------------*/
seqR1 spinseq::spintoR1(void){
  const long N=getN();
  long n;
  seqR1 y(N);
  for(n=0; n<N; n++)y.put(n,spin_spintoR1(get(n)));
  return y;
  }

/*-------------------------------------------------------------------------
 * map spin face values to R^2 sequence
 *-------------------------------------------------------------------------*/
seqR2 spinseq::spintoR2(void){
  const long N=getN();
  long n;
  seqR2 seqR2(N);
  for(n=0; n<N; n++)seqR2.put(n,spin_spintoR2(get(n)));
  return seqR2;
  }


/*=====================================
 * operators
 *=====================================*/
/*-------------------------------------------------------------------------
 * operator spinseq x = dieseq y
 *-------------------------------------------------------------------------*/
void spinseq::operator=(spinseq y){
  long n;
  const long N=getN();
  const long M=y.getN();
  char symbol;
  if(N!=M){
    fprintf(stderr,"\nERROR using spinseq x = spinseq y: size of x (%ld) is smaller than size of y (%ld)\n",N,M);
    exit(EXIT_FAILURE);
    }
  for(n=0;n<N;n++){
    symbol=y.get(n);
    put(n,symbol);
    }
  }

/*=====================================
 * external operations
 *=====================================*/
/*-------------------------------------------------------------------------
 * map spin face values to R^1
 *-------------------------------------------------------------------------*/
double spin_spintoR1(char c){
  double rval;
  switch(c){
    case 'A': rval = 1.0;  break;
    case 'B': rval = 2.0;  break;
    case 'C': rval = 3.0;  break;
    case 'D': rval = 4.0;  break;
    case 'E': rval = 5.0;  break;
    case 'F': rval = 6.0;  break;
    default:  
      fprintf(stderr,"ERROR using spin_spintoR1(c): c=%c(0x%x) is not in the valid domain {A,B,C,D,E,F}\n",c,c);
      exit(EXIT_FAILURE);
    }
  return rval;
  }

/*-------------------------------------------------------------------------
 * map spinner face values to R^2
 *-------------------------------------------------------------------------*/
vectR2 spin_spintoR2(char c){
  vectR2 xy;
  switch(c){
    case 'A': xy.put(0,         -1.0);  break;
    case 'B': xy.put(+sqrt(3)/2,-0.5);  break;
    case 'C': xy.put(+sqrt(3)/2,+0.5);  break;
    case 'D': xy.put( 0,        +1.0);  break;
    case 'E': xy.put(-sqrt(3)/2,+0.5);  break;
    case 'F': xy.put(-sqrt(3)/2,-0.5);  break;
    default:  
      fprintf(stderr,"ERROR: c=%c(0x%x) is not in the valid domain {0,A,B,C,D,E,F} in spin_spintoR2(char c)\n",c,c);
      exit(EXIT_FAILURE);
    }
  return xy;
  }

/*-------------------------------------------------------------------------
 * map R^2 values to spin face values using Lagrange Arc distance
 *-------------------------------------------------------------------------*/
spinseq spin_R2tospin_larc(seqR2 xy){
  long n;
  int m;
  long N=xy.getN();
  double d[7];
  double smallestd;
  char closestface;
  vectR2 p,q[7];
  spinseq rspin(N);

  //q[0].put(0,0,0);
  q[1]=spin_spintoR2('A');
  q[2]=spin_spintoR2('B');
  q[3]=spin_spintoR2('C');
  q[4]=spin_spintoR2('D');
  q[5]=spin_spintoR2('E');
  q[6]=spin_spintoR2('F');

  for(n=0; n<N; n++){
    p.put(xy.getx(n),xy.gety(n));
    smallestd=larc_metric(p,q[1]);
    closestface='A';
    for(m=2;m<7;m++){
      d[m] = larc_metric(p,q[m]);
      if(((m&0x01) && (d[m]<smallestd)) || ((!(m&0x01)) && (d[m]<=smallestd))){
        //----------------------------     ----------------------------------
        // bias odd samples                  bias even samples
        // towards smaller values            towards larger values
        smallestd=d[m];
        closestface='A'+m-1;
        }
      }
    rspin.put(n,closestface);
    }
  return rspin;
  }

/*-------------------------------------------------------------------------
 * map R^3 values to spin face values and (0,0,0) using Lagrange Arc metric
 *-------------------------------------------------------------------------*/
spinseq spin_R2tospin0_larc(seqR2 xy){
  long n;
  int m;
  long N=xy.getN();
  double d[7];
  double smallestd;
  char closestface;
  vectR2 p,q[7];
  spinseq rspin(N);

  q[0].put(0,0);
  q[1]=spin_spintoR2('A');
  q[2]=spin_spintoR2('B');
  q[3]=spin_spintoR2('C');
  q[4]=spin_spintoR2('D');
  q[5]=spin_spintoR2('E');
  q[6]=spin_spintoR2('F');

  for(n=0; n<N; n++){
    p.put(xy.getx(n),xy.gety(n));
    smallestd=larc_metric(p,q[0]);
    //smallestd=ae_metric(1,p,q[0]);
    closestface='0';
    for(m=1;m<7;m++){
      d[m] = larc_metric(p,q[m]);
      if(((m&0x01) && (d[m]<smallestd)) || ((!(m&0x01)) && (d[m]<=smallestd))){
        //----------------------------     ----------------------------------
        // bias odd samples                  bias even samples
        // towards smaller values            towards larger values
        smallestd=d[m];
        closestface='A'+m-1;
        }
      }
    rspin.put(n,closestface);
    }
  return rspin;
  }

/*-------------------------------------------------------------------------
 * map R^2 values to spin face values using Euclidean metric
 *   0    A    B    C    D    E   F   A+..+F
 *-------------------------------------------------------------------------*/
spinseq spin_R2tospin_euclid(seqR2 xy){
  long n;
  int m;
  long N=xy.getN();
  double d[7];
  double smallestd;
  char closestface;
  vectR2 p,q[7];
  spinseq rspin(N);

  q[1]=spin_spintoR2('A');
  q[2]=spin_spintoR2('B');
  q[3]=spin_spintoR2('C');
  q[4]=spin_spintoR2('D');
  q[5]=spin_spintoR2('E');
  q[6]=spin_spintoR2('F');

  for(n=0; n<N; n++){
    p.put(xy.getx(n),xy.gety(n));
    smallestd=ae_metric(1,p,q[1]);
    closestface='A';
    for(m=2;m<=6;m++){
      d[m] = ae_metric(1,p,q[m]);
      if(((m&0x01) && (d[m]<smallestd)) || ((!(m&0x01)) && (d[m]<=smallestd))){
        //----------------------------     ----------------------------------
        // bias towards smaller values         bias towards larger values
        smallestd=d[m];
        closestface='A'+m-1;
        }
      }
    rspin.put(n,closestface);
    }
  return rspin;
  }

/*-------------------------------------------------------------------------
 * map R^2 values to spin face and (0,0) values using Euclidean metric
 *   0    A    B    C    D    E   F   A+..+F
 *-------------------------------------------------------------------------*/
spinseq spin_R2tospin0_euclid(seqR3 xy){
  long n;
  int m;
  long N=xy.getN();
  double d[7];
  double smallestd;
  char closestface;
  vectR2 p,q[7];
  spinseq rspin(N);

  q[0]=spin_spintoR2('0');
  q[1]=spin_spintoR2('A');
  q[2]=spin_spintoR2('B');
  q[3]=spin_spintoR2('C');
  q[4]=spin_spintoR2('D');
  q[5]=spin_spintoR2('E');
  q[6]=spin_spintoR2('F');

  for(n=0; n<N; n++){
    p.put(xy.getx(n),xy.gety(n));
    smallestd=ae_metric(1,p,q[0]);
    closestface='0';
    for(m=1;m<=6;m++){
      d[m] = ae_metric(1,p,q[m]);
      //if(d[m]<smallestd)
      //if(d[m]<=smallestd)
      //if(((m&0x01) && (d[m]<smallestd)) || ((!(m&0x01)) && (d[m]<smallestd)))
      //if(((m&0x01) && (d[m]<=smallestd)) || ((!(m&0x01)) && (d[m]<=smallestd)))
      if(((m&0x01) && (d[m]<smallestd)) || ((!(m&0x01)) && (d[m]<=smallestd))){
        //----------------------------     ----------------------------------
        // bias towards smaller values         bias towards larger values
        smallestd=d[m];
        closestface='A'+m-1;
        }
      //if(m&0x01){ (alternative coding)
      //  if(d[m]<smallestd){
      //    smallestd=d[m];
      //    closestface='A'+m-1;
      //    }}
      //else
      //  if(d[m]<=smallestd){
      //    smallestd=d[m];
      //    closestface='A'+m-1;
      //    }
      }
    rspin.put(n,closestface);
    }
  return rspin;
  }

/*-------------------------------------------------------------------------
 * map R^1 values to spin face values using Euclidean metric
 *-------------------------------------------------------------------------*/
spinseq spin_R1tospin_euclid(seqR1 xy){
  long n;
  long N=xy.getN();
  char closestface;
  double p;
  spinseq rspin(N);

  for(n=0; n<N; n++){
    p = xy.get(n);
    if(p<1.5)       closestface='A';
    else if(p>=5.5) closestface='F';
    else            closestface=(char)(p+0.5-1)+'A';
    rspin.put(n,closestface);
    }
  return rspin;
  }

/*-------------------------------------------------------------------------
 * spinner metric d(a,b)
 *   d(a,b) | A    B    C    D    E     F   (b)
 *   -------|--------------------------------
 *      a= A| 0    1    2    3    2    1
 *      a= B| 1    0    1    2    3    2
 *      a= C| 2    1    0    1    2    3
 *      a= D| 3    2    1    0    1    2
 *      a= E| 2    3    2    1    0    1
 *      a= F| 1    2    3    2    1    0
 * On success return d(a,b). On error return -1.
 *-------------------------------------------------------------------------*/
double spin_metric(char a, char b){
  double ra=spin_spintoR1(a);
  double rb=spin_spintoR1(b);
  double d;

  d=(double)fabs(ra-rb);
  if     (d>3.5) d=2.0;
  else if(d>4.5) d=1.0;
  return d;
  }
