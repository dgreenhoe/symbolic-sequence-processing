/*============================================================================
 * Daniel J. Greenhoe
 * routines for Real Die rdieseqs
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
#include "realdie.h"

/*-------------------------------------------------------------------------
 * display real die metric table
 *-------------------------------------------------------------------------*/
int rdieseq::metrictbl(void){
  char a,b;
  for(a='A';a<='F';a++){
    for(b='A';b<='F';b++)printf("d(%c,%c)=%.1lf ",a,b,rdie_metric(a,b));
    printf("\n");
    }
  return 1;
  }

/*-------------------------------------------------------------------------
 * autocorrelation Rxx of a real die seqR1 x with 2N offset
 *-------------------------------------------------------------------------*/
int rdieseq::Rxxo(const seqR1 *rxx, const int showcount){
  const long N=getN();
  int rval;
  rval=Rxx(rxx,showcount);
  rxx->add(2*N);
  return rval;
  }

/*-------------------------------------------------------------------------
 * autocorrelation Rxx of a real die seqR1 x
 *-------------------------------------------------------------------------*/
int rdieseq::Rxx(const seqR1 *rxx, const int showcount){
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
double rdieseq::Rxx(const long m){
  const long mm=labs(m);
  const long N=getN();
  long n,nmm;
  double d,sum;
  char a,b;
  for(n=0,sum=0;n<(N+mm);n++){
    nmm=n-mm;
    a=(n  <0 || n  >=N)? 0.0 : get(n);
    b=(nmm<0 || nmm>=N)? 0.0 : get(nmm);
    d=(a==0 || b==0)?    1.0 : rdie_metric(a,b);
    sum+=d;
    }
  return -sum;
  }

/*=====================================
 * operators
 *=====================================*/
/*-------------------------------------------------------------------------
 * operator rdieseq x = dieseq y
 *-------------------------------------------------------------------------*/
void rdieseq::operator=(dieseq y){
  long n;
  const long N=getN();
  const long M=y.getN();
  char symbol;
  if(N!=M){
    fprintf(stderr,"ERROR using rdieseq x = rdieseq y operation: size of x (%ld) does not equal size of y (%ld)\n",N,M);
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
 * map die face values to R^1
 * A-->1 B-->2 C-->3 D-->4 E-->5 F-->6
 *-------------------------------------------------------------------------*/
int rdie_dietoR1(char c){
  int n,rval;
  char domain[6]={'A','B','C','D','E','F'};
  char element;
  for(n=0,rval=-1;n<6;n++)if(c==domain[n])rval=n;
  if(rval==-1){
    fprintf(stderr,"ERROR using rdie_dietoR1(char c): c=%c(0x%x) is not in the valid domain {0,A,B,C,D,E,F}\n",c,c);
    exit(EXIT_FAILURE);
    }
  return rval;
  }

/*-------------------------------------------------------------------------
 * map die face values to R^3
 *     |+1|     | 0|     | 0|     | 0|     | 0|     |-1|     | 0| 
 * A-->| 0| B-->|+1| C-->| 0| D-->| 0| E-->|-1| F-->| 0| 0-->| 0| 
 *     | 0|     | 0|     |+1|     |-1|     | 0|     | 0|     | 0| 
 *-------------------------------------------------------------------------*/
vectR3 rdie_dietoR3(char c){
  vectR3 xyz;
  switch(c){
    case 'A': xyz.put(+1, 0, 0);  break;
    case 'B': xyz.put( 0,+1, 0);  break;
    case 'C': xyz.put( 0, 0,+1);  break;
    case 'D': xyz.put( 0, 0,-1);  break;
    case 'E': xyz.put( 0,-1, 0);  break;
    case 'F': xyz.put(-1, 0, 0);  break;
    default:  
      fprintf(stderr,"ERROR using rdie_dietoR3(char c): c=%c(0x%x) is not in the valid domain {A,B,C,D,E,F}\n",c,c);
      exit(EXIT_FAILURE);
    }
  return xyz;
  }

/*-------------------------------------------------------------------------
 * map R^3 values to die face values using Lagrange Arc metric
 *-------------------------------------------------------------------------*/
rdieseq rdie_R3todie_larc(seqR3 xyz){
  long n;
  int m;
  long N=xyz.getN();
  double d[7];
  double smallestd;
  char closestface;
  vectR3 p,q[7];
  rdieseq rdie(N);

  //q[0].put(0,0,0);
  q[1]=rdie_dietoR3('A');
  q[2]=rdie_dietoR3('B');
  q[3]=rdie_dietoR3('C');
  q[4]=rdie_dietoR3('D');
  q[5]=rdie_dietoR3('E');
  q[6]=rdie_dietoR3('F');

  for(n=0; n<N; n++){
    p.put(xyz.getx(n),xyz.gety(n),xyz.getz(n));
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
    rdie.put(n,closestface);
    }
  return rdie;
  }

/*-------------------------------------------------------------------------
 * map R^3 values to die face values and (0,0,0) using Lagrange Arc metric
 *-------------------------------------------------------------------------*/
rdieseq rdie_R3todie0_larc(seqR3 xyz){
  long n;
  int m;
  long N=xyz.getN();
  double d[7];
  double smallestd;
  char closestface;
  vectR3 p,q[7];
  rdieseq rdie(N);

  q[0].put(0,0,0);
  q[1]=rdie_dietoR3('A');
  q[2]=rdie_dietoR3('B');
  q[3]=rdie_dietoR3('C');
  q[4]=rdie_dietoR3('D');
  q[5]=rdie_dietoR3('E');
  q[6]=rdie_dietoR3('F');

  for(n=0; n<N; n++){
    p.put(xyz.getx(n),xyz.gety(n),xyz.getz(n));
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
    rdie.put(n,closestface);
    }
  return rdie;
  }

/*-------------------------------------------------------------------------
 * map R^3 values to die face values using Euclidean metric
 *   0    A    B    C    D    E   F   A+..+F
 *-------------------------------------------------------------------------*/
rdieseq rdie_R3todie_euclid(seqR3 xyz){
  long n;
  int m;
  long N=xyz.getN();
  double d[7];
  double smallestd;
  char closestface;
  vectR3 p,q[7];
  rdieseq rdie(N);

  q[1]=rdie_dietoR3('A');
  q[2]=rdie_dietoR3('B');
  q[3]=rdie_dietoR3('C');
  q[4]=rdie_dietoR3('D');
  q[5]=rdie_dietoR3('E');
  q[6]=rdie_dietoR3('F');

  for(n=0; n<N; n++){
    p.put(xyz.getx(n),xyz.gety(n),xyz.getz(n));
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
    rdie.put(n,closestface);
    }
  return rdie;
  }

/*-------------------------------------------------------------------------
 * map R^3 values to die face and (0,0,0) values using Euclidean metric
 *   0    A    B    C    D    E   F   A+..+F
 *-------------------------------------------------------------------------*/
rdieseq rdie_R3todie0_euclid(seqR3 xyz){
  long n;
  int m;
  long N=xyz.getN();
  double d[7];
  double smallestd;
  char closestface;
  vectR3 p,q[7];
  rdieseq rdie(N);

  q[0]=rdie_dietoR3('0');
  q[1]=rdie_dietoR3('A');
  q[2]=rdie_dietoR3('B');
  q[3]=rdie_dietoR3('C');
  q[4]=rdie_dietoR3('D');
  q[5]=rdie_dietoR3('E');
  q[6]=rdie_dietoR3('F');

  for(n=0; n<N; n++){
    p.put(xyz.getx(n),xyz.gety(n),xyz.getz(n));
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
    rdie.put(n,closestface);
    }
  return rdie;
  }

/*-------------------------------------------------------------------------
 * map R^1 values to die face values using Euclidean metric
 *-------------------------------------------------------------------------*/
rdieseq rdie_R1todie_euclid(seqR1 xyz){
  long n;
  long N=xyz.getN();
  char closestface;
  double p;
  rdieseq rdie(N);

  for(n=0; n<N; n++){
    p = xyz.get(n);
    if(p<1.5)       closestface='A';
    else if(p>=5.5) closestface='F';
    else            closestface=(char)(p+0.5-1)+'A';
    rdie.put(n,closestface);
    }
  return rdie;
  }

/*-------------------------------------------------------------------------
 * real die metric d(a,b)
 *   d(a,b) | 0    A    B    C    D    E     F   (b)
 *   -------|-------------------------------------
 *      a= 0| 0    1    1    1    1    1    1   
 *      a= A| 1    0    1    1    1    1    2
 *      a= B| 1    1    0    1    1    2    1
 *      a= C| 1    1    1    0    2    1    1
 *      a= D| 1    1    1    2    0    1    1
 *      a= E| 1    1    2    1    1    0    1
 *      a= F| 1    2    1    1    1    1    0
 * On success return d(a,b). On error return -1.
 *-------------------------------------------------------------------------*/
double rdie_metric(char a, char b){
  int ra=rdie_dietoR1(a);
  int rb=rdie_dietoR1(b);
  double d;

  if(ra<0)fprintf(stderr,"a=%c(0x%x) not in domain of rdie metric d(a,b)\n",a,a);
  if(rb<0)fprintf(stderr,"b=%c(0x%x) not in domain of rdie metric d(a,b)\n",b,b);

       if(ra<0)     d=-1;
  else if(rb<0)     d=-1;
  else if(ra==rb)   {d=0;  }
  else if(ra==0)    d=1;
  else if(rb==0)    d=1;
  else if(ra+rb==7) d=2;
  else              d=1;
  return d;
  }

/*-------------------------------------------------------------------------
 * real die metric p(x,y) where x and y are rdie sequences computed as
 * p(x,y) = d(x0,y0) + d(x1,y1) + d(x2,y2) + ... + d(x{N-1},y{N-1})
 * where d(a,b) is defined above.
 * On success return d(x,y). On error return -1.
 *-------------------------------------------------------------------------*/
double rdie_metric(rdieseq x, rdieseq y){
  double rval,d;
  long n;
  long N=x.getN();
  long M=y.getN();
  long NM=(N<M)?M:N; //NM = the larger of N and M
  for(n=0,d=0;n<NM;n++){
    rval=rdie_metric(x.get(n),y.get(n));
    if(rval<0){d+=0.0; printf("rval=%lf ",rval);}
    else d+=rval;
    }
  if(N!=M){
    fprintf(stderr,"ERROR using rdie_metric(rdieseq x,rdieseq y): size of x (%ld) does not equal the size of y (%ld).\n",N,M);
    exit(EXIT_FAILURE);
    }
  return d;
  }

