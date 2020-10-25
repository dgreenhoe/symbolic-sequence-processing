/*============================================================================
 * Daniel J. Greenhoe
 * routines for Real Die fdieseqs
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
#include "c6.h"
#include "euclid.h"
#include "larc.h"
#include "die.h"
#include "fairdie.h"

/*-------------------------------------------------------------------------
 * display real die metric table
 *-------------------------------------------------------------------------*/
int fdieseq::metrictbl(void){
  char a,b;
  for(a='A';a<='F';a++){
    for(b='A';b<='F';b++)printf("d(%c,%c)=%.1lf ",a,b,fdie_metric(a,b));
    printf("\n");
    }
  return 1;
  }


/*-------------------------------------------------------------------------
 * autocorrelation Rxx of a fair die seqR1 x with 2N offset
 *-------------------------------------------------------------------------*/
int fdieseq::Rxxo(seqR1 *rxx, const int showcount) const
{
  const long N=getN();
  int rval;
  rval=Rxx(rxx,showcount);
  rxx->add(2*N);
  return rval;
}

/*-------------------------------------------------------------------------
 * autocorrelation Rxx of a fair die seqR1 x
 *-------------------------------------------------------------------------*/
int fdieseq::Rxx(seqR1 *rxx, const int showcount) const
{
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
double fdieseq::Rxx(const long m) const
{
  const long mm=labs(m);
  const long N=getN();
  long n,nmm;
  double d,sum;
  char a,b;
  for(n=0,sum=0;n<(N+mm);n++){
    nmm=n-mm;
    a = (n  <0 || n  >=N)? 0.0 : get(n);
    b = (nmm<0 || nmm>=N)? 0.0 : get(nmm);
    d = (a ==0 || b  ==0)? 1.0 : fdie_metric(a,b);
    sum+=d;
    }
  return -sum;
}


//

/*=====================================
 * operators
 *=====================================*/
/*-------------------------------------------------------------------------
 * operator fdieseq x = dieseq y
 *-------------------------------------------------------------------------*/
void fdieseq::operator=(dieseq y){
  long n;
  const long N=getN();
  const long M=y.getN();
  if(N!=M){
    fprintf(stderr,"\nERROR using fdieseq x = fdieseq y: size of x (%ld) is smaller than size of y (%ld)\n",N,M,N);
    exit(EXIT_FAILURE);
    }
  for(n=0;n<N;n++)put(n,y.get(n));
  }

/*-------------------------------------------------------------------------
 * map die face values to R^6 sequence
 *-------------------------------------------------------------------------*/
seqR6 fdieseq::dietoR6(void){
  const long N=getN();
  long n;
  char yc;
  seqR6 seq6(N);
  for(n=0; n<N; n++)seq6.put(n,die_dietoR6(get(n)));
  return seq6;
  }

/*=====================================
 * external operations
 *=====================================*/
/*-------------------------------------------------------------------------
 * map R^6 values to die face values using scaled Euclidean metric
 *-------------------------------------------------------------------------*/
fdieseq fdie_R6todie_ae(seqR6 seq){
  long n;
  int m;
  long N=seq.getN();
  double d[7];
  double smallestd;
  char closestface;
  vectR6 p,q[7];
  fdieseq fdie(N);

  q[1]=die_dietoR6('A');
  q[2]=die_dietoR6('B');
  q[3]=die_dietoR6('C');
  q[4]=die_dietoR6('D');
  q[5]=die_dietoR6('E');
  q[6]=die_dietoR6('F');

  for(n=0; n<N; n++){
    p.put(seq.get1(n),seq.get2(n),seq.get3(n),seq.get4(n),seq.get5(n),seq.get6(n));
    smallestd=ae_metric(sqrt(2.0)/2.0,p,q[1]);
    closestface='A';
    for(m=2;m<7;m++){
      d[m] = ae_metric(sqrt(2.0)/2.0,p,q[m]);
      if(((m&0x01) && (d[m]<smallestd)) || ((!(m&0x01)) && (d[m]<=smallestd))){
        //----------------------------     ----------------------------------
        // bias odd samples                  bias even samples
        // towards smaller values            towards larger values
        smallestd=d[m];
        closestface='A'+m-1;
        }
      }
    fdie.put(n,closestface);
    }
  return fdie;
  }

/*-------------------------------------------------------------------------
 * map R^6 values to die face values and (0,0,0) using Lagrange Arc metric
 *-------------------------------------------------------------------------*/
fdieseq fdie_R6todie0_ae(seqR6 xyz){
  long n;
  int m;
  long N=xyz.getN();
  double d[7];
  double smallestd;
  char closestface;
  vectR6 p,q[7];
  fdieseq fdie(N);

  q[0].put(0);
  q[1]=die_dietoR6('A');
  q[2]=die_dietoR6('B');
  q[3]=die_dietoR6('C');
  q[4]=die_dietoR6('D');
  q[5]=die_dietoR6('E');
  q[6]=die_dietoR6('F');

  for(n=0; n<N; n++){
    p.put(xyz.get1(n),xyz.get2(n),xyz.get3(n),xyz.get4(n),xyz.get5(n),xyz.get6(n));
    smallestd=ae_metric(sqrt(2.0)/2.0,p,q[0]);
    //smallestd=ae_metric(1,p,q[0]);
    closestface='0';
    for(m=1;m<7;m++){
      d[m] = ae_metric(sqrt(2.0)/2.0,p,q[m]);
      if(((m&0x01) && (d[m]<smallestd)) || ((!(m&0x01)) && (d[m]<=smallestd))){
        //----------------------------     ----------------------------------
        // bias odd samples                  bias even samples
        // towards smaller values            towards larger values
        smallestd=d[m];
        closestface='A'+m-1;
        }
      }
    fdie.put(n,closestface);
    }
  return fdie;
  }

/*-------------------------------------------------------------------------
 * map R^6 values to die face values using Euclidean metric
 *-------------------------------------------------------------------------*/
fdieseq fdie_R6todie_euclid(seqR6 xyz){
  const long N=xyz.getN();
  long n;
  int m;
  double d[7];
  double smallestd;
  char closestface;
  vectR6 p,q[7];
  fdieseq fdie(N);


  //q[0]=die_dietoR6('0'); use by fdie_R6todie0_euclid(seqR6 xyz)
  q[1]=die_dietoR6('A');
  q[2]=die_dietoR6('B');
  q[3]=die_dietoR6('C');
  q[4]=die_dietoR6('D');
  q[5]=die_dietoR6('E');
  q[6]=die_dietoR6('F');

  for(n=0; n<N; n++){
    p.put(xyz.get1(n),xyz.get2(n),xyz.get3(n),xyz.get4(n),xyz.get5(n),xyz.get6(n));
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
    fdie.put(n,closestface);
    }
  return fdie;
  }

/*-------------------------------------------------------------------------
 * map R^6 values to die face values using Euclidean metric
 *-------------------------------------------------------------------------*/
fdieseq fdie_R6todie_larc(seqR6 xyz){
  long n;
  int m;
  long N=xyz.getN();
  double d[7];
  double smallestd;
  char closestface;
  vectR6 p,q[7];
  fdieseq fdie(N);

  q[1]=die_dietoR6('A');
  q[2]=die_dietoR6('B');
  q[3]=die_dietoR6('C');
  q[4]=die_dietoR6('D');
  q[5]=die_dietoR6('E');
  q[6]=die_dietoR6('F');

  for(n=0; n<N; n++){
    p.put(xyz.get1(n),xyz.get2(n),xyz.get3(n),xyz.get4(n),xyz.get5(n),xyz.get6(n));
    smallestd=larc_metric(p,q[1]);
    closestface='A';
    for(m=2;m<=6;m++){
      d[m] = larc_metric(p,q[m]);
      if(((m&0x01) && (d[m]<smallestd)) || ((!(m&0x01)) && (d[m]<=smallestd))){
        //----------------------------     ----------------------------------
        // bias towards smaller values         bias towards larger values
        smallestd=d[m];
        closestface='A'+m-1;
        }
      }
    fdie.put(n,closestface);
    }
  return fdie;
  }

/*-------------------------------------------------------------------------
 * map R^6 values to die face and (0,0,0) values using Euclidean metric
 *-------------------------------------------------------------------------*/
fdieseq fdie_R6todie0_euclid(seqR6 xyz){
  long n;
  int m;
  long N=xyz.getN();
  double d[7];
  double smallestd;
  char closestface;
  vectR6 p,q[7];
  fdieseq fdie(N);

  q[0]=die_dietoR6('0');
  q[1]=die_dietoR6('A');
  q[2]=die_dietoR6('B');
  q[3]=die_dietoR6('C');
  q[4]=die_dietoR6('D');
  q[5]=die_dietoR6('E');
  q[6]=die_dietoR6('F');

  for(n=0; n<N; n++){
    p.put(xyz.get1(n),xyz.get2(n),xyz.get3(n),xyz.get4(n),xyz.get5(n),xyz.get6(n));
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
    fdie.put(n,closestface);
    }
  return fdie;
  }


/*-------------------------------------------------------------------------
 * map R^1 values to die face values using Euclidean metric
 *-------------------------------------------------------------------------*/
fdieseq fdie_R1todie_euclid(seqR1 xyz){
  long n;
  long N=xyz.getN();
  char closestface;
  double p;
  fdieseq fdie(N);

  for(n=0; n<N; n++){
    p = xyz.get(n);
    if(p<1.5)       closestface='A';
    else if(p>=5.5) closestface='F';
    else            closestface=(char)(p+0.5-1)+'A';
    fdie.put(n,closestface);
    }
  return fdie;
  }

/*-------------------------------------------------------------------------
 * real die metric d(a,b)
 *   d(a,b) | 0    A    B    C    D    E     F   (b)
 *   -------|-------------------------------------
 *      a= 0| 0    1    1    1    1    1    1
 *      a= A| 1    0    1    1    1    1    1
 *      a= B| 1    1    0    1    1    1    1
 *      a= C| 1    1    1    0    1    1    1
 *      a= D| 1    1    1    1    0    1    1
 *      a= E| 1    1    1    1    1    0    1
 *      a= F| 1    1    1    1    1    1    0
 * On success return d(a,b). On error return -1.
 *-------------------------------------------------------------------------*/
double fdie_metric(char a, char b){
  int ra=die_dietoR1(a);
  int rb=die_dietoR1(b);
  double d;

  if(ra<0)fprintf(stderr,"a=%c(0x%x) not in domain of fdie metric d(a,b)\n",a,a);
  if(rb<0)fprintf(stderr,"b=%c(0x%x) not in domain of fdie metric d(a,b)\n",b,b);

       if(ra<0)     d=-1.0;
  else if(rb<0)     d=-1.0;
  else if(ra==rb)   d= 0.0;
  else if(ra==0)    d= 1.0;
  else if(rb==0)    d= 1.0;
  else              d= 1.0;
    
  return d;
  }

/*-------------------------------------------------------------------------
 * real die metric p(x,y) where x and y are fdie sequences computed as
 * p(x,y) = d(x0,y0) + d(x1,y1) + d(x2,y2) + ... + d(x{N-1},y{N-1})
 * where d(a,b) is defined above.
 * On success return d(x,y). On error return -1.
 *-------------------------------------------------------------------------*/
double fdie_metric(fdieseq x, fdieseq y){
  double rval,d;
  long n;
  long N=x.getN();
  long M=y.getN();
  long NM=(N<M)?N:M; //NM = the smaller of N and M
  for(n=0,d=0;n<NM;n++){
    rval=fdie_metric(x.get(n),y.get(n));
    if(rval<0){d+=0.0; printf("rval=%lf ",rval);}
    else d+=rval;
    }
  if(N!=M){
    fprintf(stderr,"ERROR using fdie_metric(x,y): size of x (%ld) does not equal the size of y (%ld).\n",N,M);
    exit(EXIT_FAILURE);
    }
  return d;
  }



