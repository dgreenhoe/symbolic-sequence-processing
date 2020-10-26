/*============================================================================
 * Daniel J. Greenhoe
 * routines for Real gsp dnanseqs
 *============================================================================*/
/*=====================================
 * headers
 *=====================================*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
#include "dnan.h"

/*=====================================
 * prototypes
 *=====================================*/
void dnan_phistogram(seqR1 *data, const long start, const long end, FILE *ptr);

/*-------------------------------------------------------------------------
 * copy a dnaseq <y> into dnaseq x starting at location <n>
 * and fill any remaining locations with <c>
 *-------------------------------------------------------------------------*/
void dnanseq::put(dnanseq *y, const long n, const char symbol){
  const long N=getN();
  long i,j;
  long M=y->getN();
  if(n>=N){
    fprintf(stderr,"\nERROR using dnaseq:put(y,n,symbol): n=%ld outside sequence domain [0:%ld]\n",n,N-1);
    exit(EXIT_FAILURE);
    }
  if(!dnan_domain(symbol)){
    fprintf(stderr,"\nERROR using dnaseq:put(y,n,symbol): symbol='%c'=0x%x not in sequence range {A,T,C,G}\n",symbol,symbol);
    exit(EXIT_FAILURE);
    }
  else{
    for(i=0;i<n;i++)    put(i,symbol);
    for(j=0;j<M;j++,i++)put(i,y->get(j));
    for(   ;i<N;i++)    put(i,symbol);
    }
  }

/*-------------------------------------------------------------------------
 * put the value <c> into the sequence x from location <start> to <end>
 *-------------------------------------------------------------------------*/
void dnanseq::put(const long start, const long end, const char symbol){
  const long N=getN();
  long n;
  if(start<0||end>=N||start>end){
    fprintf(stderr,"ERROR using dnaseq::put(%ld,%ld,'%c')\n",start,end,symbol);
    exit(EXIT_FAILURE);
    }
  for(n=start;n<=end;n++) put(n,symbol);
  }

/*-------------------------------------------------------------------------
 * fill the dnanseq with pseudo-random DNA values
 * using seed value <seed> 
 * distributed with the weight values <wA,wT,wC,wG,wN>
 * where each weight value wX in an integer in the closed interval [0,100]
 * and where the sum of the intervals must be 100.
 *-------------------------------------------------------------------------*/
int dnanseq::randomize(long start, long finish, int wA,int wT,int wC,int wG,int wN){
  int r,u;
  long n;
  int sum=wA+wT+wC+wG+wN;
  char symbol;
  if(sum!=100){
    fprintf(stderr,"ERROR using dnanseq::randomize(start,finish,wA,wT,wC,wG,wN): sum of weight values = %d != 100\n",sum);
    exit(EXIT_FAILURE);
    }
  for(n=start; n<=finish; n++){
    r=rand();
    u = r%100;
    if     (u<wA)          symbol='A';
    else if(u<wA+wT)       symbol='T';
    else if(u<wA+wT+wC)    symbol='C';
    else if(u<wA+wT+wC+wG) symbol='G';
    else                   symbol='N';
    put(n,symbol);
    }
  return 0;
  }

/*-------------------------------------------------------------------------
 * map gsp face values to R^1
 * A-->1  T-->2  C-->3  G-->4
 * all other values --> 0
 *-------------------------------------------------------------------------*/
seqR1 dnanseq::dnantoR1(void){
  const long N=getN();
  long n;
  char symbol;
  seqR1 seqR1(N);
  for(n=0; n<N; n++){
    symbol=get(n);
    switch(symbol){
      case 'A': seqR1.put(n,1);  break;
      case 'T': seqR1.put(n,2);  break;
      case 'C': seqR1.put(n,3);  break;
      case 'G': seqR1.put(n,4);  break;
      case 'N': seqR1.put(n,0);  break;
      default:  
        fprintf(stderr,"\nERROR using dnaseq:dnatoR1(): symbol='%c'=0x%x not in sequence range {A,T,C,G,N}\n",symbol,symbol);
        exit(EXIT_FAILURE);
      }
    }
  return seqR1;
  }

/*-------------------------------------------------------------------------
 * map gsp face values to R^2 sequence
 *-------------------------------------------------------------------------*/
seqR2 dnanseq::dnantoR2(void){
  const long N=getN();
  long n;
  char dnasym;
  vectR2 r2sym;
  seqR2 seqR2(N);
  for(n=0; n<N; n++){
    dnasym=get(n);
    r2sym =dnan_dnatoR2(dnasym);
    seqR2.put(n,r2sym);
    }
  return seqR2;
  }

/*-------------------------------------------------------------------------
 * downsample seqR1 by a factor of <factor>
 *-------------------------------------------------------------------------*/
dnanseq dnanseq::downsample(int factor){
  const long N=getN();
  long n,m;
  long M;
  char symbol;
  if(factor<1){
    fprintf(stderr,"\nERROR using dnanseq::downsample: factor=%d must be at least 1\n",factor);
    exit(EXIT_FAILURE);
    }
  M=N/factor;
  dnanseq newseq(M);
  for(n=0,m=0; m<M; n+=factor,m++){
    symbol=get(n);
    newseq.put(m,symbol);
    }
  return newseq;
  }

/*-------------------------------------------------------------------------
 * compute histogram of dna sequence
 * return seqR1 y of length 6 where 
 * y[1]-->number of dna 'A' symbols, 
 * y[2]-->number of dna 'T' symbols, 
 * y[3]-->number of dna 'C' symbols, 
 * y[4]-->number of dna 'G' symbols, 
 * y[5]-->number of dna 'N' symbols, 
 * y[0]-->number of all other values
 * y[6]-->total number of symbols y[1],y[2],...,y[5]
 *-------------------------------------------------------------------------*/
seqR1 dnanseq::histogram(const long start, const long end, int display, FILE *fptr){
  seqR1 data(7);
  long n;
  long bin;
  char symbol;
  data.clear();
  for(n=start;n<=end;n++){
    symbol=get(n);
    switch(symbol){
      case 'A': bin=1; break;
      case 'T': bin=2; break;
      case 'C': bin=3; break;
      case 'G': bin=4; break;
      case 'N': bin=5; break;
      default : bin=0; break;
      }
    if(bin!=0) data.increment(6);
    data.increment(bin);
    }
  if(display)   dnan_phistogram(&data,start,end,stdout);
  if(fptr!=NULL)dnan_phistogram(&data,start,end,fptr  );
  return data;
  }

/*-------------------------------------------------------------------------
 * print DNA histogram with data pointed to by <data>
 * to stream pointed to by ptr
 *-------------------------------------------------------------------------*/
void dnan_phistogram(seqR1 *data, const long start, const long end, FILE *ptr){
  long N=end-start+1;
  long bin;
  fprintf(ptr,"\n");
  fprintf(ptr," -------------------------------------------------------------------------\n");
  fprintf(ptr,"|  Histogram for dnan sequence [x_n|n=%7ld-%7ld] (length %7ld)    |\n|",start,end,N);
  fprintf(ptr,"         A         T         C         G        AT        CG       N     |\n|");
  for(bin=1;bin<=4;bin++)fprintf(ptr,"%10.0lf",data->get(bin));  
  fprintf(ptr,"%10.0lf%10.0lf%10.0lf   |\n|",data->get(1)+data->get(2),data->get(3)+data->get(4),data->get(5)); 
  for(bin=1;bin<=4;bin++)fprintf(ptr," (%6.2lf%%)",data->get(bin)/(double)N*100.0);  
  fprintf(ptr," (%6.2lf%%) (%6.2lf%%)  (%6.2lf%%)  |\n",(data->get(1)+data->get(2))/(double)N*100.0,(data->get(3)+data->get(4))/(double)N*100.0,data->get(5)/(double)N*100.0); 
  fprintf(ptr," -------------------------------------------------------------------------\n");
  }
/*=====================================
 * operators
 *=====================================*/
/*-------------------------------------------------------------------------
 * operator dnanseq x = dnanseq y
 *-------------------------------------------------------------------------*/
void dnanseq::operator=(dnanseq y){
  const long N=getN();
  const long M=y.getN();
  long n;
  char symbol;
  if(N!=M){
    fprintf(stderr,"\nERROR using dnanseq::operator=: length of x (%ld) differs from length of y (%ld).\n",N,M);
    exit(EXIT_FAILURE);
    }
  for(n=0;n<N;n++){
    symbol = y.get(n);
    put(n,symbol);
    }
  }

/*=====================================
 * external operations
 *=====================================*/
/*-------------------------------------------------------------------------
 * return 1 if in domain
 * return 0 if not in domain 
 *-------------------------------------------------------------------------*/
int dnan_domain(const char c){
  int rval;
  rval=(dnan_dnatoR1(c)>=0)? 1 : 0;
  return rval;
  }

/*-------------------------------------------------------------------------
 * map gsp face values to R^1
 * A-->1 T-->2 C-->3 G-->4 N-->5 0-->0  other-->-1
 *-------------------------------------------------------------------------*/
double dnan_dnatoR1(const char symbol){
  double r;
  switch(symbol){
    case 'N': r=0.0;  break;
    case 'A': r=1.0;  break;
    case 'T': r=2.0;  break;
    case 'C': r=3.0;  break;
    case 'G': r=4.0;  break;
    default : 
    fprintf(stderr,"ERROR using dnan_dnatoR1(symbol): symbol='%c' (0x%x) is not in the valid domain {A,T,C,G,N}\n",symbol,symbol);
    exit(EXIT_FAILURE);
    }
  return r;
  }

/*-------------------------------------------------------------------------
 * map dna values to R^2
 *-------------------------------------------------------------------------*/
vectR2 dnan_dnatoR2(const char c){
  vectR2 xy;
  switch(c){
    case '0': xy.put( 0, 0);  break;
    case 'A': xy.put( 1.0,       0  );  break;
    case 'T': xy.put(-1.0,       0  );  break;
    case 'C': xy.put( 0,        +1.0);  break;
    case 'G': xy.put( 0,        -1.0);  break;
    case 'N': xy.put( 0,           0);  break;
    default: 
      fprintf(stderr,"ERROR: c='%c' (0x%x) is not in the valid domain {0,A,T,C,G,N} in dnan_dnatoR2(char c)\n",c,c);
      exit(EXIT_FAILURE);
    }
  return xy;
  }

/*-------------------------------------------------------------------------
 * map R^2 values to dna values using Euclidean metric
 *   0    A    B    C    D    E   F   A+..+F
 *-------------------------------------------------------------------------*/
dnanseq dnan_R2todnan_euclid(seqR2 xy){
  long n;
  int m;
  long N=xy.getN();
  double d[5];
  double smallestd;
  char closestface;
  vectR2 p,q[6];
  dnanseq rdna(N);

  //q[0].put(0,0,0);
  q[1]=dnan_dnatoR2('A');
  q[2]=dnan_dnatoR2('T');
  q[3]=dnan_dnatoR2('C');
  q[4]=dnan_dnatoR2('G');
  q[5]=dnan_dnatoR2('N');

  for(n=0; n<N; n++){
    p.put(xy.getx(n),xy.gety(n));
    smallestd=ae_metric(1,p,q[1]);
    closestface='A';
    for(m=2;m<=5;m++){
      d[m] = ae_metric(1,p,q[m]);
      if(((m&0x01) && (d[m]<smallestd)) || ((!(m&0x01)) && (d[m]<=smallestd))){
        //----------------------------     ----------------------------------
        // bias odd samples                  bias even samples
        // towards smaller values            towards larger values
        smallestd=d[m];
        switch(m){
          case 1: closestface = 'A'; break;
          case 2: closestface = 'T'; break;
          case 3: closestface = 'C'; break;
          case 4: closestface = 'G'; break;
          case 5: closestface = 'N'; break;
          default: fprintf(stderr,"Error in dnan_R2todnan_larc(seqR2 xy)\n");
          }
        }
      }
    rdna.put(n,closestface);
    }
  return rdna;
  }

/*-------------------------------------------------------------------------
 * map R^1 values to gsp face values using Euclidean metric
 *-------------------------------------------------------------------------*/
dnanseq dnan_R1todnan_euclid(seqR1 xy){
  long n;
  long N=xy.getN();
  char closestface;
  double p;
  dnanseq rgsp(N);

  for(n=0; n<N; n++){
    p = xy.get(n);
    if(p<0.5)       closestface='N';
    else if(p<1.5)  closestface='A';
    else if(p>=3.5) closestface='G';
    else if(p>=2.5) closestface='C';
    else            closestface='T';
    rgsp.put(n,closestface);
    }
  return rgsp;
  }

/*-------------------------------------------------------------------------
 * dnan metric d(a,b)
 *    d(a,b)| 0    A    T    C    G    N  (b)
 *   -------|------------------------------
 *      a= 0| 0    1    1    1    1    1  
 *      a= A| 1    0    1    1    1    s  
 *      a= T| 1    1    0    1    1    s  
 *      a= C| 1    1    1    0    1    s  
 *      a= G| 1    1    1    1    0    s  
 *      a= N| 1    s    s    s    s    0  
 * where s = 1/sqrt(2) = 0.70710678118654752440084436210485...
 * On success return d(a,b). On error return -1.
 *-------------------------------------------------------------------------*/
double dnan_metric(const char a, const char b){
  const int ra=dnan_dnatoR1(a);
  const int rb=dnan_dnatoR1(b);
  const double s=1.0/sqrt(2.0);
  double d;

  if(ra<0)fprintf(stderr,"a=%c(0x%x) not in domain of gsp metric d(a,b)\n",a,a);
  if(rb<0)fprintf(stderr,"b=%c(0x%x) not in domain of gsp metric d(a,b)\n",b,b);

       if(ra<0)     d=-1;
  else if(rb<0)     d=-1;
  else if(ra==rb)   d= 0;
  else if(ra==0)    d= 1;
  else if(rb==0)    d= 1;
  else if(ra==5)    d= s;
  else if(rb==5)    d= s;
  else              d= 1;
  return d;
  }

/*-------------------------------------------------------------------------
 * real gsp metric p(x,y) where x and y are rgsp sequences computed as
 * p(x,y) = d(x0,y0) + d(x1,y1) + d(x2,y2) + ... + d(x{N-1},y{N-1})
 * where d(a,b) is defined above.
 * On success return d(x,y). On error return -1.
 *-------------------------------------------------------------------------*/
double dnan_metric(dnanseq x, dnanseq y)
{
  double rval,d;
  long n;
  long N=x.getN();
  long M=y.getN();
  if(N!=M){
    fprintf(stderr,"\nERROR using dnan_metric(x,y): size of x (%ld) does not equal the size of y (%ld)\n",N,M);
    exit(EXIT_FAILURE);
    }
  for(n=0,d=0;n<N;n++){
    rval=dnan_metric(x.get(n),y.get(n));
    if(rval<0){d+=0.0; printf("rval=%lf ",rval);}
    else d+=rval;
    }
  return d;
}

/*-------------------------------------------------------------------------
 * autocorrelation Rxx of a dna sequence x with 2N offset
 *-------------------------------------------------------------------------*/
int dnanseq::Rxxo(seqR1 *rxx, const int showcount) const
{
  const long N=getN();
  int rval;
  rval=Rxx(rxx,showcount);
  rxx->add(2*N);
  return rval;
}

/*-------------------------------------------------------------------------
 * autocorrelation Rxx of a dna sequence x
 *-------------------------------------------------------------------------*/
int dnanseq::Rxx(seqR1 *rxx, const int showcount) const
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
double dnanseq::Rxx(const long m) const
{
  const long mm=labs(m);
  const long N=getN();
  long n,nmm;
  double d,sum;
  char a,b;
  for(n=0,sum=0;n<(N+mm);n++){
    nmm=n-mm;
    a=(n  <0 || n  >=N)? 0.0 : get(n);
    b=(nmm<0 || nmm>=N)? 0.0 : get(nmm);
    d=(a==0 || b==0)?    1.0 : dnan_metric(a,b);
    sum+=d;
    }
  return -sum;
}

/*-------------------------------------------------------------------------
 * read dna seqR1 from file in FASTA format
 * reference: https://www.genomatix.de/online_help/help/sequence_formats.html
 *-------------------------------------------------------------------------*/
int read_fasta_file(const char *filename, char *description, dnanseq *x){
  FILE *fptr;
  int bufN,i;
  long n;
  char symbol;
  char buffer[1024];

  if(filename==NULL)fptr=stdout;
  else         fptr=fopen(filename,"r");
  if(fptr==NULL){
    fprintf(stderr,"\nERROR using read_fasta_file(%s,...): unable to open file.\n",filename);
    exit(EXIT_FAILURE);
    }
  n=0;
  sprintf(description,"No description line found in FASTA file %s.",filename);//default description
  while(fgets(buffer,1024,fptr)!=NULL){
    if(buffer[0]=='>')strcpy(description,buffer);
    else{
      bufN = strlen(buffer);
      for(i=0;i<bufN-1;i++){
        switch(buffer[i]){
          case 'A': symbol='A'; break;
          case 'T': symbol='T'; break;
          case 'C': symbol='C'; break;
          case 'G': symbol='G'; break;
          case 'N': symbol='N'; break;
          case 'a': symbol='A'; break;
          case 't': symbol='T'; break;
          case 'c': symbol='C'; break;
          case 'g': symbol='G'; break;
          case 'n': symbol='N'; break;
          default:  symbol='x'; fprintf(stderr,"unknown character %c (%02x) in dnanseq\n",buffer[i],buffer[i]);
          }
        x->put(n,symbol);
        n++;
        }
      }
    }
  fclose(fptr);
  return 0;
  }


