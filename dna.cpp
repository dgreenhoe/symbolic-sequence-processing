/*============================================================================
 * Daniel J. Greenhoe
 * routines for Real gsp dnaseqs
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
#include "dna.h"

/*=====================================
 * prototypes
 *=====================================*/
void dna_phistogram(seqR1 *data, const long start, const long end, FILE *ptr);


/*-------------------------------------------------------------------------
 * copy a dnaseq <y> into dnaseq x starting at location <n>
 * and fill any remaining locations with <c>
 *-------------------------------------------------------------------------*/
void dnaseq::put(dnaseq *y, const long n, const char symbol){
  const long N=getN();
  long i,j;
  long M=y->getN();
  if(n>=N){
    fprintf(stderr,"\nERROR using dnaseq:put(y,n,symbol): n=%ld outside sequence domain [0:%ld]\n",n,N-1);
    exit(EXIT_FAILURE);
    }
  if(!dna_domain(symbol)){
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
void dnaseq::put(const long start, const long end, const char symbol){
  const long N=getN();
  long n;
  if(start<0||end>=N||start>end){
    fprintf(stderr,"ERROR using dnaseq::put(%ld,%ld,'%c')\n",start,end,symbol);
    exit(EXIT_FAILURE);
    }
  for(n=start;n<=end;n++) put(n,symbol);
  }

/*-------------------------------------------------------------------------
 * fill the dnaseq with pseudo-random DNA values
 * using seed value <seed> 
 * distributed with the weight values <wA,wT,wC,wG>
 * where each weight value wX in an integer in the closed interval [0,100]
 * and where the sum of the intervals must be 100.
 *-------------------------------------------------------------------------*/
int dnaseq::randomize(long start, long finish, int wA,int wT,int wC,int wG){
  int r,u;
  long n;
  int sum=wA+wT+wC+wG;
  char symbol;
  if(sum!=100){
    fprintf(stderr,"ERROR using dnaseq::randomize(start,finish,wA,wT,wC,wG): sum of weight values = %d != 100\n",sum);
    exit(EXIT_FAILURE);
    }
  for(n=start; n<=finish; n++){
    r=rand();
    u = r%100;
    if     (u<wA)             symbol='A';
    else if(u<wA+wT)          symbol='T';
    else if(u<wA+wT+wC)       symbol='C';
    else                      symbol='G';
    put(n,symbol);
    }
  return 0;
  }

/*-------------------------------------------------------------------------
 * map dna face values to R^1 using PAM scheme
 * A-->-1.5  T-->-0.5  C-->0.5  G-->0.5
 * all other values --> 0
 *-------------------------------------------------------------------------*/
seqR1 dnaseq::dnatoR1pam(void){
  const long N=getN();
  long n;
  char symbol;
  seqR1 seqR1(N);
  for(n=0; n<N; n++){
    symbol=get(n);
    switch(symbol){
      case 'A': seqR1.put(n,-1.5);  break;
      case 'C': seqR1.put(n,-0.5);  break;
      case 'T': seqR1.put(n, 0.5);  break;
      case 'G': seqR1.put(n, 1.5);  break;
      default: 
        fprintf(stderr,"\nERROR using dnaseq:dnatoR1pam(): symbol='%c'=0x%x not in sequence range {A,T,C,G}\n",symbol,symbol);
        exit(EXIT_FAILURE);
      }
    }
  return seqR1;
  }

/*-------------------------------------------------------------------------
 * map dna face values to R^1 using AT/CG binary scheme
 * A-->1  T-->1  C-->-1  G-->-1
 * all other values --> 0
 *-------------------------------------------------------------------------*/
seqR1 dnaseq::dnatoR1bin(void){
  const long N=getN();
  long n;
  char symbol;
  seqR1 seqR1(N);
  for(n=0; n<N; n++){
    symbol=get(n);
    switch(symbol){
      case 'A': seqR1.put(n, 1);  break;
      case 'C': seqR1.put(n,-1);  break;
      case 'T': seqR1.put(n, 1);  break;
      case 'G': seqR1.put(n,-1);  break;
      default : 
        fprintf(stderr,"\nERROR using dnaseq:dnatoR1bin(): symbol='%c'=0x%x not in sequence range {A,T,C,G}\n",symbol,symbol);
        exit(EXIT_FAILURE);
      }
    }
  return seqR1;
  }

/*-------------------------------------------------------------------------
 * map dna face values to C^1
 *-------------------------------------------------------------------------*/
seqC1 dnaseq::dnatoC1(void){
  const long N=getN();
  long n;
  char symbol;
  complex yy;
  seqC1 y(N);
  for(n=0; n<N; n++){
    symbol = get(n);
    yy = dnatoC1c(symbol);
    y.put(n,yy);
    }
  return y;
  }

/*-------------------------------------------------------------------------
 * map dna face values to complex plane C^1
 *                                                        
 *                 imaginary axis                         
 *                       |                                
 *                       |                                
 *                       |                                
 * (cos135,sin135)=C     |     A=(cos45,sin45)            
 *                       |                                
 *     ------------------|----------------- real axis     
 *                       |                                
 * (cos225,sin225)=T     |     G=(cos315,sin315)          
 *                       |                                
 *                       |                                
 *                       |                                
 *                                                        
 *-------------------------------------------------------------------------*/
complex dnatoC1c(char c){
  complex rc;
  switch(c){
    case 'A': rc = expi( 45.0/180.0*PI);  break;
    case 'C': rc = expi(135.0/180.0*PI);  break;
    case 'T': rc = expi(225.0/180.0*PI);  break;
    case 'G': rc = expi(315.0/180.0*PI);  break;
    //
    //case 'A': rc = expi( 45.0/180.0*PI);  break;  //Gilleans 2007 mapping
    //case 'G': rc = expi(135.0/180.0*PI);  break;
    //case 'C': rc = expi(225.0/180.0*PI);  break;
    //case 'T': rc = expi(315.0/180.0*PI);  break;
    //
    //case 'A': rc.put(+1,+1);  break;
    //case 'G': rc.put(-1,+1);  break;
    //case 'C': rc.put(-1,-1);  break;
    //case 'T': rc.put(+1,-1);  break;
    case '0': rc.put( 0, 0);  break;
    default:  rc.put( 0, 0);
    fprintf(stderr,"ERROR using dnatoC1(char c): c=%c(0x%x) is not in the valid domain {0,A,C,T,G}. Returning (0,0).\n",c,c);
    }
  return rc;
  }

/*-------------------------------------------------------------------------
 * map dna values to R^4 sequence
 *-------------------------------------------------------------------------*/
seqR4 dnaseq::dnatoR4(void){
  const long N=getN();
  long n;
  char yc;
  seqR4 seq4(N);
  for(n=0; n<N; n++)seq4.put(n,dnatoR4c(get(n)));
  return seq4;
  }

/*-------------------------------------------------------------------------
 * map dna face values to R^4
 * A-->(1,0,0,0)    T-->(0,0,1,0)
 * C-->(0,1,0,0)    G-->(0,0,0,1)
 * 0-->(0,0,0,0)    
 * on ERROR return (0,0,0,0)    
 *-------------------------------------------------------------------------*/
vectR4 dnatoR4c(char c){
  vectR4 rfour;
  switch(c){
    case '0': rfour.put(0,0,0,0);  break;
    case 'A': rfour.put(1,0,0,0);  break;
    case 'C': rfour.put(0,1,0,0);  break;
    case 'T': rfour.put(0,0,1,0);  break;
    case 'G': rfour.put(0,0,0,1);  break;
    default:  rfour.put(0,0,0,0);
    fprintf(stderr,"ERROR using dnatoR4(char c): c=%c(0x%x) is not in the valid domain {A,C,T,G}. Returning (0,0,0,0).\n",c,c);
    }
  return rfour;
  }



/*-------------------------------------------------------------------------
 * map gsp face values to R^1
 * A-->1  T-->2  C-->3  G-->4
 *-------------------------------------------------------------------------*/
seqR1 dnaseq::dnatoR1(void){
  const long N=getN();
  long n;
  char symbol;
  seqR1 seqR1(N);
  for(n=0; n<N; n++){
    symbol = get(n);
    switch(symbol){
      case 'A': seqR1.put(n,1);  break;
      case 'T': seqR1.put(n,2);  break;
      case 'C': seqR1.put(n,3);  break;
      case 'G': seqR1.put(n,4);  break;
      default:  
        fprintf(stderr,"\nERROR using dnaseq:dnatoR1(): symbol='%c'=0x%x not in sequence range {A,T,C,G}\n",symbol,symbol);
        exit(EXIT_FAILURE);
      }
    }
  return seqR1;
  }

/*-------------------------------------------------------------------------
 * map gsp face values to R^2 sequence
 *-------------------------------------------------------------------------*/
seqR2 dnaseq::dnatoR2(void){
  const long N=getN();
  long n;
  char symbol;
  vectR2 yy;
  seqR2 seqR2(N);
  for(n=0; n<N; n++){
    symbol=get(n);
    yy=dna_dnatoR2(symbol);
    seqR2.put(n,yy);
    }
  return seqR2;
  }

/*-------------------------------------------------------------------------
 * downsample seqR1 by a factor of <factor>
 *-------------------------------------------------------------------------*/
dnaseq dnaseq::downsample(int factor){
  const long N=getN();
  long n,m;
  long M;
  char symbol;
  if(factor<1){
    fprintf(stderr,"\nERROR using dnaseq::downsample: factor=%d must be at least 1\n",factor);
    exit(EXIT_FAILURE);
    }
  M=N/factor;
  dnaseq newseq(M);
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
 * y[0]-->number of all other values
 * y[5]-->total number of symbols y[1],y[2],...,y[5]
 *-------------------------------------------------------------------------*/
seqR1 dnaseq::histogram(const long start, const long end, int display, FILE *fptr){
  seqR1 data(6);
  long n;
  long bin;
  double p;
  int i;
  char symbol;
  FILE *ptr;
  data.clear();
  for(n=start;n<=end;n++){
    symbol=get(n);
    switch(symbol){
      case 'A': bin=1; break;
      case 'T': bin=2; break;
      case 'C': bin=3; break;
      case 'G': bin=4; break;
      default : bin=0; break;
      }
    if(bin!=0) data.increment(5);
    data.increment(bin);
    }
  if(display)   dna_phistogram(&data,start,end,stdout);
  if(fptr!=NULL)dna_phistogram(&data,start,end,fptr  );
  return data;
  }

/*-------------------------------------------------------------------------
 * print DNA histogram with data pointed to by <data>
 * to stream pointed to by ptr
 *-------------------------------------------------------------------------*/
void dna_phistogram(seqR1 *data, const long start, const long end, FILE *ptr){
  const long N=end-start+1;
  long bin;
  fprintf(ptr,"\n");
  fprintf(ptr," -------------------------------------------------------------------------\n");
  fprintf(ptr,"|  Histogram for dna sequence [x_n|n=%7ld-%7ld] (length %7ld)    |\n|",start,end,N);
  fprintf(ptr,"         A         T         C         G        AT        CG     other   |\n|");
  for(bin=1;bin<=4;bin++)fprintf(ptr,"%10.0lf",data->get(bin));  
  fprintf(ptr,"%10.0lf%10.0lf%10.0lf   |\n|",data->get(1)+data->get(2),data->get(3)+data->get(4),data->get(0)); 
  for(bin=1;bin<=4;bin++)fprintf(ptr," (%6.2lf%%)",data->get(bin)/(double)N*100.0);  
  fprintf(ptr," (%6.2lf%%) (%6.2lf%%)  (%6.2lf%%)  |\n",(data->get(1)+data->get(2))/(double)N*100.0,(data->get(3)+data->get(4))/(double)N*100.0,data->get(0)/(double)N*100.0); 
  fprintf(ptr," -------------------------------------------------------------------------\n");
  }

/*=====================================
 * operators
 *=====================================*/
/*-------------------------------------------------------------------------
 * operator dnaseq x = dnaseq y
 *-------------------------------------------------------------------------*/
void dnaseq::operator=(dnaseq y){
  const long N=getN();
  const long M=y.getN();
  long n;
  char symbol;
  if(N!=M){
    fprintf(stderr,"\nERROR using dnaseq::operator=: length of x (%ld) differs from length of y (%ld).\n",N,M);
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
 * map dna symbols to R^1
 * A-->1 T-->2 C-->3 G-->4 0-->0  other-->-1
 *-------------------------------------------------------------------------*/
double dna_dnatoR1(char symbol){
  double r;
  switch(symbol){
    case 'A': r=1.0;  break;
    case 'T': r=2.0;  break;
    case 'C': r=3.0;  break;
    case 'G': r=4.0;  break;
    default : 
    fprintf(stderr,"ERROR using dna_dnatoR1(symbol): symbol='%c' (0x%x) is not in the valid domain {A,T,C,G}\n",symbol,symbol);
    exit(EXIT_FAILURE);
    }
  return r;
  }

/*-------------------------------------------------------------------------
 * map dna symbols to R^2
 *-------------------------------------------------------------------------*/
vectR2 dna_dnatoR2(char symbol){
  vectR2 r;
  switch(symbol){
    case 'A': r.put( 1.0,       0  );  break;
    case 'T': r.put(-1.0,       0  );  break;
    case 'C': r.put( 0,        +1.0);  break;
    case 'G': r.put( 0,        -1.0);  break;
    default : 
    fprintf(stderr,"ERROR using dna_dnatoR2(symbol): symbol='%c' (0x%x) is not in the valid domain {A,T,C,G}\n",symbol,symbol);
    exit(EXIT_FAILURE);
    }
  return r;
  }

/*-------------------------------------------------------------------------
 * map R^2 values to dna values using Euclidean metric
 *   0    A    B    C    D    E   F   A+..+F
 *-------------------------------------------------------------------------*/
dnaseq dna_R2todna_euclid(seqR2 xy){
  long n;
  int m;
  long N=xy.getN();
  double d[5];
  double smallestd;
  char closestface;
  vectR2 p,q[5];
  dnaseq rdna(N);

  //q[0].put(0,0,0);
  q[1]=dna_dnatoR2('A');
  q[2]=dna_dnatoR2('T');
  q[3]=dna_dnatoR2('C');
  q[4]=dna_dnatoR2('G');

  for(n=0; n<N; n++){
    p.put(xy.getx(n),xy.gety(n));
    smallestd=ae_metric(1,p,q[1]);
    closestface='A';
    for(m=2;m<5;m++){
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
          default: fprintf(stderr,"Error in dna_R2todna_larc(seqR2 xy)\n");
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
dnaseq dna_R1todna_euclid(seqR1 xy){
  long n;
  long N=xy.getN();
  char closestface;
  double p;
  dnaseq rgsp(N);

  for(n=0; n<N; n++){
    p = xy.get(n);
    if(p<1.5)       closestface='A';
    else if(p>=3.5) closestface='G';
    else if(p>=2.5) closestface='C';
    else            closestface='T';
    rgsp.put(n,closestface);
    }
  return rgsp;
  }

/*-------------------------------------------------------------------------
 * dna metric d(a,b)
 *    d(a,b)| 0    A    T    C    G  (b)
 *   -------|-------------------------
 *      a= 0| 0    1    1    1    1  
 *      a= A| 1    0    1    1    1  
 *      a= T| 1    1    0    1    1  
 *      a= C| 1    1    1    0    1  
 *      a= G| 1    1    1    1    0  
 * On success return d(a,b). On error return -1.
 *-------------------------------------------------------------------------*/
double dna_metric(char a, char b){
  int ra=dna_dnatoR1(a);
  int rb=dna_dnatoR1(b);
  double d;

  if(ra<0)fprintf(stderr,"a=%c(0x%x) not in domain of gsp metric d(a,b)\n",a,a);
  if(rb<0)fprintf(stderr,"b=%c(0x%x) not in domain of gsp metric d(a,b)\n",b,b);

       if(ra<0)     d=-1.0;
  else if(rb<0)     d=-1.0;
  else if(ra==rb)   d= 0.0;
  else              d= 1.0;
  return d;
  }

/*-------------------------------------------------------------------------
 * real gsp metric p(x,y) where x and y are rgsp sequences computed as
 * p(x,y) = d(x0,y0) + d(x1,y1) + d(x2,y2) + ... + d(x{N-1},y{N-1})
 * where d(a,b) is defined above.
 * On success return d(x,y). On error return -1.
 *-------------------------------------------------------------------------*/
double dna_metric(dnaseq x, dnaseq y){
  double rval,d;
  long n;
  long N=x.getN();
  long M=y.getN();
  long NM=(N<M)?N:M; //NM = the smaller of N and M
  for(n=0,d=0;n<NM;n++){
    rval=dna_metric(x.get(n),y.get(n));
    if(rval<0){d+=0.0; printf("rval=%lf ",rval);}
    else d+=rval;
    }
  if(N!=M){
    fprintf(stderr,"ERROR using dna_metric(x,y): size of x (%ld) does not equal the size of y (%ld)\n",N,M);
    exit(EXIT_FAILURE);
    }
  return d;
  }

/*-------------------------------------------------------------------------
 * autocorrelation Rxx of a dna sequence x with 2N offset
 *-------------------------------------------------------------------------*/
int dnaseq::Rxxo(const seqR1 *rxx, const int showcount){
  const long N=getN();
  int rval;
  rval=Rxx(rxx,showcount);
  rxx->add(2*N);
  return rval;
  }

/*-------------------------------------------------------------------------
 * autocorrelation Rxx of a dna sequence x
 *-------------------------------------------------------------------------*/
int dnaseq::Rxx(const seqR1 *rxx, const int showcount){
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
double dnaseq::Rxx(const long m){
  const long mm=labs(m);
  const long N=getN();
  long n,nmm;
  double d,sum;
  char a,b;
  for(n=0,sum=0;n<(N+mm);n++){
    nmm=n-mm;
    a=(n  <0 || n  >=N)? 0.0 : get(n);
    b=(nmm<0 || nmm>=N)? 0.0 : get(nmm);
    d=(a==0 || b==0)?    1.0 : dna_metric(a,b);
    sum+=d;
    }
  return -sum;
  }

/*-------------------------------------------------------------------------
 * read dnan sequence from FASTA formatted file
 * and return how many symbols are in it.
 * reference: https://www.genomatix.de/online_help/help/sequence_formats.html
 *-------------------------------------------------------------------------*/
long numsym_fasta_file(const char *filename){
  FILE *fptr;
  int bufN;
  long N=0;
  char buffer[1024];

  if(filename==NULL)fptr=stdout;
  else              fptr=fopen(filename,"r");
  if(fptr==NULL){
    fprintf(stderr,"Unable to open file %s for reading.\n",filename);
    return -1;
    }
  while(fgets(buffer,1024,fptr)!=NULL){
    if(buffer[0]=='>'); //printf("description: %s",buffer);
    else{ 
      bufN = strlen(buffer);
      N += bufN-1;
      }
    }
  return N;
  }

/*-------------------------------------------------------------------------
 * read dnan sequence into <x> from FASTA formatted file
 * reference: https://www.genomatix.de/online_help/help/sequence_formats.html
 *-------------------------------------------------------------------------*/
int read_fasta_file(const char *filename, char *description, dnaseq *x){
  FILE *fptr;
  char buffer[1024];
  int bufN,i;
  long n;
  char symbol;

  if(filename==NULL)fptr=stdout;
  else              fptr=fopen(filename,"r");
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
          case 'a': symbol='A'; break;
          case 't': symbol='T'; break;
          case 'c': symbol='C'; break;
          case 'g': symbol='G'; break;
          default:  symbol='x'; fprintf(stderr,"unknown character %c (%02x) in dnaseq\n",buffer[i],buffer[i]);
          }
        x->put(n,symbol);
        n++;
        }
      }
    }
  fclose(fptr);
  return 0;
  }


/*-------------------------------------------------------------------------
 * return 1 if in domain
 * return 0 if not in domain 
 *-------------------------------------------------------------------------*/
int dna_domain(const char c){
  int rval;
  rval=(dna_dnatoR1(c)>=0)? 1 : 0;
  return rval;
  }

