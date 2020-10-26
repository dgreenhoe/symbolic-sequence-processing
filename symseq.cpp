/*============================================================================
 * Daniel J. Greenhoe
 * routines for Real Die symseqs
 *============================================================================*/
/*=====================================
 * headers
 *=====================================*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "symseq.h"

/*-------------------------------------------------------------------------
 * constructor initializing symseq to '.'
 *-------------------------------------------------------------------------*/
symseq::symseq(long M){
  long n;
  void *memptr;
  N=M;
  memptr=malloc(N*sizeof(char));
  if(memptr==NULL){
    fprintf(stderr,"symseq::symseq memory allocation for %ld elements failed\n",M);
    exit(EXIT_FAILURE);
    }
  x = (char *)memptr;
  for(n=0; n<N; n++)x[n]='.';
  }

/*-------------------------------------------------------------------------
 * constructor initializing symseq to '.'
 *-------------------------------------------------------------------------*/
symseq::symseq(const long M,const unsigned seed,const char *symbols){
  void *memptr;
  N=M;
  memptr=malloc(N*sizeof(char));
  if(memptr==NULL){
    fprintf(stderr,"symseq::symseq memory allocation for %ld elements failed\n",M);
    exit(EXIT_FAILURE);
    }
  x = (char *)memptr;
  randomize(seed,symbols);
  }

/*-------------------------------------------------------------------------
 * fill the symseq with the value '.'
 *-------------------------------------------------------------------------*/
void symseq::clear(void){
  long n;
  for(n=0; n<N; n++) x[n] = '.';
  }

/*-------------------------------------------------------------------------
 * get a symbol from the symseq x at location n
 *-------------------------------------------------------------------------*/
char symseq::get(const long n) const 
{
  if(n<0 || n>=N){//domain check
    fprintf(stderr,"ERROR using symseq::get(n): n=%ld outside the domain [0:%ld] of the sequence.\n",n,N);
    exit(EXIT_FAILURE);
    }
  return x[n];
}

/*-------------------------------------------------------------------------
 * get a symbol from the symseq x at location n
 * but exit if symbol is not in the string <range>
 * example: symbol=get(n,"ABCDEF");
 *-------------------------------------------------------------------------*/
char symseq::get(const long n, const char *range) const
{
  int match;
  char *sptr;
  char symbol;
  if(n<0 || n>=N){//domain check
    fprintf(stderr,"\nERROR using symseq::get(n): n=%ld outside the domain [0:%ld] of the sequence.\n",n,N);
    exit(EXIT_FAILURE);
    }
  symbol = x[n];
  for(match=0,sptr=(char*)range;*sptr!='\0';sptr++) if(symbol==*sptr) match=1;
  if(!match){//range check
    fprintf(stderr,"\nERROR using symbol=symseq::get(%ld): symbol='%c' (0x%02X) outside the range {%s} of the sequence.\n",n,symbol,symbol,range);
    exit(EXIT_FAILURE);
    }
  return symbol;
}

/*-------------------------------------------------------------------------
 * put a single value from the symseq x at location n
 *-------------------------------------------------------------------------*/
const void symseq::put(const long n, const char symbol){
  if(n<0 || n>=N){//domain check
    fprintf(stderr,"ERROR using symseq::put(n): n=%ld outside the domain [0:%ld] of the sequence.\n",n,N-1);
    exit(EXIT_FAILURE);
    }
  x[n]=symbol;
  }

/*-------------------------------------------------------------------------
 * put the symbol <symbol> into the sequence <x> 
 * from location <start> to location <end>
 *-------------------------------------------------------------------------*/
void symseq::put(const long start, const long end, const char symbol){
  long n;
  for(n=start;n<=end;n++) put(n,symbol);
  }

/*-------------------------------------------------------------------------
 * list contents of dieseq
 *-------------------------------------------------------------------------*/
void symseq::list(const long start, const long end, const char *str1, const char *str2, const int display, FILE *fptr){
  long n,m;
  if(strlen(str1)>0){
    if(display)printf("%s",str1);
    if(fptr!=NULL)fprintf(fptr,"%s",str1);
    }
  for(n=start,m=1; n<=end; n++,m++){
    if(display)printf("%c",get(n));
    if(m%50==0&&display)printf("\n");
    else if(m%10==0&&display)printf(" ");
    if(fptr!=NULL){
      fprintf(fptr,"%c",get(n));
      if(m%50==0)fprintf(fptr,"\n");
      else if(m%10==0)fprintf(fptr," ");
      }
    }
  if(strlen(str2)>0){
    if(display)printf("%s",str2);
    if(fptr!=NULL)fprintf(fptr,"%s",str2);
    }
  }

/*-------------------------------------------------------------------------
 * shift symseq n elements to the right inserting zeros on the left
 * example: if x = [ a b c d e f ] (N=6), then shiftR(2) results in
 *             x = [ 0 0 a b c d ] (N=6).
 *-------------------------------------------------------------------------*/
void symseq::shiftR(long n){
  long m;
  for(m=N-1;m-n>=0;m--)x[m]=x[m-n];
  for(m=0;m<n;m++)     x[m]='.';
  }

/*-------------------------------------------------------------------------
 * shift symseq n elements to the left inserting zeros on the right
 * example: if x = [ a b c d e f ] (N=6), then shiftL(2) results in
 *             x = [ c d e f 0 0 ] (N=6).
 *-------------------------------------------------------------------------*/
void symseq::shiftL(long n){
  long m;
  for(m=0;m<N-n;m++) x[m]=x[m+n];
  for(m=N-2;m<N;m++) x[m]='.';
  }

/*-------------------------------------------------------------------------
 * fill the sequence with uniformly distributed pseudo-random symbols
 * from the string <symbols>
 *-------------------------------------------------------------------------*/
void symseq::randomize(const char *symbols){
  const long N=getN();
  const int  M=strlen(symbols);
  int r,i;
  long n;
  for(n=0; n<N; n++){
    r=rand();
    i = r%M;
    put(n,symbols[i]);
    }
  }

/*=====================================
 * operators
 *=====================================*/
/*-------------------------------------------------------------------------
 * operator symseq x = symseq y
 *-------------------------------------------------------------------------*/
void symseq::operator=(symseq y){
  const long M=y.getN();
  long n;
  if(N!=M){
    fprintf(stderr,"ERROR using symseq x = symseq y operation: size of x (%ld) does not equal size of y (%ld)\n",N,M);
    exit(EXIT_FAILURE);
    }
  for(n=0;n<N;n++)x[n]=y.get(n);
  }

/*-------------------------------------------------------------------------
 * operator symseq x == symseq y
 * compare x and y; return 1 if the same, 0 if different.
 *-------------------------------------------------------------------------*/
int symseq::operator==(symseq y){
  const long M=y.getN();
  long n;
  int retval;
  char xsym,ysym;

  if(N!=M){
    fprintf(stderr,"ERROR using symseq x == symseq y operation: size of x (%ld) does not equal size of y (%ld)\n",N,M);
    exit(EXIT_FAILURE);
    }
  for(n=0,retval=1;n<N;n++){
    xsym=  get(n);
    ysym=y.get(n);
    if(xsym!=ysym)retval=0;
    }
  return retval;
  }

/*=====================================
 * external functions
 *=====================================*/
/*-------------------------------------------------------------------------
 * compare dieseq x and dieseq y
 * return the number of locations in which the two sequences are different
 * return 0 if the same
 *-------------------------------------------------------------------------*/
long cmp(const symseq *x, const symseq *y, int showdiff, FILE *fptr){
  const long N=x->getN();
  const long M=y->getN();
  char xsym,ysym;
  long n;
  long count;
  if(N!=M){
    fprintf(stderr,"\nERROR using cmp(symseq *x,symseq *y): size of x (%ld) != size of y (%ld).\n",N,M);
    exit(EXIT_FAILURE);
    }
  for(n=0,count=0;n<N;n++){
    xsym=x->get(n);
    ysym=y->get(n);
    if(xsym!=ysym){
      count++;
      if(showdiff)  fprintf(stdout,"%6ld: x[%6ld]=%c(0x%02x)  y[%6ld]=%c(0x%02x)\n",count,n,xsym,xsym,n,ysym,ysym);
      if(fptr!=NULL)fprintf(fptr,  "%6ld: x[%6ld]=%c(0x%02x)  y[%6ld]=%c(0x%02x)\n",count,n,xsym,xsym,n,ysym,ysym);
      }
    }
  return count;
  }

/*-------------------------------------------------------------------------
 * copy the sequence <*x> = [x_start ... x_end] 
 * into the sequence <*y> = [y_0     ... y_{N-1}] 
 * where N = end-start+1
 *-------------------------------------------------------------------------*/
void copy(const long start, const long end, const symseq *x, symseq *y){
  const long Nx = x->getN();
  const long Ny = y->getN();
  long n,m;
  double xx;
  y->clear();
  if(end>=Nx){
    fprintf(stderr,"ERROR using copy(start,end,seqR1 *x, seqR1 *y): <end>=%ld is too large.\n",end);
    exit(EXIT_FAILURE);
    }
  if((end-start+1)!=Ny){
    fprintf(stderr,"ERROR using copy(start,seqR1 *x, seqR1 *y): length of [x_start...x_end] (%ld) != length of y(%ld).\n",Nx-start,Ny);
    exit(EXIT_FAILURE);
    }
  for(n=start,m=0;n<=end;n++,m++){
    xx = x->get( n );
    y->put( m, xx );
    }
  }

/*-------------------------------------------------------------------------
 * downsample sequence by a factor of <M>
 * and write to sequence pointed to by <y>
 *-------------------------------------------------------------------------*/
void downsample(int M, symseq *x, symseq *y){
  const long Nx = x->getN();
  const long Ny = y->getN();
  char symbol;
  long n,m;
  if(M<1){//check validity of factor <M>
    fprintf(stderr,"ERROR using symseq::downsample(M,y): factor=%d must be at least 1\n",M);
    exit(EXIT_FAILURE);
    }
  if(Ny != Nx/M){//check validity of the length of output sequence <y>
    fprintf(stderr,"ERROR using symseq::downsample(M,y): length %ld of output sequence y must be N/M = %ld/%d = %ld\n",Ny,Nx,M,Nx/M);
    exit(EXIT_FAILURE);
    }
  for(n=0,m=0; m<Ny; n+=M,m++){
    symbol = x->get(n);
    y->put(m,symbol);
    }
  }