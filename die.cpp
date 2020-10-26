/*============================================================================
 * Daniel J. Greenhoe
 * routines for Real Die dieseqs
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
#include "r6.h"
#include "c1.h"
#include "die.h"

/*=====================================
 * prototypes
 *=====================================*/
void phistogram(seqR1 *data, const long start, const long end, FILE *ptr);

/*-------------------------------------------------------------------------
 * fill the dieseq with weighted pseudo-random die face values
 *-------------------------------------------------------------------------*/
int dieseq::randomize(long start, long end, int wA,int wB,int wC,int wD,int wE,int wF){
  int r,u;
  long n;
  char symbol;
  int sum=wA+wB+wC+wD+wE+wF;
  if(sum!=100){
    fprintf(stderr,"dieseq::randomize error: sum of weight values = %d != 100\n",sum);
    return -1;
    }
  //printf("start=%ld end=%ld  weights=(%03d %03d %03d %03d %03d %03d)\n", start,end,wA,wB,wC,wD,wE,wF);
  for(n=start; n<=end; n++){
    r=rand();
    u = r%100;
    if     (u<wA)             symbol='A';
    else if(u<wA+wB)          symbol='B';
    else if(u<wA+wB+wC)       symbol='C';
    else if(u<wA+wB+wC+wD)    symbol='D';
    else if(u<wA+wB+wC+wD+wE) symbol='E';
    else                      symbol='F';
    put(n,symbol);
    }
  return 0;
  }

/*-------------------------------------------------------------------------
 * map die face values to R^1
 * A-->1  B-->2  C-->3  D-->4  E-->5  F-->6
 * all other values --> 0
 *-------------------------------------------------------------------------*/
seqR1 dieseq::dietoR1(void){
  const long N=getN();
  long n;
  char sym;
  double xR1;
  seqR1 y(N);
  for(n=0; n<N; n++){
    sym = get(n);
    xR1 = die_dietoR1(sym);
    y.put(n,xR1);
    }
  return y;
  }

/*-------------------------------------------------------------------------
 * map die face values to R^1 using PAM scheme
 * A-->-2.5  B-->-1.5  C-->-0.5  D-->0.5  E-->1.5  F-->2.5
 * all other values --> 0
 *-------------------------------------------------------------------------*/
seqR1 dieseq::dietoR1pam(void){
  const long N=getN();
  long n;
  char sym;
  double xR1;
  seqR1 y(N);
  for(n=0; n<N; n++){
    sym = get(n);
    xR1 = die_dietoR1pam(sym);
    y.put(n,xR1);
    }
  return y;
  }

/*-------------------------------------------------------------------------
 * map die face values to C^1
 *-------------------------------------------------------------------------*/
seqC1 dieseq::dietoC1(void){
  const long N=getN();
  long n;
  char sym;
  complex xC1;
  seqC1 y(N);
  for(n=0; n<N; n++){
    sym = get(n);
    xC1 = die_dietoC1c(sym);
    y.put(n,xC1);
    }
  return y;
  }

/*-------------------------------------------------------------------------
 * map die face values to R^3 sequence
 *-------------------------------------------------------------------------*/
seqR3 dieseq::dietoR3(void){
  const long N=getN();
  long n;
  char sym;
  vectR3 xR3;
  seqR3 y(N);
  for(n=0; n<N; n++){
    sym = get(n);
    xR3 = die_dietoR3(sym);
    y.put(n,xR3);
    }
  return y;
  }

/*-------------------------------------------------------------------------
 * compute histogram of dna sequence
 * return seqR1 y of length 6 where 
 * y[1]-->number of dna 'A' symbols, 
 * y[2]-->number of dna 'B' symbols, 
 * y[3]-->number of dna 'C' symbols, 
 * y[4]-->number of dna 'D' symbols, 
 * y[5]-->number of dna 'D' symbols, 
 * y[6]-->number of dna 'D' symbols, 
 * y[0]-->number of all other values
 * y[7]-->total number of symbols y[1],y[2],...,y[6]
 *-------------------------------------------------------------------------*/
seqR1 dieseq::histogram(const long start, const long end, int display, FILE *fptr){
  seqR1 data(8);
  long n;
  long bin;
  char symbol;
  data.clear();
  for(n=start;n<=end;n++){
    symbol=get(n);
    switch(symbol){
      case 'A': bin=1; break;
      case 'B': bin=2; break;
      case 'C': bin=3; break;
      case 'D': bin=4; break;
      case 'E': bin=5; break;
      case 'F': bin=6; break;
      default : bin=0; break;
      }
    if(bin!=0) data.increment(7);
    data.increment(bin);
    }
  if(display)   phistogram(&data,start,end,stdout);
  if(fptr!=NULL)phistogram(&data,start,end,fptr  );
  return data;
  }

/*-------------------------------------------------------------------------
 * print die sequence histogram with data pointed to by <data>
 * to stream pointed to by ptr
 *-------------------------------------------------------------------------*/
void phistogram(seqR1 *data, const long start, const long end, FILE *ptr){
  const long N=end-start+1;
  long bin;
  fprintf(ptr,"\n");
  fprintf(ptr," -------------------------------------------------------------------------\n");
  fprintf(ptr,"|  Histogram for sequence [x_n|n=%7ld-%7ld] (length %7ld)        |\n|",start,end,N);
  for(bin=1;bin<=6;bin++)fprintf(ptr,"         %c",'A'+(char)bin-1);  
  fprintf(ptr,"     extra   |\n|");
  for(bin=1;bin<=6;bin++)fprintf(ptr,"%10.0lf",data->get(bin));  
  fprintf(ptr,"%10.0lf   |\n|",data->get(0)); 
  for(bin=1;bin<=6;bin++)fprintf(ptr," (%6.2lf%%)",data->get(bin)/(double)N*100.0);  
  fprintf(ptr,"  (%6.2lf%%)  |\n",data->get(0)/(double)N*100.0); 
  fprintf(ptr," -------------------------------------------------------------------------\n");
  }


/*=====================================
 * operators
 *=====================================*/
/*-------------------------------------------------------------------------
 * operator dieseq x = dieseq y
 *-------------------------------------------------------------------------*/
void dieseq::operator=(dieseq y){
  const long N=  getN();
  const long M=y.getN();
  long n;
  char symbol;
  if(N!=M){
    fprintf(stderr,"ERROR using dieseq x = dieseq y operation: size of x (%ld) does not equal size of y (%ld)\n",N,M);
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
 * map die face values to R^1
 *-------------------------------------------------------------------------*/
double die_dietoR1(char c){
  double rval;
  switch(c){
    case 'A': rval = 1.0;  break;
    case 'B': rval = 2.0;  break;
    case 'C': rval = 3.0;  break;
    case 'D': rval = 4.0;  break;
    case 'E': rval = 5.0;  break;
    case 'F': rval = 6.0;  break;
    default:  
      fprintf(stderr,"ERROR using die_dietoR1(c): c=%c(0x%x) is not in the valid domain {A,B,C,D,E,F}\n",c,c);
      exit(EXIT_FAILURE);
    }
  return rval;
  }

/*-------------------------------------------------------------------------
 * map die face values to R^1 PAM
 *-------------------------------------------------------------------------*/
double die_dietoR1pam(char c){
  double rval;
  switch(c){
    case 'A': rval = -2.5;  break;
    case 'B': rval = -1.5;  break;
    case 'C': rval = -0.5;  break;
    case 'D': rval =  0.5;  break;
    case 'E': rval =  1.5;  break;
    case 'F': rval =  2.5;  break;
    default:  
      fprintf(stderr,"ERROR using die_dietoR1pam(c): c=%c(0x%x) is not in the valid domain {A,B,C,D,E,F}\n",c,c);
      exit(EXIT_FAILURE);
    }
  return rval;
  }

/*-------------------------------------------------------------------------
 * map die face values to complex plane C^1
 *                                                        
 *                 imaginary axis                         
 *                       |                                
 *                       B=(cos90,sin90)                  
 *                       |                                
 *                       |                                
 * (cos150,sin150)=C     |     A=(cos30,sin30)            
 *                       |                                
 *     ------------------|----------------- real axis     
 *                       |                                
 * (cos210,sin210)=D     |     F=(cos330,sin330)          
 *                       |                                
 *                       |                                
 *                       E=(cos270,sin270)                
 *                       |                                
 *                                                        
 *-------------------------------------------------------------------------*/
complex die_dietoC1c(char c){
  complex rc;
  switch(c){
    case 'A': rc = expi( 30.0/180.0*PI);  break;
    case 'B': rc = expi( 90.0/180.0*PI);  break;
    case 'C': rc = expi(150.0/180.0*PI);  break;
    case 'D': rc = expi(210.0/180.0*PI);  break;
    case 'E': rc = expi(270.0/180.0*PI);  break;
    case 'F': rc = expi(330.0/180.0*PI);  break;
    case '0': rc.put(0,0);                break;
    default:  rc.put(0,0);
    fprintf(stderr,"ERROR using dietoC1(char c): c=%c(0x%x) is not in the valid domain {0,A,B,C,D,E,F}. Returning (0,0).\n",c,c);
    }
  return rc;
  }

/*-------------------------------------------------------------------------
 * map die face values to R^3
 *     |+1|     | 0|     | 0|     | 0|     | 0|     |-1| 
 * A-->| 0| B-->|+1| C-->| 0| D-->| 0| E-->|-1| F-->| 0| 
 *     | 0|     | 0|     |+1|     |-1|     | 0|     | 0| 
 *-------------------------------------------------------------------------*/
vectR3 die_dietoR3(char c){
  vectR3 xyz;
  switch(c){
    case 'A': xyz.put(+1, 0, 0);  break;
    case 'B': xyz.put( 0,+1, 0);  break;
    case 'C': xyz.put( 0, 0,+1);  break;
    case 'D': xyz.put( 0, 0,-1);  break;
    case 'E': xyz.put( 0,-1, 0);  break;
    case 'F': xyz.put(-1, 0, 0);  break;
    default:  
      fprintf(stderr,"ERROR using die_dietoR3(c): c=%c(0x%x) is not in the valid domain {A,B,C,D,E,F}\n",c,c);
      exit(EXIT_FAILURE);
    }
  return xyz;
  }

/*-------------------------------------------------------------------------
 * map die face values to R^6
 * A-->(1,0,0,0,0,0)    D-->(0,0,0,1,0,0)
 * B-->(0,1,0,0,0,0)    E-->(0,0,0,0,1,0)
 * C-->(0,0,1,0,0,0)    F-->(0,0,0,0,0,1)
 * 0-->(0,0,0,0,0,0)    
 * on ERROR return (0,0,0,0,0,0)    
 *-------------------------------------------------------------------------*/
vectR6 die_dietoR6(char c){
  vectR6 rsix;
  switch(c){
    case 'A': rsix.put(1,0,0,0,0,0);  break;
    case 'B': rsix.put(0,1,0,0,0,0);  break;
    case 'C': rsix.put(0,0,1,0,0,0);  break;
    case 'D': rsix.put(0,0,0,1,0,0);  break;
    case 'E': rsix.put(0,0,0,0,1,0);  break;
    case 'F': rsix.put(0,0,0,0,0,1);  break;
    default:  
      fprintf(stderr,"ERROR using dietoR6(char c): c=%c(0x%x) is not in the valid domain {0,A,B,C,D,E,F}.\n",c,c);
      exit(EXIT_FAILURE);
    }
  return rsix;
  }



