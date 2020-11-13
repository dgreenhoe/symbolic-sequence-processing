/*============================================================================
 * Daniel J. Greenhoe
 * elliptic routines
 * "ellipse" here is defined as all the points (x,y) in R^2 that satisfy
 *   x^2    y^2   
 *   ---  + ---  = 1  
 *   a^2    b^2
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
#include "elliptic.h"

/*-------------------------------------------------------------------------
 * normalize a point <p> with respect to the ellipse(a,b,phi,xo,yo).
 * In particular, return a point q that is at the same position relative 
 * to the ellipse(a,b,0,0,0) as <p> is to the ellipse(a,b,phi,xo,yo).
 *-------------------------------------------------------------------------*/
vectR2 ellipsec::normalize(vectR2 p){
  vectR2 xyo(xo,yo);
  p -= xyo;    // remove offset (xo,yo)
  p &= (-phi); // remove rotation <phi>
  //printf("normalize(): p=(%lf,%lf) phi=%lf \n",p.getx(),p.gety(),phi);
  return p;
  }

/*-------------------------------------------------------------------------
 * starting at the point p on ellipse(a,b,phi,xo,yo) at parameter <s>,
 * search in the direction <direction> for the point q on the ellipse at parameter <t>
 * that is a distance <d> from p as measured along the ellipse.
 * <direction> is either +1 (forward search direction) or -1 (reverse search direction)
 * <N> = resolution control (higher N means higher resolution)
 * <errord> = error measurement (lower errord means "better" estimate)
 *-------------------------------------------------------------------------*/
int ellipsec::findt_dfroms(double s, double d, int direction, long int N, double *t, double *errord){
  double ta, tb;
  long int n;
  double delta=PI/N;
  double dd;
  double besttb, smalleste;
  double prm;

  prm = perimeter(N);
  if(prm<2*d) return 0; // if ellipse too small for search, return error
  if(direction==-1) delta *= -1;
  smalleste = prm;
  ta = s;
  tb = ta;
  for(n=0;n<N;n++){
    tb += delta;
    dd = pathlength(ta,tb,N);
    if(fabs(d-dd)<smalleste){besttb = tb; smalleste=fabs(d-dd);}
    }
  *t = besttb;
  *errord = smalleste;
  return 1;
  }

/*-------------------------------------------------------------------------
 * compute a given x, y, and b
 *
 *   x^2    y^2                        x
 *   ---  + ---  = 1  ==>  a = ----------------
 *   a^2    b^2                 sqrt(1-y^2/b^2)     
 *
 * return 1 if successful, 0 if error
 *-------------------------------------------------------------------------*/
int ellipsec::setab_givenxyb(vectR2 p, double bb){
  vectR2 q = normalize(p);
  double x=q.getx();
  double y=q.gety();
  b = bb;
  if(b<=0) return 0;
  if(x==0) return 0;
  if(y==0) return 0;
  if(fabs(y)>=fabs(b)) return 0;
  a = fabs(sqrt((b*b*x*x)/(b*b-y*y)));
  return 1;
  }

int ellipsec::setab_givenxyb(double x, double y, double bb){
  vectR2 p(x,y);
  return setab_givenxyb(p,bb);
  }

int ellipsec::setab_givenxya(vectR2 p, double aa){
  vectR2 q = normalize(p);
  double x=q.getx();
  double y=q.gety();
  //printf("p=(%lf,%lf) q=(%lf,%lf) phi=%lf\n",p.getx(),p.gety(),q.getx(),q.gety(),phi);
  a = aa;
  if(a<=0) return 0;
  if(x==0) return 0;
  if(y==0) return 0;
  if(fabs(x)>=fabs(a)) return 0;
  b = fabs(sqrt((a*a*y*y)/(a*a-x*x)));
  return 1;
  }

int ellipsec::setab_givenxya(double x, double y, double aa){
  vectR2 p(x,y);
  return setab_givenxya(p,aa);
  }

/*-------------------------------------------------------------------------
 * compute parameter t for parametric representation of ellipse(a,b) 
 * for a point p=(x,y) on the ellipse.
 *-------------------------------------------------------------------------*/
int ellipsec::tgivenxy( const double x, const double y, double *t )
{
  double tt,xx,yy;
  vectR2 xy(x,y); // xy = (x,y)
  ellipsec en(a,b,0,0,0); // normalized ellipse(a,b)
  if(a==0) return 0;
  if(b==0) return 0;
  if((x==0)&&(y==0)) return 0; // at origin

  //xy=xy.offset(-xo,-yo);
  //xy=xy.rotate(-phi);
  xy=normalize(xy);
  xx = xy.getx();
  yy = xy.gety();
  *t=2*PI;

  if(xx==0){// on y-axis
    if(yy>0) tt= PI/2;
    else     tt=-PI/2;
    }
  else if(yy==0){// on x-axis
    if(xx>0) tt= 0;
    else     tt= PI;
    }
  else{
    //printf("2. x=%lf y=%lf xx=%.12lf yy=%lf a=%.12lf b=%lf phi=%lf\n",x,y,xx,yy,a,b,phi);
    if(fabs(xx)>fabs(yy)) tt = acos(fabs(xx)/a);
    else                  tt = asin(fabs(yy)/b);
    }

  if(xx<0 && yy>0) tt = PI-tt;   // 2nd quadrant
  if(xx<0 && yy<0) tt = tt-PI;   // 3rd quadrant
  if(xx>0 && yy<0) tt = -tt;     // 4th quadrant

  *t = tt;
  return 1;
}

int ellipsec::tgivenxy(const vectR2 p, double *t)
{
  const double x = p.getx();
  const double y = p.gety();
  return tgivenxy(x,y,t);
}

/*-------------------------------------------------------------------------
 * find the point (x(t),y(t)) on ellipse(a,b,phi,xo,yo) at parameter <t>
 *-------------------------------------------------------------------------*/
vectR2  ellipsec::xy(double t){
  vectR2 p(a*cos(t),b*sin(t)); 
  vectR2 xyo(xo,yo);
  p &= phi; // rotate counter-clockwise by <phi> radians
  p += xyo; // offset by (xo,yo)
  return p;
  }

/*-------------------------------------------------------------------------
 * compute the value x(t) on ellipse(a,b,phi,xo,yo) at parameter <t>
 *-------------------------------------------------------------------------*/
double ellipsec::x(double t){
  vectR2 p(a*cos(t),b*sin(t));   // normalized point value
  vectR2 r(cos(phi), -sin(phi)); // rotation vector
  return (r^p) + xo; 
}

/*-------------------------------------------------------------------------
 * compute the value y(t) on ellipse(a,b,phi,xo,yo) at parameter <t>
 *-------------------------------------------------------------------------*/
double ellipsec::y(double t){
  vectR2 p(a*cos(t),b*sin(t));   // normalized point value
  vectR2 r(sin(phi),  cos(phi)); // rotation vector
  return (r^p) + yo; 
  }

/*-------------------------------------------------------------------------
 * estimate of perimeter due to Ramanujan
 *-------------------------------------------------------------------------*/
double ellipsec::estimate(void){
  return PI*(3*(a+b)-sqrt(10*a*b+3*(a*a+b*b)));
  }

/*-------------------------------------------------------------------------
 * path length of ellipse(a,b,phi,xo,yo)
 *-------------------------------------------------------------------------*/
double ellipsec::pathlength(double ta, double tb, long int N){
  double sum=0;
  double t1=ta, t2;
  const double delta=(tb-ta)/(double)N;
  long int n;
  vectR2 p,q;

  for (n=0; n<N; n++){
    t2=t1+delta;
    p = xy(t1);
    q = xy(t2);
    sum += chordlength(p,q);
    t1=t2;
    }
  return sum;
  }

/*-------------------------------------------------------------------------
 * Balloon metric from <p> to <q>
 *-------------------------------------------------------------------------*/
double metric_balloon(vectR2 p, vectR2 q){
  double rp=p.mag(),  thetap=p.theta();
  double rq=q.mag(),  thetaq=q.theta();
  ellipsec ellipse;
  double tp,tq;
  double d,ds;
  //printf("rp=%lf  rq=%lf  thetap=%lf  thetaq=%lf\n",rp,rq,thetap,thetaq);
  /*----------------------------------------
   * define ellipse centered at (0,0) for metric
   *----------------------------------------*/
  if(rp>=rq){
    ellipse.setphi(thetap);
    if(ellipse.setab_givenxya(q,rp)==0)printf("ERROR computing elliptic parameter b\n");
    tp = 0;
    if(ellipse.tgivenxy(q,&tq)==0)printf("ERROR finding parameter in ellipse\n");
    //printf("p=(%lf,%lf) q=(%lf,%lf) a=%lf b=%lf phi=%lf tp=%lf tq=%lf\n",p.getx(),p.gety(),q.getx(),q.gety(), ellipse.geta(),ellipse.getb(),ellipse.getphi(),tp,tp);
    }
  else{
    ellipse.setphi(thetaq);
    if(ellipse.setab_givenxya(p,rq)==0)printf("ERROR computing elliptic parameter b\n");
    tq=0;
    if(ellipse.tgivenxy(p,&tp)==0)printf("ERROR finding parameter in ellipse\n");
    }
  /*----------------------------------------
   * measure pathlength on ellipse
   *----------------------------------------*/
  //printf("p=(%lf %lf) q=(%lf %lf) ellipse(%lf,%lf,%lf) tp=%lf  tq=%lf\n",p.getx(),p.gety(),q.getx(),q.gety(),ellipse.geta(),ellipse.getb(),ellipse.getphi(),tp,tq);
  d=ellipse.pathlength(tp,tq,1000);
  ds=d*2/PI;
  return ds;
  }




