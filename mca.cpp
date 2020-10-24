/*============================================================================
 * Daniel J. Greenhoe
 * routines for Mean Circular arcs
 *         y
 *         |   o p        Let (rp,tp) be the polar location of point p.
 *         |  /           where rp is the Euclidean distance from (0,0) to p 
 *         | /            and tp is radian measure from the x-axis to p.
 *         |/tp           Let (rq,tq) be the polar location of point q.
 * --------|---------- x  The "Mean Circular arc" r(theta) is defined here as
 *         |\tq                      rp + rq
 *         | \            r(theta) = -------  for  tq <= theta <= tp
 *         |  o q                       2
 *
 * note: this metric is *deprecated* because it has been found to be not so 
 *       "useful" for "practical" applications.
 * For example, the Mean Circular Arc metric from (0,1/2) to (1,100) is
 *   99.505 for the Euclidean metric
 *   63.347 for the 2/pi scaled Euclidean metric
 *   63.348 for the Lagrange Arc metric  BUT only
 *    0.319907 for the Mean Circular Arc metric
 *============================================================================*/
/*=====================================
 * headers
 *=====================================*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<main.h>
#include<r1.h>
#include<r2.h>
#include<r3.h>
#include<mca.h>

/*-------------------------------------------------------------------------
 * path length s of Mean Circular Arc from a point p at polar coordinate (rp,tp)
 *                                   to point q at polar coordinate (rq,tq).
 *          rp+rq
 *    s =  ------- tdiff
 *            2
 *-------------------------------------------------------------------------*/
double mca_arclength(double rp, double rq, double tdiff){
  if(rp==0)    return -1; 
  if(rq==0)    return -2; 
  if(tdiff==0) return -3;
  if(tdiff>PI) tdiff = 2*PI-tdiff;
  double y = (rp+rq)/2*tdiff;
  return y;
  }

/*-------------------------------------------------------------------------
 * Mean Circular arc metric from <p> to <q> in R^2
 *-------------------------------------------------------------------------*/
double mca_metric(vectR2 p, vectR2 q){
  double rp=p.mag(), rq=q.mag();
  double tdiff = pqtheta(p,q);
  vectR2  pq=p-q;
  double d;
  if(rp==0 || rq==0 || tdiff<=0) d = pq.mag();
  else if(rp==rq)                d = rp*tdiff;
  else                           d = mca_arclength(rp, rq, tdiff);
  //printf("rp=%lf rq=%lf tdiff=%lf d=%lf ds=%lf\n",rp,rq,tdiff,d,d*2/PI);
  return d*2/PI;
  }

/*-------------------------------------------------------------------------
 * Mean Circular arc metric from <p> to <q> in R^3
 *-------------------------------------------------------------------------*/
double mca_metric(vectR3 p, vectR3 q){
  double rp=p.mag(), rq=q.mag();
  double tdiff = pqtheta(p,q);
  vectR3  pq=p-q;
  double d;
  if(rp==0 || rq==0 || tdiff<=0) d = pq.mag();
  else if(rp==rq)                d = rp*tdiff;
  else                           d = mca_arclength(rp, rq, tdiff);
  return d*2/PI;
  }

/*-------------------------------------------------------------------------
 * Find the polar length of a point q with radial measure tq that is a 
 * distance <d> from the point <p> with polar coordinates (rp,tp)
 * using search resolution <N>
 *-------------------------------------------------------------------------*/
vectR3 mca_findq(vectR3 p, double theta, double phi, double d, long int N){
  double smallesterror=100000,rq,dd,errord;
  vectR3 bestq;
  vectR3 q;
  vectR3 pq;

  for(rq=0; rq<=5; rq+=5.0/N){
    q.polartoxyz(rq,theta,phi);
    dd=mca_metric(p,q);
    errord=fabs(d-dd);
    if(errord<smallesterror){
      bestq=q;
      smallesterror=errord;
      }
    }
  return bestq;
  }

