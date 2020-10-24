/*============================================================================
 * Daniel J. Greenhoe
 *============================================================================*/
extern int plot_balloon_metric_ball(vectR2 p, double r, FILE *fptr);
extern int plot_balloon_metric_ball(double x, double y, double r, FILE *fptr){vectR2 p(x,y); return plot_balloon_metric_ball(p,r,fptr);}
extern int plot_balloon_metric_ball(vectR2 p, double r){return plot_balloon_metric_ball(p,r,stdout);}
extern int plot_balloon_metric_ball(double x, double y, double r){return plot_balloon_metric_ball(x,y,r,stdout);}

extern int plot_larc_ball(vectR2 p, double r, FILE *fptr);
extern int plot_larc_ball(double x, double y, double r, FILE *fptr){vectR2 p(x,y); return plot_larc_ball(p,r,fptr);}
extern int plot_larc_ball(vectR2 p, double r){return plot_larc_ball(p,r,stdout);}
extern int plot_larc_ball(double x, double y, double r){return plot_larc_ball(x,y,r,stdout);}

extern int plot_larc_ball(vectR3 p, double r, FILE *fptr);
extern int plot_larc_ball(double x, double y, double z, double r, FILE *fptr){vectR3 p(x,y,z); return plot_larc_ball(p,r,fptr);}

extern int plot_euclidean_metric_ball(double alpha, vectR3 p, double r, FILE *fptr);
extern int plot_euclidean_metric_ball(double alpha, double x, double y, double z, double r, FILE *fptr){vectR3 p(x,y,z); return plot_euclidean_metric_ball(alpha,p,r,fptr);}

extern int plot_mca_metric_ball(vectR2 p, double r, FILE *fptr);
extern int plot_mca_metric_ball(double x, double y, double r, FILE *fptr){vectR2 p(x,y); return plot_mca_metric_ball(p,r,fptr);}


