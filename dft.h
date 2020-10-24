/*============================================================================
 * Daniel J. Greenhoe
 * DSP operations 
 *============================================================================*/
extern complex dftn(seqR1 *x, const long n);
extern complex dftn(seqC1 *x, const long n);
extern vectC4  dftn(seqR4 *x, const long n);
extern vectC6  dftn(seqR6 *x, const long n);
extern int dft(seqR1 *x, seqC1 *y);
extern int dft(seqC1 *x, seqC1 *y);
extern int dft(seqR4 *x, seqC4 *y);
extern int dft(seqR6 *x, seqC6 *y);





