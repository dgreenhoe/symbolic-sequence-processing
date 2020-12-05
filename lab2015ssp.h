/*============================================================================
 * Daniel J. Greenhoe
 *============================================================================*/
extern int lab_fdie_ocs (const unsigned seed, const long N, const char *basefilename);
extern int lab_rdie_ocs (const unsigned seed, const long N, const char *basefilename);
extern int lab_spin_ocs (const unsigned seed, const long N, const char *basefilename);
extern int lab_wdie_ocs (const unsigned seed, const long N, const char *basefilename);
extern int lab_wrdie_ocs(const unsigned seed, const long N, const char *basefilename);
extern int lab_wspin_ocs(const unsigned seed, const long N, const char *basefilename);
extern int lab_dna_ocs  (const unsigned seed, const long N, const char *basefilename);
extern int lab_dna_ocs  (const char *datafilename, const char *basefilename);
extern int lab_dnan_ocs (const char *datafilename, const char *basefilename);

extern int lab_rdie_auto(FILE *fptr);
extern int lab_rdie_averaging_auto(FILE *fptr);
extern int lab_rdie_averaging_R3_auto(FILE *fptr);

extern int lab_rdie_auto(const unsigned seed, const long N, const char *filename);
extern int lab_fdie_auto(const unsigned seed, const long N, const char *filename);
extern int lab_wdie_auto(const unsigned seed, const long N, const char *filename);
extern int lab_spin_auto(const unsigned seed, const long N, const char *filename);
extern int lab_wspin_auto(const unsigned seed, const long N, const char *filename);

extern int lab_rdie_averaging_auto(const char *filename);
extern int lab_rdie_averaging_R3_auto(const char *filename);

extern int lab_wdie_seq(const unsigned seed, long start, long finish, int wA,int wB,int wC,int wD,int wE,int wF, const char *filename);
extern int lab_spin_seq(const unsigned seed, const long start, const long finish, const char *filename);
extern int lab_wspin_seq(const unsigned seed, const long start, const long finish, int wA,int wB,int wC,int wD,int wE,int wF, const char *filename);

extern int lab_dna_seq(const long start, const long finish, const char *datafilename, const char *outfilename);
extern int lab_dna_seq(const unsigned seed, const long start, const long finish, const char *filename);
extern int lab_rdie_lp (const unsigned seed, const long N, const long M, const char *filename);
extern int lab_wrdie_hp(const unsigned seed, const long N, const long M, const char *basefilename);
extern int lab_wspin_hp(const unsigned seed, const long N, const long M, const char *basefilename);
extern int lab_spin_lp (const unsigned seed, const long N, const long M, const char *basefilename);
extern int lab_fdie_lp (const unsigned seed, const long N, const long M, const char *basefilename);
extern int lab_wdie_hp(const unsigned seed, const long N, const long M, const char filter, const char *basefilename);
extern int lab_die_nonstat34(const unsigned seed, const long N, const long M, const int vx, const char *basefilename);
extern int lab_fdie_seq(const unsigned seed, const long start, const long finish, const char *filename);
extern int lab_die_edge(const unsigned seed, const long N, const long M, const long Mh, const int vx, const char *basefilename);
extern int lab_dna_dft(const char *datafilename, const char *basefilename);
extern int lab_dna_nonstatCT(const unsigned seed, const long N, const long M, const int vx, const double plotmin, const double plotmax, const char *basefilename);
extern int lab_dna_edge(const unsigned seed, const long N, const long M, const long Mh, const int vx, const char *basefilename);
extern int lab_dna_edge(const long Mh, const char *datafilename, const char *basefilename);
extern int lab_dna_averaging(const long Mh, const char *datafilename, const char *basefilename);


