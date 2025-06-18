/*
 * Header file for C/C++ wrapper for Perple_X-meemum.
 * - Perple_X-meemum by Jamie Connolly (see http://www.perplex.ethz.ch/)
 * - This header was created automatically using create_C_header.py.
 * - It is designed for use with the C wrapper written by Lars Kaislaniemi (lars.kaislaniemi@iki.fi)
 */

/* Changes:
 * - 2018-06-19: Automatic header creation implemented for use with Perple_X 6.8.3+ (Bob Myhill)
 * - 2018-03-15: Changes in parameters to accomodate changes in Perple_X 6.8.1 (Julien Cornet)
 * - 2018-03-15: Change in common block cxt15: K19 replaced k5 in the declaration (Julien Cornet)
 * - 2016-10-24: Changes in parameters according to changes in Perple_X 6.7.4 (Julien Cornet)
 * - 2013-06-11: Changes to parameters and common block definitions to accommodate changes in PerpleX 6.6.8.
 * - 2012-03-23: Changed p_k21 from 500000 to 1000000, according to changes in Perple_X on 2012-03-07.
*/

/* PerpleX constant parameters (see explanation below or in Perple_X
 * sources perplex_parameters.h. These MUST HAVE the same values as
 * in Perple_X code (perplex_parameters.h), otherwise resulting
 * in array size mismatch and overall havoc.
 */
#define p_n10 7
#define p_n11 8
#define p_n1 11
#define p_n2 12
#define p_n3 13
#define p_n4 14
#define p_n5 15
#define p_n6 16
#define p_n7 17
#define p_n8 18
#define p_n9 19
#define p_n12 20
#define p_n13 21
#define p_n0 30
#define p_mst 3
#define p_mdim 8
#define p_msp p_mdim+2
#define p_ms1 p_msp-1
#define p_h5 5
#define p_h6 500
#define p_h8 200
#define p_h9 30
#define p_i6 2
#define p_i7 20
#define p_i8 28
#define p_i9 200
#define p_i10 50
#define p_i11 100
#define p_j3 4
#define p_j4 6
#define p_j5 8
#define p_j6 8
#define p_j9 160000
#define p_k0 25
#define p_k1 1500000
#define p_k2 100000
#define p_k3 2000
#define p_k4 32
#define p_k5 12
#define p_k7 p_k5+1
#define p_k8 p_k5+2
#define p_k9 30
#define p_k10 400
#define p_k14 18
#define p_k15 6
#define p_k16 100
#define p_k17 7
#define p_k18 29*p_k1
#define p_k19 3*p_k5
#define p_k21 2000000
#define p_k20 (p_mdim+3)*p_k21
#define p_k22 p_mdim*p_k19
#define p_k23 25
#define p_k13 p_mdim*p_k21
#define p_k24 p_k13*(p_mst-1)
#define p_l2 5
#define p_l3 p_l2+2
#define p_l5 1000
#define p_l6 1000
#define p_l7 2048
#define p_l8 10
#define p_l9 100
#define p_nsp 17
#define p_l10 p_nsp+p_l9+4
#define p_m0 8
#define p_m1 60
#define p_m2 8
#define p_m3 3
#define p_m4 20
#define p_m6 6
#define p_m7 15
#define p_m8 9
#define p_m9 10
#define p_m10 6
#define p_m11 6
#define p_m12 4
#define p_m13 8
#define p_m14 2
#define p_m15 12
#define p_m16 6
#define p_m17 5
#define p_m18 6
#define p_nx 500
#define p_ny 500
#define p_lchar 400
#define p_kd2 p_k8*35

/* Maximum length of the filename passed to ini_phaseq / in
 * cst228.prject. Including trailing \0.
 */
#define p_max_filename_len 100

/* Phase name length. Not defined in PerpleX, but used as 
 * numeric constant 14. We add one to include the 
 * trailing \0 in C 
 */
#define p_pname_len 15

/* Component name length. Not defined in Perplex, but used as
 * numeric constant 5. We add one to include the 
 * trailing \0 in C
 */
#define p_cname_len 6

#define p_size_sysprops p_i8   /* num of sysprops */
#define p_size_phases p_k0     /* (max) num of phases */
#define p_size_components p_k5 /* (max) num of components */


/* PerpleX (Fortran) subroutines that read initial setup
 * and initialize arrays for calculations
 */
void input1_(int *first, int *output, int *err);
void input2_(int *first);
void setau1_(int *output);
void input9_(int *first, int *output);
void initlp_();
/* void vrsion_(); // not used */ 
/* void fopen1_(); // not used */
/* void iniprp_(); // not used */


/* PerpleX subroutines for phase equilibria calculations, 
 * derivative property calculations, and printing out
 * the results (calpr0, can be used for debugging here)
 */
void lpopt0_(int *result);
void getloc_(int *itri, int *jtri, int *ijpt, double *wt, int *nodata);
void calpr0_(int *lu);
/* void rebulk_(int *isstatic); // not used */
/* void reopt_(int *result);    // not used */
/* void yclos1_(double *clamda, double *x, int *is, int *jphct, int *quit); // not used */


/* Functions of the wrapper */
int ini_phaseq(char const *);
int phaseq(double, double, int, double*, int*, double*, double*, double*, char*, int);
void print_comp_order();
int get_comp_order(char **order);
int number_of_components();

/* Variables of the wrapper */
// int lpopt_warmstart;   /* set cst111_.istart to 1 before calling lpopt */




/* The common block definitions of Perple_X. Always keep the order
 * of variables inside one struct the same as in Perple_X. 
 * The total size of struct should match the common block size in Fortran.
 * This may cause troubles between platforms (e.g. integer size 32 / 64 bit?)
 */
 
/* Modified 2013-06-11: Here, in principle, we only need to define those structs
 * that are read/written in the C code. Rest of them are now commented out.
 */

extern struct {
	int iam;
} cst4_;

extern struct {
	int gflu, aflu, fluid[p_k5], shear, lflu, volume, rxn;
} cxt20_;

extern struct {
	double nopt[p_i10];
	int iopt[p_i10];
	int lopt[p_i10];
} opts_;
      
extern struct {
	char prject[100];
	char tfname[100];
} cst228_;

extern struct {
      int icomp,istct,iphct,icp;
} cst6_;

extern struct {
	char xname[p_k5][8];
	char vname[p_l2][8];
} csta2_;

extern struct {
      double cblk[p_k5];
      int jbulk;
} cst300_;

extern struct {
      double atwt[p_k0];
} cst45_;

extern struct {
      double a[p_k1][p_k5], b[p_k5], c[p_k1];
} cst313_;

extern struct {
	int ipot;
	int jv[p_l2];
	int iv[p_l2];
} cst24_;

extern struct {
      double v[p_l2];
      double tr,pr,r,ps;
} cst5_;

extern struct {
	int jphct, istart;
} cst111_;

extern struct {						/* k19 replace k5 in common block cxt15 in ppx681 (JCornet 2018)*/
	double cp3[p_k19][p_k0];
	double amt[p_k19];
	int kkp[p_k19], np, ncpd, ntot;
} cxt15_;

extern struct {
	char pname[p_k5][14];
} cxt21a_;

extern struct {
	double props[p_k5][p_i8], psys[p_i8], psys1[p_i8], pgeo[p_i8], pgeo1[p_i8];
} cxt22_; 

extern struct {
	double pcomp[p_k5][p_k0];
} cst324_;

extern struct {
	char cname[p_k5][5];
} csta4_;

extern struct {
	double gtot, fbulk[p_k0], gtot1, fbulk1[p_k0];
} cxt81_;
