#include <stdio.h>
#include <math.h>

/*
 * This program solves numerically the RG equations for all
 * LV operators \eta, \xi and \kappa. We thus develop those
 * parameters from some high scale above the EW breaking scale
 * (where we suppose new physics should arise) down to the
 * EW breaking scale. We assume that in this range of energies
 * all particles of the Standard Model come to play, thus the
 * equations are written in the EW-basis. 
 *
 * We have the following LV parameters:
 *       \eta_L, \eta_Q, \eta_e, \eta_u, \eta_d
 *       \xi', \xi, \xi_3, \kappa
 *
 * All \eta's are in fact matrices in flavor space.
 *
 */


/* This is the high scale we go down from */
#define HIGH_SCALE  1.0e+15     /* GeV */

/* We take this as the scale of EW symmetry breaking */
#define EW_SCALE    100.0       /* GeV */

typedef struct 
{
    double r,i;
} tcomplex;
typedef tcomplex *pcomplex;

/*
 * Now, we'll need some notion of matrices (in the flavor space).
 */
#define N_FLAVORS   3                /* total number of flavors */
typedef tcomplex  tmatrix[N_FLAVORS][N_FLAVORS];
typedef tmatrix  *pmatrix;

/**************************************************************
 *             Physical quantities and parameters
 **************************************************************/

/*
 * Current scale: either \mu or t = log \mu
 */
double   mu;
double   t;

/*
 * Running SM couplings
 * Should be recomputed for each new scale
 */
double   g1, g, g3;                        /* for U(1), SU(2), SU(3) */
tmatrix  tlambda_e, tlambda_u, tlambda_d;  /* Yukawa */
pmatrix  lambda_e = &tlambda_e, 
         lambda_u = &tlambda_u, 
         lambda_d = &tlambda_d;

/*
 * Running LV parameters
 */
tmatrix teta_L, teta_Q, teta_e, teta_u, teta_d;
pmatrix eta_L = &teta_L,
        eta_Q = &teta_Q,
        eta_e = &teta_e,
        eta_u = &teta_u,
        eta_d = &teta_d;

/**************************************************************
 * Basic operations we'd like to do with the complex numbers
 **************************************************************/
/**
 * Sets a complex number to a specific value
 */
void
cset(pcomplex num, double re, double im)
{
    num->i = im;
    num->r = re;
}

/**
 * Sets a complex number zero
 */
void 
czero(pcomplex num)
{
    memset(num, 0, sizeof(tcomplex));
}

/**
 * Copies a complex number
 */
ccpy(pcomplex dst, pcomplex src)
{
    memcpy(dst, src, sizeof(tcomplex));
}

/**
 * Prints a complex number
 */
cprint(pcomplex num)
{
    printf("{%g,%g}", num->r, num->i);
}

/**
 * Adds two complex numbers
 */
void
cadd(pcomplex c1, pcomplex c2, pcomplex dst)
{
    dst->r = c1->r + c2->r;
    dst->i = c1->i + c2->i;
}

/**
 * Multiplies two complex numbers
 */
void
cmult(pcomplex c1, pcomplex c2, pcomplex dst)
{
    double tr, ti;

    tr = c1->r * c2->r  -  c1->i * c2->i;
    ti = c1->r * c2->i  +  c1->i * c2->r;

    dst->r = tr;
    dst->i = ti;
}

/**************************************************************
 * And basic operations we'd like to do with the matrices 
 **************************************************************/
/**
 * Gets access to a specific element in the complex matrix
 * as a "pcomplex".
 *  
 *    m  -- assumed to be pmatrix
 */
#define ELT(m, i, k) (&((*m)[(i)][(k)]))

/**
 * Sets a matrix to zero
 */
void
mzero(pmatrix src)
{
    memset(src, 0, sizeof(tmatrix));
}

/**
 * Copies one matrix to another
 */
void
mcpy(pmatrix dst, pmatrix src)
{
    memcpy(dst, src, sizeof(tmatrix));
}

/**
 * Prints a matrix
 */
void
mprint(pmatrix src)
{
    int i, k;

    for (i = 0; i < N_FLAVORS; i++)
    {
	printf("( ");
	for (k = 0; k < N_FLAVORS; k++)
	{
	    printf("{%g,%g} ",(*src)[i][k].r, (*src)[i][k].i);
	}
	printf(" )\n");
    }
}

/** 
 * This is hermitean conjugate of a matrix
 */
void 
mhconj(pmatrix src, pmatrix dst)
{
    int i,k;

    for (i = 0; i < N_FLAVORS; i++)
    {
	for (k = 0; k < N_FLAVORS; k++)
	{
	    (*dst)[k][i].r = (*src)[i][k].r;
	    (*dst)[k][i].i = - (*src)[i][k].i;
	}
    }
}

/**
 * This is the sum of two matrices
 */
void
madd(pmatrix src1, pmatrix src2, pmatrix dst)
{
    int i, k;

    for (i = 0; i < N_FLAVORS; i++)
    {
	for (k = 0; k < N_FLAVORS; k++)
	{
	    cadd(ELT(src1, i, k), ELT(src2, i, k), ELT(dst, i, k));
	}
    }
}

/**
 * This is the product of two matrices
 */
void
mmult(pmatrix src1, pmatrix src2, pmatrix dst)
{
    int      i, k, l;
    tcomplex tmp, tmp1;

    for (i = 0; i < N_FLAVORS; i++)
    {
	for (k = 0; k < N_FLAVORS; k++)
	{
	    czero(&tmp1);
	    for (l = 0; l < N_FLAVORS; l++)
	    {
		cmult(ELT(src1,i,l), ELT(src2,l,k), &tmp);
		cadd(&tmp1, &tmp, &tmp1);
	    }
	    ccpy(ELT(dst, i, k), &tmp1);
	}
    }
}

int 
main(int argc, char **argv)
{
  pmatrix p, ptr, p1;
  tmatrix a, a1;
  tcomplex c1, c2, c3;
  int     t1, t2;
  tmatrix  m1 = { { {1,1},{1,2},{1,3} },
		  { {2,1},{2,2},{2,3} },
		  { {3,1},{3,2},{3,3} } };
  tmatrix  m2 = { { {5,1},{6,2},{7,3} },
		  { {5,1},{6,2},{7,3} },
		  { {5,1},{6,2},{7,3} } };
  tmatrix  m3;

  /* 
   * the result of the product of m1 and m2 is
   *        [9+33*I,  6+42*I,  3+51*I ] 
   *        [24+36*I, 24+48*I, 24+60*I]
   *        [39+39*I, 42+54*I, 45+69*I]
   * 
   * the result of the sum of m1 and m2 is
   *        [6+2*I, 7+4*I, 8+6*I ] 
   *        [7+2*I, 8+4*I, 9+6*I ] 
   *        [8+2*I, 9+4*I, 10+6*I]
   */
  printf("---------- sizeof test ---------------------\n");
  printf("size of tcomplex = %d\n", sizeof(tcomplex));
  printf("size of tmatrix = %d\n", sizeof(tmatrix));
  printf("size of *pmatrix = %d\n", sizeof(*p));


  printf("---------- complex product test ------------\n");
  cset(&c1, 1.0, 2.5);
  cset(&c2, 3.0, 4);  
  cmult(&c1, &c2, &c3);  /* the result should be "-7.0+11.5*I" */
  printf("the product of "); cprint(&c1); printf(" and ");
  cprint(&c2); printf("is: \n");
  cprint(&c3); printf("\n");

  printf("---------- matrix fillup test --------------\n");
  mprint(&m1);

  printf("---------- matrix conjugate test -----------\n");
  ptr = &a;
  for (t1 = 0; t1 < N_FLAVORS; t1++)
  {
      for (t2 = 0; t2 < N_FLAVORS; t2++)
      {
	  cset(ELT(ptr, t1, t2), (double)(t1+1), (double)(t2+1));
      }
  }
  printf("the original matrix:\n");
  mprint(ptr);
  p1 = &a1;
  mhconj(ptr, p1);
  printf("the hermitean conjugated matrix:\n");
  mprint(p1);
  printf("---------- matrix sum test ---------------\n");
  madd(&m1, &m2, &m3);
  printf("the sum of m1 + m2 is\n");
  mprint(&m3);
  printf("---------- matrix product test -----------\n");
  mmult(&m1, &m2, &m3);
  printf("the product of m1 * m2 is\n");
  mprint(&m3);
  
}
