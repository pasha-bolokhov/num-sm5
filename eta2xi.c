#include <stdio.h>
#include <math.h>


/* 
 *  The mixing matrix is N:
 *
 *          / A  B \
 *     N = ||      ||;
 *          \ C  D /
 */
double  A, B, C, D;
/* 
 * This is the eigenvalues of the matrix N
 */
double  kappa1, kappa2;
/*
 * This is the eigenvectors of the matrix N
 */
double  X1, X2;    /* first eigen vector -- corr. to kappa1  */
double  Y1, Y2;    /* second eigen vector -- corr. to kappa2 */

double  otherX1, otherX2;   /* for checking the eigenvectors */
double  otherY1, otherY2;  

/* The squares of the vectors, and their mutual scalar product */
double  Xsqr, Ysqr, XY;
/* the determinant of the linear system */
double  delta;

/*
 *  The "inverse" mixing matrix: which relates the "unmixing" couplings
 *  "a" and "b" to the old (mixing) ones ("eta2" and "xi" correspondingly)
 *  This corresponds to R^{-1}
 */
double  R11, R12;
double  R21, R22;


/**
 ** Here we are not certain yet how the number of flavors affects
 **/

/*
 * The number of flavours 
 */
double N_f = 1.0;

/*
 * The first coefficient of the QED beta function (without pi^2)
 */
double beta_0 = 1.0 / 12.0;

/*
 * QED charge we accept at scale of 1 GeV
 */
double e_0_2 = 4.0 * M_PI / 137.0;

/* 
 * The ratio of the coupling constants;
 */
double alpha = 0.0;
double beta = 0.0;
/* assisting variable */
double quot = 0.0;

/*
 * The power of the enhancement/suppression
 */
double power_a, power_b;

/****
 ** Here we store scaled vectors and all related stuff
 ****/
double chi_eta2_1, chi_eta2_2;   /* transformed unit-length eta2 vector */
double chi_xi_1,   chi_xi_2;     /* transformed unit-length xi vector   */

/*
 * Lengths of the above vectors 
 */
double len_chi_eta2 = 0.0;             
double len_chi_xi   = 0.0;

/*
 * Angles of the RG transformation of unit-length basis vectors 
 */
double phi_eta2;
double phi_xi;
double theta;     /* angle in between \chi^\xi and \chi^{\eta_2} */

/*
 * Enhancement of width of the box aligned along the \xi-axis
 */
double width = 0.0;

/*
 * This routine does the actual job: scales the two
 * coupling constants from \mu = 1 GeV to \mu = 10^{19} GeV.
 */
void
scaleup_couplings(double  init_eta2, double init_xi,
		  double *final_eta2, double *final_xi)
{
    double a, b;
    double a1, b1;

    /* multiply \chi by R^{-1} */
    a = R11 * init_eta2  +  R12 * init_xi;
    b = R21 * init_eta2  +  R22 * init_xi;

    /* then enhance them by \alpha and \beta */
    a1 = a * alpha;
    b1 = b * beta;

    /* then rotate them back to \chi_1 */
    *final_eta2 = X1 * a1  +  Y1 * b1;
    *final_xi   = X2 * a1  +  Y2 * b1;
}

int
main(int argc, char *argv[])
{

#if 0


  N_f = 1.0;

  quot0 = 1.0 - N_f * 2.0 * beta_0 * (e_0_2 / (M_PI * M_PI)) *
    19.0 * log (10.0);

  quot = pow(quot0, power);

#endif

  printf("--- preliminary check ---\n");
  printf("pi = %g, pi^2 = %g\n", M_PI, M_PI * M_PI);

  printf("\n");
  printf("--- matrix eignevalues ---\n");
  A  =   50.0/96.0;     /* 96-th   =  25/48  */
  B  = - 10.0/96.0;     /* 96-th   =  -5/48  */
  C  =   2.0 /96.0;     /* 96-th   =   1/48  */
  D  =   16.0/96.0;     /* 96-th   =   1/6   */
  printf(" Matrix:  %g   %g\n", A, B);
  printf("          %g   %g\n", C, D);

  kappa1 = ( (A + D) + sqrt ( (A - D) * (A - D)   +   4.0 * B * C ) ) / 2.0;
  kappa2 = ( (A + D) - sqrt ( (A - D) * (A - D)   +   4.0 * B * C ) ) / 2.0;
  printf("kappa1 = %g, kappa2 = %g\n", kappa1, kappa2);

  printf("\n");
  printf("--- calculating eigenvectors ---\n");
#if 0
  X1 = - B / ( A - kappa1 );
  X2 =   1.0;
#endif
  X1 = 1.0;
  X2 = - C / ( D - kappa1 );
  printf("X = (%g, %g)\n", X1, X2);

  Y1 = - B / ( A - kappa2 );
  Y2 =   1.0;
  printf("Y = (%g, %g)\n", Y1, Y2);

  /* now we check that X and Y are indeed eigen vectors */
  otherX1 = A * X1  +  B * X2;
  otherX2 = C * X1  +  D * X2;
  otherY1 = A * Y1  +  B * Y2;
  otherY2 = C * Y1  +  D * Y2;
  printf("checking the X eigenvector (kappa1):\n");
  printf("  X[1]'/ X[1] = %g, X[2]'/ X[2] = %g\n", 
	 otherX1 / X1, otherX2 / X2);
  printf("checking the Y eigenvector (kappa2):\n");
  printf("  Y[1]'/ Y[1] = %g, Y[2]'/ Y[2] = %g\n", 
	 otherY1 / Y1, otherY2 / Y2);

  /* now these are not used, but let them stay */
  Xsqr = X1 * X1 + X2 * X2;
  Ysqr = Y1 * Y1 + Y2 * Y2;
  XY   = X1 * Y1 + X2 * Y2;
  printf("Xsqr = %g, Ysqr = %g, XY = %g\n", Xsqr, Ysqr, XY);

  delta = X1 * Y2 - X2 * Y1;
  R11 =  Y2 / delta;
  R12 = -Y1 / delta;
  R21 = -X2 / delta;
  R22 =  X1 / delta;
  printf("delta = %g\n", delta);
  printf("The inverse mixing matrix (R^{-1}) is:\n");
  printf("   %g   %g  \n", R11, R12);
  printf("   %g   %g  \n", R21, R22);

  power_a = - kappa1 / (2.0 * beta_0);
  power_b = - kappa2 / (2.0 * beta_0);
  printf("power_a = %g, power_b = %g\n", power_a, power_b);

  quot = 1.0 - 2.0 * beta_0 * (e_0_2 / (M_PI * M_PI)) *
    19.0 * log (10.0);
  printf("quot = %g\n", quot);

  printf("\n");
  printf("--- the enhancement coefficients for the eigen-vectors ---\n");
  alpha = pow(quot, power_a);
  beta  = pow(quot, power_b);
  printf("alpha = %g, beta = %g\n", alpha, beta);

  /*
   * Now we have to get the transformation of unit-length basis
   * vectors
   */
  printf("\n");
  printf("--- RG-scaled unit-length(formerly) basis vectors ---\n");
  scaleup_couplings(1.0, 0.0, &chi_eta2_1, &chi_eta2_2);
  scaleup_couplings(0.0, 1.0, &chi_xi_1, &chi_xi_2);
  printf("\\chi^{\\eta_2} = (%g, %g)\n", chi_eta2_1, chi_eta2_2);
  printf("\\chi^\\xi      = (%g, %g)\n", chi_xi_1, chi_xi_2);

  /*
   * Now we compute their lengths
   */
  len_chi_eta2 = sqrt ( chi_eta2_1 * chi_eta2_1   +  chi_eta2_2 * chi_eta2_2 );
  len_chi_xi   = sqrt (   chi_xi_1 * chi_xi_1     +    chi_xi_2 * chi_xi_2   );
  printf("length of \\chi^{\\eta_2} = %g\n", len_chi_eta2);
  printf("length of \\chi^\\xi = %g\n", len_chi_xi);

  /*
   * Now we obtain their angles
   */
  printf("\n");
  printf("--- lengths, rotations and width enhancement ---\n");
  phi_eta2 = atan2 ( chi_eta2_2 , chi_eta2_1 );
  phi_xi   = atan2 ( - chi_xi_1 , chi_xi_2   );
  theta    = phi_xi  +  0.5 * M_PI - phi_eta2;
  width    = len_chi_eta2 * sin ( theta );
  printf("\\phi^{\\eta_2} = %g degrees\n", phi_eta2 * 180.0 / M_PI);
  printf("\\phi^{\\xi}    = %g degrees\n", phi_xi * 180.0 / M_PI);
  printf("\\theta        = %g degrees\n",  theta * 180.0 / M_PI);
  printf("width enhancement of a box aligned along \\xi-axis is %g\n",
	 width);
}
