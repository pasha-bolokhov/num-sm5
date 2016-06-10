#include <stdio.h>
#include <math.h>

/*
 * GNUPLOT testing commands:             e^2
 *                                        
 * print (-5.0/48.0)*1.0*log(0.1/0.001)*(0.0917)/(pi*pi)
 * 
 * print (-5.0/48.0)*1.0*log(2.0/0.1) * (0.0917)/(pi*pi)
 *
 * print (-5.0/48.0)*1.0*log(100.0/2.0)*(0.09547100)/(pi*pi)
 */

/*
 * How many intervals to break the scale into,
 * for the solution of the differential equation
 */
#define NUMBER_INTERVALS          10
double number_intervals = (double)NUMBER_INTERVALS;

#define LEAP1_from    0.001    /*   1 MeV */
#define LEAP1_to      0.100    /* 100 MeV */
#define LEAP2_from    0.100    /* 100 MeV */
#define LEAP2_to      2.0      /*   2 GeV */
#define LEAP3_from    2.0      /*   2 GeV */
#define LEAP3_to      100.0    /* 100 GeV */

double e02 = 0.0;         /* to be defined at each leap */
double mu0 = 0.0;         /* to be defined at each leap */

/*
 * Divergency coefficients of diagrams.
 * For convenience, we've dragged all common factors
 * of 1/\pi^2 and e^2 outside the matrix
 */

double cAA =   19.0 / 48.0;
double cGG =   1.0  / 8.0;
double cLL = - 1.0  / 48.0;
double cKK = - 5.0  / 48.0;

double calphaQED = 1.0 / 6.0;    /* what we call \alpha_1^{QED} */

double pi = M_PI;
double pi2 = -1000.0;


/* 
 * number of leptonic flavors 
 */
double num_lept_flav = 1.0;  

/*
 * LEAP 1
 * This is the correspondence of indices to the operators:
 *     0,1 = \eta^e_1   \eta^e_2
 *     2   = \xi
 */
double chi_LEAP1[3];

/*
 * RG Mixing matrix for the LEAP 1:
 * includes \eta^e_1, \eta^e_2, \xi
 */
double MM1[3][3];

double init_MM1(void)
{
  /* first line */
  MM1[0][0] = cAA + cGG;
  MM1[0][1] = 0.0;
  MM1[0][2] = 0.0;

  /* second line */
  MM1[1][0] = 0.0;
  MM1[1][1] = cAA + cGG;
  MM1[1][2] = cKK;

  /* third line */
  MM1[2][0] = 0.0;
  MM1[2][1] = - cLL;
  MM1[2][2] = num_lept_flav * calphaQED;
}


/*
 * LEAP 2
 * This is the correspondence of indices to the operators:
 *     0,1 = \eta^e_1   \eta^e_2
 *     2,3 = \eta^\mu_1 \eta^\mu_2
 *     4   = \xi
 */
double chi_LEAP2[5];

/*
 * RG Mixing matrix for the LEAP 2:
 * includes \eta^e_1, \eta^e_2, \eta^\mu_1, \eta^\mu_2, \xi
 */
double MM2[5][5];

double init_MM2(void)
{
  /* first line  --  \eta^e_1 */
  MM2[0][0] = cAA + cGG;
  MM2[0][1] = 0.0;
  MM2[0][2] = 0.0;
  MM2[0][3] = 0.0;
  MM2[0][4] = 0.0;

  /* second line -- \eta^e_2 */
  MM2[1][0] = 0.0;
  MM2[1][1] = cAA + cGG;
  MM2[1][2] = 0.0;
  MM2[1][3] = 0.0;
  MM2[1][4] = cKK;

  /* third line -- \eta^\mu_1 */
  MM2[2][0] = 0.0;
  MM2[2][1] = 0.0;
  MM2[2][2] = cAA + cGG;
  MM2[2][3] = 0.0;
  MM2[2][4] = 0.0;

  /* fourth line -- \eta^\mu_2 */
  MM2[3][0] = 0.0;
  MM2[3][1] = 0.0;
  MM2[3][2] = 0.0;
  MM2[3][3] = cAA + cGG;
  MM2[3][4] = cKK;

  /* fifth line -- \xi */
  MM2[4][0] = 0.0;
  MM2[4][1] = - cLL;
  MM2[4][2] = 0.0;
  MM2[4][3] = - cLL;
  MM2[4][4] = num_lept_flav * calphaQED;
}

/*
 * LEAP 3
 * This is the correspondence of indices to the operators:
 *     0,1 = \eta^e_1   \eta^e_2
 *     2,3 = \eta^\mu_1 \eta^\mu_2
 *     4,5 = \eta^\tau_1 \eta^\tau_2
 *     6   = \xi
 */
double chi_LEAP3[7];

/*
 * RG Mixing matrix for the LEAP 3:
 * includes \eta^e_1, \eta^e_2, \eta^\mu_1, \eta^\mu_2, 
 *          \eta^\tau_1, \eta^\tau_2,\xi
 */
double MM3[7][7];

double init_MM3(void)
{
  /* first line  --  \eta^e_1 */
  MM3[0][0] = cAA + cGG;
  MM3[0][1] = 0.0;
  MM3[0][2] = 0.0;
  MM3[0][3] = 0.0;
  MM3[0][4] = 0.0;
  MM3[0][5] = 0.0;
  MM3[0][6] = 0.0;


  /* second line -- \eta^e_2 */
  MM3[1][0] = 0.0;
  MM3[1][1] = cAA + cGG;
  MM3[1][2] = 0.0;
  MM3[1][3] = 0.0;
  MM3[1][4] = 0.0;
  MM3[1][5] = 0.0;
  MM3[1][6] = cKK;

  /* third line -- \eta^\mu_1 */
  MM3[2][0] = 0.0;
  MM3[2][1] = 0.0;
  MM3[2][2] = cAA + cGG;
  MM3[2][3] = 0.0;
  MM3[2][4] = 0.0;
  MM3[2][5] = 0.0;
  MM3[2][6] = 0.0;

  /* fourth line -- \eta^\mu_2 */
  MM3[3][0] = 0.0;
  MM3[3][1] = 0.0;
  MM3[3][2] = 0.0;
  MM3[3][3] = cAA + cGG;
  MM3[3][4] = 0.0;
  MM3[3][5] = 0.0;
  MM3[3][6] = cKK;

  /* fifth line -- \eta^\tau_1 */
  MM3[4][0] = 0.0;
  MM3[4][1] = 0.0;
  MM3[4][2] = 0.0;
  MM3[4][3] = 0.0;
  MM3[4][4] = cAA + cGG;
  MM3[4][5] = 0.0;
  MM3[4][6] = 0.0;

  /* sixth line -- \eta^\tau_2 */
  MM3[5][0] = 0.0;
  MM3[5][1] = 0.0;
  MM3[5][2] = 0.0;
  MM3[5][3] = 0.0;
  MM3[5][4] = 0.0;
  MM3[5][5] = cAA + cGG;
  MM3[5][6] = cKK;

  /* seventh line -- \xi */
  MM3[6][0] = 0.0;
  MM3[6][1] = - cLL;
  MM3[6][2] = 0.0;
  MM3[6][3] = - cLL;
  MM3[6][4] = 0.0;
  MM3[6][5] = - cLL;
  MM3[6][6] = num_lept_flav * calphaQED;
}

/*
 * Calculates the value of the electromagnetic coupling
 * at the scale \mu
 */
double
e2_coupling(double mu, double num_lept_flavors)
{

  return e02 / ( 1.0 - 2.0 * num_lept_flavors * e02 * log(mu / mu0) /
		 (12.0 * pi2));
}

int
main(int argc, char *argv[])
{
  int    i, intno;
  double start_mu = 0.0;
  double end_mu   = 0.0;
  double delta_mu = 0.0;

  pi2 = pi * pi;
  printf("----------- preliminary check -----------\n");
  printf("pi = %g\n", pi);

  {
    double tmp;

    e02           = 4.0 * pi * (1./137.); /* define it at the scale \mu0 */
    mu0           = LEAP1_from;
    tmp =  e2_coupling(0.1, 1.0) / (4.0 * pi);
    printf("e2_coupling(100 MeV, 1 flavor) / (4 pi) = 1/%g\n",
	   1.0/tmp);

    tmp =  e2_coupling(100.0, 3.0) / (4.0 * pi); /* actually we get more flavors */
    printf("e2_coupling(100 GeV, 3 flavors) / (4 pi) = 1/%g\n",
	   1.0/tmp);
  }


  /**********
   *  LEAP 1
   **********/

  printf("----------- LEAP 1 ----------\n");
  num_lept_flav = 1.0;
  e02           = 4.0 * pi * (1./137.); /* define it at the scale \mu0 */
  mu0           = LEAP1_from;
  init_MM1();
  printf("setting all \\eta's to zero, \\xi to 1.0\n");
  chi_LEAP1[0] = 0.0;
  chi_LEAP1[1] = 0.0;
  chi_LEAP1[2] = 1.0;
  /* Lazy computation */
  {
    double tmp[3];

    memset(tmp, 0, sizeof(tmp));
    /* direct try - no solving */
    for (i = 0; i < 3; i++)
      {
	int k;

	tmp[i] = chi_LEAP1[i];
	for (k = 0; k < 3; k++)
	  {
	    tmp[i] += MM1[i][k] * chi_LEAP1[k] * log(LEAP1_to/LEAP1_from) *
	      e2_coupling(LEAP1_from, 1.0) / pi2;
	  }
      }

    printf("lazy computation gives at 100 MeV:\n");
    printf("  \\eta^e_1 = %f, \\eta^e_2 = %f, \\xi = %f\n",
	   tmp[0], tmp[1], tmp[2]);
  }

  /* Breaking into intervals */
  start_mu = log10(LEAP1_from);
  end_mu   = log10(LEAP1_to);
  delta_mu = (end_mu - start_mu) / number_intervals;
  printf("delta_log mu = %f\n", delta_mu);
  for (intno = 1; intno <= number_intervals; intno++)
    {
      double tmp[3];
      
      double dintno = (double)intno;
      double mu = pow(10.0, start_mu  + delta_mu * intno);

      double e2 = e2_coupling(mu, num_lept_flav);

      printf("GG scale %f, e^2/4*pi = 1/%g\n", mu, 1.0/(e2/(4.0*pi)));
      memset(tmp, 0, sizeof(tmp));

      for (i = 0; i < 3; i++)
	{
	  int k;

	  tmp[i] = chi_LEAP1[i];
	  for (k = 0; k < 3; k++)
	    {
	      tmp[i] += MM1[i][k] * chi_LEAP1[k] *
		(delta_mu * log(10.0)) * e2 / pi2;
	    }
	}
      memcpy(chi_LEAP1, tmp, sizeof(tmp));
    }
  printf("***************** at %f GeV\n", LEAP1_to);
  printf("Continuous solving in %d steps gives:\n", NUMBER_INTERVALS);
  printf("  \\eta^e_1 = %f, \\eta^e_2 = %f, \\xi = %f\n",
	 chi_LEAP1[0], chi_LEAP1[1], chi_LEAP1[2]);

  /**********
   *  LEAP 2
   **********/

  printf("----------- LEAP 2 ----------\n");
  num_lept_flav = 2.0;
  e02 = e2_coupling(LEAP2_from, 1.0);
  mu0 = LEAP2_from;
  init_MM2();
  /* copy \chi from the previous leap */
  memset(chi_LEAP2, 0, sizeof(chi_LEAP2));
  chi_LEAP2[0] = chi_LEAP1[0];
  chi_LEAP2[1] = chi_LEAP1[1];
  chi_LEAP2[2] = 0.0;
  chi_LEAP2[3] = 0.0;
  chi_LEAP2[4] = chi_LEAP1[2];

  /* Breaking into intervals */
  start_mu = log10(LEAP2_from);
  end_mu   = log10(LEAP2_to);
  delta_mu = (end_mu - start_mu) / number_intervals;
  for (intno = 1; intno <= number_intervals; intno++)
    {
      double tmp[5];
      
      double dintno = (double)intno;
      double mu = pow(10.0, start_mu  + delta_mu * intno);

      double e2 = e2_coupling(mu, num_lept_flav);

      printf("GG scale %f, e^2/4*pi = 1/%g\n", mu, 1.0/(e2/(4.0*pi)));
      memset(tmp, 0, sizeof(tmp));

      for (i = 0; i < 5; i++)
	{
	  int k;

	  tmp[i] = chi_LEAP2[i];
	  for (k = 0; k < 5; k++)
	    {
	      tmp[i] += MM2[i][k] * chi_LEAP2[k] *
		(delta_mu * log(10.0)) * e2 / pi2;
	    }
	}
      memcpy(chi_LEAP2, tmp, sizeof(tmp));
    }
  printf("***************** at %f GeV\n", LEAP2_to);
  printf("Continuous solving in %d steps gives:\n", NUMBER_INTERVALS);
  printf("  \\eta^e_1 = %f, \\eta^e_2 = %f,\n", chi_LEAP2[0], chi_LEAP2[1]);
  printf("  \\eta^\\mu_1 = %f, \\eta^\\mu_2 = %f,\n", chi_LEAP2[2], chi_LEAP2[3]);
  printf("  \\xi = %f\n", chi_LEAP2[4]);


  /**********
   *  LEAP 3
   **********/

  printf("----------- LEAP 3 ----------\n");
  num_lept_flav = 3.0;
  e02 = e2_coupling(LEAP3_from, 2.0);
  mu0 = LEAP3_from;
  init_MM3();
  /* copy \chi from the previous leap */
  memset(chi_LEAP3, 0, sizeof(chi_LEAP3));
  chi_LEAP3[0] = chi_LEAP2[0];
  chi_LEAP3[1] = chi_LEAP2[1];
  chi_LEAP3[2] = chi_LEAP2[2];
  chi_LEAP3[3] = chi_LEAP2[3];
  chi_LEAP3[4] = 0;
  chi_LEAP3[5] = 0;
  chi_LEAP3[6] = chi_LEAP2[4];

  /* Breaking into intervals */
  start_mu = log10(LEAP3_from);
  end_mu   = log10(LEAP3_to);
  delta_mu = (end_mu - start_mu) / number_intervals;
  for (intno = 1; intno <= number_intervals; intno++)
    {
      double tmp[7];
      
      double dintno = (double)intno;
      double mu = pow(10.0, start_mu  + delta_mu * intno);

      double e2 = e2_coupling(mu, num_lept_flav);

      printf("GG scale %f, e^2/4*pi = 1/%g\n", mu, 1.0/(e2/(4.0*pi)));
      memset(tmp, 0, sizeof(tmp));

      for (i = 0; i < 7; i++)
	{
	  int k;

	  tmp[i] = chi_LEAP3[i];
	  for (k = 0; k < 7; k++)
	    {
	      tmp[i] += MM3[i][k] * chi_LEAP3[k] *
		(delta_mu * log(10.0)) * e2 / pi2;
	    }
	}
      memcpy(chi_LEAP3, tmp, sizeof(tmp));
    }
  printf("***************** at %f GeV\n", LEAP3_to);
  printf("Continuous solving in %d steps gives:\n", NUMBER_INTERVALS);
  printf("  \\eta^e_1 = %f, \\eta^e_2 = %f,\n", chi_LEAP3[0], chi_LEAP3[1]);
  printf("  \\eta^\\mu_1 = %f, \\eta^\\mu_2 = %f,\n", chi_LEAP3[2], chi_LEAP3[3]);
  printf("  \\eta^\\tau_1 = %f, \\eta^\\tau_2 = %f,\n", chi_LEAP3[4], chi_LEAP3[5]);
  printf("  \\xi = %f\n", chi_LEAP3[6]);
}
