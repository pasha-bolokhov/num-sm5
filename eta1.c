
#include <stdio.h>
#include <math.h>

int
main()
{
  double N_f = 3.0;

  double beta_0 = 1.0 / 12.0;

  double e_0_2 = 4.0 * M_PI / 137.0;

  double quot0 = 0.0;
  double quot;

  double power = -25.0 / (8.0 * N_f);

  double mu0 = 1.0; /* GeV */
  double mu_new = 0.0;
  double e_new2;

  quot0 = 1.0 - N_f * 2.0 * beta_0 * (e_0_2 / (M_PI * M_PI)) *
    19.0 * log (10.0);

  quot = pow(quot0, power);

  printf("pi = %g, pi^2 = %g\n", M_PI, M_PI * M_PI);
  printf(" for N_f = %g, quotient (scaling) of the eta_1 coupling constant is"
	 " %g\n", N_f, quot);
  printf(" The power is %g, the quot0 is %g\n", power, quot0);


  mu_new = mu0 * exp ( - (M_PI * M_PI) / ( 2.0 * N_f * beta_0 * e_0_2));
  printf(" new scale is %g\n", mu_new);

  e_new2 = e_0_2 / (1.0 - 2.0 * N_f * beta_0 * (e_0_2 / (M_PI * M_PI)) *
		   log (mu_new/mu0));
  printf(" e^2(mu_new) = %g, e_0^2 = %g\n", e_new2, e_0_2);

}
