// HARDIK DAVE
//LID DRIVEN CAVITY PROBLEM using Stream-vorticity PDE Formulation
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void main()
{
  int i, j, m, n;
  double dx, dy, Re, beta;
  printf("Input the value of m\n");
  scanf("%d", &m);

  printf("Input the value of n\n");
  scanf("%d", &n);

  printf("Input the value of Re (Reynolds No)\n");
  scanf("%lf", &Re);
  dx = pow((m - 1), -1);
  dy = pow((n - 1), -1);

  beta = dx / dy;

  double u[m][n], v[m][n];
  double psi[m][n], psi_old[m][n], omega[m][n], omega_old[m][n];
  double error_psi, error_omega;
  int iteration = 1;

  // Intialization and Boundary conditions
  for (j = 0; j < n; j++)
  {
    for (i = 0; i < m; i++)
    {
      if (i == 0) // Left wall bc
      {
        u[i][j] = 0;
        v[i][j] = 0;
        psi[i][j] = 0;
      }
      else if (i == (m - 1)) // Right wall bc
      {
        u[i][j] = 0;
        v[i][j] = 0;
        psi[i][j] = 0;
      }
      else if (j == 0) // Bottom wall bc
      {
        u[i][j] = 0;
        v[i][j] = 0;
        psi[i][j] = 0;
      }
      else if (j == (n - 1)) // Top wall bc
      {
        u[i][j] = 1;
        v[i][j] = 0;
        psi[i][j] = 0;
      }
      else // Interior points
      {
        u[i][j] = 0;
        v[i][j] = 0;
        psi[i][j] = 0;
      }
    }
  }
  // Omega intialisation
  for (j = 0; j < n; j++)
  {
    for (i = 0; i < m; i++)
    {
      if (i == 0) // Omega at Left boundary
      {
        omega[i][j] = 2 * (psi[i + 1][j] - psi[i][j]);
        omega[i][j] = -omega[i][j] / pow(dx, 2);
      }

      else if (i == (m - 1)) // Omega at Right wall
      {
        omega[i][j] = 2 * (psi[i - 1][j] - psi[i][j]);
        omega[i][j] = -omega[i][j] / pow(dx, 2);
      }
      else if (j == 0) // Omega at Bottom wall
      {
        omega[i][j] = 2 * (psi[i][j + 1] - psi[i][j]);
        omega[i][j] = -omega[i][j] / pow(dy, 2);
      }
      else if (j == (n - 1)) // Omega at Top wall
      {
        omega[i][j] = 2 * (psi[i][j - 1] - psi[i][j] + (u[i][j] * dy));
        omega[i][j] = -omega[i][j] / pow(dy, 2);
      }
      else // Interior points
      {
        omega[i][j] = 0;
      }
    }
  }
  // Gauss Sidiel Method
  do
  {
    for (j = 0; j < n; j++)
    {
      for (i = 0; i < m; i++) // assign to old values
      {
        psi_old[i][j] = psi[i][j];
        omega_old[i][j] = omega[i][j];
      }
    }
    // Calculation for vorticity "omega"
    for (j = 1; j < (n - 1); j++)
    {
      for (i = 1; i < (m - 1); i++)
      {
        omega[i][j] = ((1.0 / (2.0 * (1.0 + pow(beta, 2)))) * (((1.0 - ((psi[i][j + 1] - psi[i][j - 1]) * ((beta * Re) / 4.0))) * omega[i + 1][j]) + ((1.0 + ((psi[i][j + 1] - psi[i][j - 1]) * ((beta * Re) / 4.0))) * omega[i - 1][j]) + ((1.0 + ((psi[i + 1][j] - psi[i - 1][j]) * (Re / (4.0 * beta)))) * (pow(beta, 2) * omega[i][j + 1])) + ((1.0 - ((psi[i + 1][j] - psi[i - 1][j]) * (Re / (4.0 * beta)))) * (pow(beta, 2) * omega[i][j - 1]))));
      }
    }
    // Calculation for stream function "psi"
    for (j = 1; j < (n - 1); j++)
    {
      for (i = 1; i < (m - 1); i++)
      {
        psi[i][j] = ((1.0 / (2.0 * (1 + pow(beta, 2)))) * (psi[i + 1][j] + psi[i - 1][j] + (pow(beta, 2) * (psi[i][j + 1] + psi[i][j - 1])) + (pow(dx, 2) * omega[i][j])));
      }
    }
    // Update vorticity after final stream function calculation
    for (j = 0; j < n; j++)
    {
      for (i = 0; i < m; i++)
      {
        if (i == 0) // Omega at Left wall
        {
          omega[i][j] = 2 * (psi[i + 1][j] - psi[i][j]);
          omega[i][j] = -omega[i][j] / pow(dx, 2);
        }

        else if (i == (m - 1)) // Omega at Right wall
        {
          omega[i][j] = 2 * (psi[i - 1][j] - psi[i][j]);
          omega[i][j] = -omega[i][j] / pow(dx, 2);
        }
        else if (j == 0) // Omega at Bottom wall
        {
          omega[i][j] = 2 * (psi[i][j + 1] - psi[i][j]);
          omega[i][j] = -omega[i][j] / pow(dy, 2);
        }
        else if (j == (n - 1)) // Omega at Top wall
        {
          omega[i][j] = 2 * (psi[i][j - 1] - psi[i][j] + (u[i][j] * dy));
          omega[i][j] = -omega[i][j] / pow(dy, 2);
        }
      }
    }
    // Error calculation
    error_psi = 0;
    error_omega = 0;

    for (j = 1; j < (n - 1); j++)
    {
      for (i = 1; i < (m - 1); i++)
      {
        error_psi = error_psi + (pow((psi[i][j] - psi_old[i][j]), 2));
        error_omega = error_omega + (pow((omega[i][j] - omega_old[i][j]), 2));
      }
    }
    error_psi = sqrt((error_psi / ((m - 2) * (n - 2))));
    error_omega = sqrt((error_omega / ((m - 2) * (n - 2))));

    printf("iteration=%d\t", iteration);
    printf("error=%.20lf\n", error_psi);
    printf("error=%.20lf\n", error_omega);
    printf("\n\n");

    iteration++;
  } while (error_psi > 1e-6 || error_omega > 1e-6);

  // Update velocity u and v
  for (j = 1; j < (n - 1); j++)
  {
    for (i = 1; i < (m - 1); i++)
    {
      // central difference discretization
      u[i][j] = (psi[i][j + 1] - psi[i][j - 1]) / (2 * dy);
      v[i][j] = -(psi[i + 1][j] - psi[i - 1][j]) / (2 * dx);
    }
  }
  FILE *fp1 ;
  double x,y;
  fp1 = fopen("All DATA.dat", "w");

  fprintf(fp1, "ZONE I=%d, J=%d\n",m,n);
  for (j = 0; j < n; j++)
  {
    y = j*dy;
    for (i = 0; i < m; i++)
    {
      x= i*dx;
      fprintf(fp1, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", x,y,u[i][j],v[i][j], psi[i][j],omega[i][j]);
    }
  }
 FILE *fp2;
	fp2 = fopen("u_centerline.dat", "w");
	for (j = 0; j < n; j++)
		fprintf(fp2, "%lf \t%lf \n", u[50][j], j * dy);

	FILE *fp3;
	fp3 = fopen("v_centerline.dat", "w");
	for (i = 0; i < m; i++)
		fprintf(fp3, "%lf \t%lf \n", v[i][50], i * dx);
    
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
}
