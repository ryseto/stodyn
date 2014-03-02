/* BiCGSTAB - Weiss, Algorithm 12
 * Copyright (C) 2006-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bicgstab.c,v 2.5 2007/11/25 18:43:11 kichiki Exp $
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
#include "config.h"

#include <cstdio>
#include <cstdlib>
#include "libiter.h"
#include "memory-check.h" // CHECK_MALLOC


/* CBLAS is a bit unstable for this scheme,
 * so that here CBLAS is disabled.
 */
#ifdef HAVE_CBLAS_H
#undef HAVE_CBLAS_H
#endif /* HAVE_CBLAS_H */


#ifdef HAVE_CBLAS_H
/* use ATLAS' CBLAS routines */
#ifdef OSX
#include <Accelerate/Accelerate.h>
#else
extern "C"{
#include <atlas/cblas.h>
}
#endif

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
/* use Fortran BLAS routines */

double
ddot_(int* N, 
      double* X, int* incX, 
      double* Y, int* incY);
double
dnrm2_(int* N, 
       double* X, int* incX);
int
dcopy_(int* N,
       double* X, int* incX,
       double* Y, int* incY);
int
daxpy_(int* N,
       double* alpha,
       double* X, int* incX,
       double* Y, int* incY);
int
dscal_(int* N,
       double* alpha,
       double* X, int* incX);

# else // !HAVE_BLAS_H
/* use local BLAS routines */

#include "myblas.h"

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H


/* Ref: Weiss, Algorithm 12 BiCGSTAB
 * INPUT
 *   n : dimension of the problem
 *   b [n] : r-h-s vector
 *   atimes (int n, static double *x, double *b, void *param) :
 *        calc matrix-vector product A.x = b.
 *   atimes_param : parameters for atimes().
 *   it : struct iter. following entries are used
 *        it->max = kend : max of iteration
 *        it->eps = eps  : criteria for |r^2|/|b^2|
 * OUTPUT
 *   returned value : 0 == success, otherwise (-1) == failed
 *   x [n] : solution
 *   it->niter : # of iteration
 *   it->res2  : |r^2| / |b^2|
 */
int
bicgstab (int n, const double *b, double *x,
	  void (*atimes) (int, const double *, double *, void *),
	  void *atimes_param,
	  struct iter *it)
{
#ifndef HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /* use Fortran BLAS routines */

  int i_1 = 1;
  double d_1 = 1.0;
  double d_m1 = -1.0;

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H

  int ret = -1;
  double eps2 = it->eps * it->eps;
  int itmax = it->max;

  double *r  = (double *)malloc (sizeof (double) * n);
  double *rs = (double *)malloc (sizeof (double) * n);
  double *p  = (double *)malloc (sizeof (double) * n);
  double *ap = (double *)malloc (sizeof (double) * n);
  double *s  = (double *)malloc (sizeof (double) * n);
  double *t  = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (r,  "bicgstab");
  CHECK_MALLOC (rs, "bicgstab");
  CHECK_MALLOC (p,  "bicgstab");
  CHECK_MALLOC (ap, "bicgstab");
  CHECK_MALLOC (s,  "bicgstab");
  CHECK_MALLOC (t,  "bicgstab");

  double rsap; // (r*, A.p)
  double st;
  double t2;

  double rho, rho1;
  double delta;
  double gamma;
  double beta;

  double res2 = 0.0;

#ifdef HAVE_CBLAS_H
  /**
   * ATLAS version
   */

  double b2 = cblas_ddot (n, b, 1, b, 1); // (b,b)
  eps2 *= b2;

  atimes (n, x, r, atimes_param);    // r = A.x ...
  cblas_daxpy (n, -1.0, b, 1, r, 1); //         - b

  cblas_dcopy (n, r, 1, rs, 1); // r* = r
  cblas_dcopy (n, r, 1, p, 1);  // p  = r

  rho = cblas_ddot (n, rs, 1, r, 1); // rho = (r*, r)

  int i;
  for (i = 0; i < itmax; i ++)
    {
      atimes (n, p, ap, atimes_param); // ap = A.p
      rsap = cblas_ddot (n, rs, 1, ap, 1); // rsap = (r*, A.p)
      delta = - rho / rsap;

      cblas_dcopy (n, r, 1, s, 1);         // s = r ...
      cblas_daxpy (n, delta, ap, 1, s, 1); //   + delta A.p
      atimes (n, s, t, atimes_param); // t = A.s

      st = cblas_ddot (n, s, 1, t, 1); // st = (s, t)
      t2 = cblas_ddot (n, t, 1, t, 1); // t2 = (t, t)
      gamma = - st / t2;

      cblas_dcopy (n, s, 1, r, 1);        // r = s ...
      cblas_daxpy (n, gamma, t, 1, r, 1); //   + gamma t

      cblas_daxpy (n, delta, p, 1, x, 1); // x = x + delta p...
      cblas_daxpy (n, gamma, s, 1, x, 1); //       + gamma s

      res2 = cblas_ddot (n, r, 1, r, 1);
      if (it->debug == 2)
	{
	  fprintf (it->out,
		   "libiter-bicgstab(cblas) %d %e\n",
		   i, res2 / b2);
	}
      if (res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}

      rho1 = cblas_ddot (n, rs, 1, r, 1); // rho = (r*, r)
      beta = rho1 / rho * delta / gamma;
      rho = rho1;

      cblas_daxpy (n, gamma, ap, 1, p, 1); // p = p + gamma A.p
      cblas_dscal (n, beta, p, 1);         // p = beta (p + gamma A.p)
      cblas_daxpy (n, 1.0, r, 1, p, 1); // p = r + beta(p + gamma A.p)
    }

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /**
   * BLAS version
   */

  double b2 = ddot_ (&n, b, &i_1, b, &i_1); // (b,b)
  eps2 *= b2;

  atimes (n, x, r, atimes_param);       // r = A.x ...
  daxpy_ (&n, &d_m1, b, &i_1, r, &i_1); //         - b

  dcopy_ (&n, r, &i_1, rs, &i_1); // r* = r
  dcopy_ (&n, r, &i_1, p, &i_1);  // p  = r

  rho = ddot_ (&n, rs, &i_1, r, &i_1); // rho = (r*, r)

  int i;
  for (i = 0; i < itmax; i ++)
    {
      atimes (n, p, ap, atimes_param); // ap = A.p
      rsap = ddot_ (&n, rs, &i_1, ap, &i_1); // rsap = (r*, A.p)
      delta = - rho / rsap;

      dcopy_ (&n, r, &i_1, s, &i_1);          // s = r ...
      daxpy_ (&n, &delta, ap, &i_1, s, &i_1); //   + delta A.p
      atimes (n, s, t, atimes_param); // t = A.s

      st = ddot_ (&n, s, &i_1, t, &i_1); // st = (s, t)
      t2 = ddot_ (&n, t, &i_1, t, &i_1); // t2 = (t, t)
      gamma = - st / t2;

      dcopy_ (&n, s, &i_1, r, &i_1);         // r = s ...
      daxpy_ (&n, &gamma, t, &i_1, r, &i_1); //   + gamma t

      daxpy_ (&n, &delta, p, &i_1, x, &i_1); // x = x + delta p...
      daxpy_ (&n, &gamma, s, &i_1, x, &i_1); //       + gamma s

      res2 = ddot_ (&n, r, &i_1, r, &i_1);
      if (it->debug == 2)
	{
	  fprintf (it->out,
		   "libiter-bicgstab(blas) %d %e\n",
		   i, res2 / b2);
	}
      if (res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}
      if (res2 > 1.0e20)
	{
	  // already too big residual
	  break;
	}

      rho1 = ddot_ (&n, rs, &i_1, r, &i_1); // rho = (r*, r)
      beta = rho1 / rho * delta / gamma;
      rho = rho1;

      daxpy_ (&n, &gamma, ap, &i_1, p, &i_1); // p = p + gamma A.p
      dscal_ (&n, &beta, p, &i_1);            // p = beta (p + gamma A.p)
      daxpy_ (&n, &d_1, r, &i_1, p, &i_1); // p = r + beta(p + gamma A.p)
    }

# else // !HAVE_BLAS_H
  /**
   * local BLAS version
   */

  double b2 = my_ddot (n, b, 1, b, 1); // (b,b)
  eps2 *= b2;

  atimes (n, x, r, atimes_param);    // r = A.x ...
  my_daxpy (n, -1.0, b, 1, r, 1); //         - b

  my_dcopy (n, r, 1, rs, 1); // r* = r
  my_dcopy (n, r, 1, p, 1);  // p = r

  rho = my_ddot (n, rs, 1, r, 1); // rho = (r*, r)

  int i;
  for (i = 0; i < itmax; i ++)
    {
      atimes (n, p, ap, atimes_param); // ap = A.p
      rsap = my_ddot (n, rs, 1, ap, 1); // rsap = (r*, A.p)
      delta = - rho / rsap;

      my_dcopy (n, r, 1, s, 1);         // s = r ...
      my_daxpy (n, delta, ap, 1, s, 1); //   + delta A.p
      atimes (n, s, t, atimes_param); // t = A.s

      st = my_ddot (n, s, 1, t, 1); // st = (s, t)
      t2 = my_ddot (n, t, 1, t, 1); // t2 = (t, t)
      gamma = - st / t2;

      my_dcopy (n, s, 1, r, 1);        // r = s ...
      my_daxpy (n, gamma, t, 1, r, 1); //   + gamma t

      my_daxpy (n, delta, p, 1, x, 1); // x = x + delta p...
      my_daxpy (n, gamma, s, 1, x, 1); //       + gamma s

      res2 = my_ddot (n, r, 1, r, 1);
      if (it->debug == 2)
	{
	  fprintf (it->out,
		   "libiter-bicgstab(myblas) %d %e\n",
		   i, res2 / b2);
	}
      if (res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}

      rho1 = my_ddot (n, rs, 1, r, 1); // rho = (r*, r)
      beta = rho1 / rho * delta / gamma;
      rho = rho1;

      my_daxpy (n, gamma, ap, 1, p, 1); // p = p + gamma A.p
      my_dscal (n, beta, p, 1);         // p = beta (p + gamma A.p)
      my_daxpy (n, 1.0, r, 1, p, 1); // p = r + beta(p + gamma A.p)
    }

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H

  free (r);
  free (rs);
  free (p);
  free (ap);
  free (s);
  free (t);

  if (it->debug == 1)
    {
      fprintf (it->out, "libiter-bicgstab %d %e\n", i, res2 / b2);
    }

  it->niter = i;
  it->res2  = res2 / b2;
  return (ret);
}
