/* CGS -- Weiss, Algorithm 11
 * Copyright (C) 2006-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: cgs.c,v 2.5 2007/11/25 18:49:02 kichiki Exp $
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


/* Ref: Weiss, Algorithm 11 CGS
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
cgs (int n, const double *b, double *x,
     void (*atimes) (int, const double *, double *, void *),
     void *atimes_param,
     struct iter *it)
{
#ifndef HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /* use Fortran BLAS routines */

  int i_1 = 1;
  double d_m1 = -1.0;
  double d_2 = 2.0;

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H

  int ret = -1;
  double eps2 = it->eps * it->eps;
  int itmax = it->max;

  double *r  = (double *)malloc (sizeof (double) * n);
  double *r0 = (double *)malloc (sizeof (double) * n);
  double *p  = (double *)malloc (sizeof (double) * n);
  double *u  = (double *)malloc (sizeof (double) * n);
  double *ap = (double *)malloc (sizeof (double) * n);
  double *q  = (double *)malloc (sizeof (double) * n);
  double *t  = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (r,  "cgs");
  CHECK_MALLOC (r0, "cgs");
  CHECK_MALLOC (p,  "cgs");
  CHECK_MALLOC (u,  "cgs");
  CHECK_MALLOC (ap, "cgs");
  CHECK_MALLOC (q,  "cgs");
  CHECK_MALLOC (t,  "cgs");


  double r0ap;
  double rho, rho1;
  double delta;
  double beta;

  double res2 = 0.0;

#ifdef HAVE_CBLAS_H
  /**
   * ATLAS version
   */

  double b2 = cblas_ddot (n, b, 1, b, 1); // (b,b)
  eps2 *= b2;

  // initial residue
  atimes (n, x, r, atimes_param); // r = A.x
  cblas_daxpy (n, -1.0, b, 1, r, 1); // r = A.x - b

  cblas_dcopy (n, r, 1, r0, 1); // r0* = r
  cblas_dcopy (n, r, 1, p, 1); // p = r
  cblas_dcopy (n, r, 1, u, 1); // u = r

  rho = cblas_ddot (n, r0, 1, r, 1); // rho = (r0*, r)

  int i;
  for (i = 0; i < itmax; i ++)
    {
      atimes (n, p, ap, atimes_param); // ap = A.p
      r0ap = cblas_ddot (n, r0, 1, ap, 1); // r0ap = (r0*, A.p)
      delta = - rho / r0ap;

      cblas_dcopy (n, u, 1, q, 1); // q = u
      cblas_dscal (n, 2.0, q, 1); // q = 2 u
      cblas_daxpy (n, delta, ap, 1, q, 1); // q = 2 u + delta A.p

      atimes (n, q, t, atimes_param); // t = A.q

      cblas_daxpy (n, delta, t, 1, r, 1); // r = r + delta t
      cblas_daxpy (n, delta, q, 1, x, 1); // x = x + delta q

      res2 = cblas_ddot (n, r, 1, r, 1);
      if (it->debug == 2)
	{
	  fprintf (it->out, "libiter-cgs %d %e\n", i, res2 / b2);
	}
      if (res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}

      rho1 = cblas_ddot (n, r0, 1, r, 1); // rho = (r0*, r)
      beta = rho1 / rho;
      rho = rho1;

      // here t is not used so that this is used for working area.
      double *qu = t;
      cblas_dcopy (n, q, 1, qu, 1); // qu = q
      cblas_daxpy (n, -1.0, u, 1, qu, 1); // qu = q - u
      cblas_dcopy (n, r, 1, u, 1); // u = r
      cblas_daxpy (n, beta, qu, 1, u, 1); // u = r + beta (q - u)

      cblas_daxpy (n, beta, p, 1, qu, 1); // qu = q - u + beta * p
      cblas_dcopy (n, u, 1, p, 1); // p = u
      cblas_daxpy (n, beta, qu, 1, p, 1); // p = u + beta (q - u + b * p)
    }

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /**
   * BLAS version
   */

  double b2 = ddot_ (&n, b, &i_1, b, &i_1); // (b,b)
  eps2 *= b2;

  // initial residue
  atimes (n, x, r, atimes_param); // r = A.x
  daxpy_ (&n, &d_m1, b, &i_1, r, &i_1); // r = A.x - b

  dcopy_ (&n, r, &i_1, r0, &i_1); // r0* = r
  dcopy_ (&n, r, &i_1, p, &i_1); // p = r
  dcopy_ (&n, r, &i_1, u, &i_1); // u = r

  rho = ddot_ (&n, r0, &i_1, r, &i_1); // rho = (r0*, r)

  int i;
  for (i = 0; i < itmax; i ++)
    {
      atimes (n, p, ap, atimes_param); // ap = A.p
      r0ap = ddot_ (&n, r0, &i_1, ap, &i_1); // r0ap = (r0*, A.p)
      delta = - rho / r0ap;

      dcopy_ (&n, u, &i_1, q, &i_1); // q = u
      dscal_ (&n, &d_2, q, &i_1); // q = 2 u
      daxpy_ (&n, &delta, ap, &i_1, q, &i_1); // q = 2 u + delta A.p

      atimes (n, q, t, atimes_param); // t = A.q

      daxpy_ (&n, &delta, t, &i_1, r, &i_1); // r = r + delta t
      daxpy_ (&n, &delta, q, &i_1, x, &i_1); // x = x + delta q

      res2 = ddot_ (&n, r, &i_1, r, &i_1);
      if (it->debug == 2)
	{
	  fprintf (it->out, "libiter-cgs %d %e\n", i, res2 / b2);
	}
      if (res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}

      rho1 = ddot_ (&n, r0, &i_1, r, &i_1); // rho = (r0*, r)
      beta = rho1 / rho;
      rho = rho1;

      // here t is not used so that this is used for working area.
      double *qu = t;
      dcopy_ (&n, q, &i_1, qu, &i_1); // qu = q
      daxpy_ (&n, &d_m1, u, &i_1, qu, &i_1); // qu = q - u
      dcopy_ (&n, r, &i_1, u, &i_1); // u = r
      daxpy_ (&n, &beta, qu, &i_1, u, &i_1); // u = r + beta (q - u)

      daxpy_ (&n, &beta, p, &i_1, qu, &i_1); // qu = q - u + beta * p
      dcopy_ (&n, u, &i_1, p, &i_1); // p = u
      daxpy_ (&n, &beta, qu, &i_1, p, &i_1); // p = u + beta (q - u + b * p)
    }

# else // !HAVE_BLAS_H
  /**
   * local BLAS version
   */

  double b2 = my_ddot (n, b, 1, b, 1); // (b,b)
  eps2 *= b2;

  // initial residue
  atimes (n, x, r, atimes_param); // r = A.x
  my_daxpy (n, -1.0, b, 1, r, 1); // r = A.x - b

  my_dcopy (n, r, 1, r0, 1); // r0* = r
  my_dcopy (n, r, 1, p, 1); // p = r
  my_dcopy (n, r, 1, u, 1); // u = r

  rho = my_ddot (n, r0, 1, r, 1); // rho = (r0*, r)

  int i;
  for (i = 0; i < itmax; i ++)
    {
      atimes (n, p, ap, atimes_param); // ap = A.p
      r0ap = my_ddot (n, r0, 1, ap, 1); // r0ap = (r0*, A.p)
      delta = - rho / r0ap;

      my_dcopy (n, u, 1, q, 1); // q = u
      my_dscal (n, 2.0, q, 1); // q = 2 u
      my_daxpy (n, delta, ap, 1, q, 1); // q = 2 u + delta A.p

      atimes (n, q, t, atimes_param); // t = A.q

      my_daxpy (n, delta, t, 1, r, 1); // r = r + delta t
      my_daxpy (n, delta, q, 1, x, 1); // x = x + delta q

      res2 = my_ddot (n, r, 1, r, 1);
      if (it->debug == 2)
	{
	  fprintf (it->out, "libiter-cgs %d %e\n", i, res2 / b2);
	}
      if (res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}

      rho1 = my_ddot (n, r0, 1, r, 1); // rho = (r0*, r)
      beta = rho1 / rho;
      rho = rho1;

      // here t is not used so that this is used for working area.
      double *qu = t;
      my_dcopy (n, q, 1, qu, 1); // qu = q
      my_daxpy (n, -1.0, u, 1, qu, 1); // qu = q - u
      my_dcopy (n, r, 1, u, 1); // u = r
      my_daxpy (n, beta, qu, 1, u, 1); // u = r + beta (q - u)

      my_daxpy (n, beta, p, 1, qu, 1); // qu = q - u + beta * p
      my_dcopy (n, u, 1, p, 1); // p = u
      my_daxpy (n, beta, qu, 1, p, 1); // p = u + beta (q - u + b * p)
    }

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H

  free (r);
  free (r0);
  free (p);
  free (u);
  free (ap);
  free (q);
  free (t);

  if (it->debug == 1)
    {
      fprintf (it->out, "libiter-cgs it= %d res^2= %e\n", i, res2);
    }

  it->niter = i;
  it->res2  = res2 / b2;
  return (ret);
}
