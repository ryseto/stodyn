/* Classical CG method -- Weiss' Algorithm 2
 * Copyright (C) 2006-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: cg.c,v 2.5 2007/11/25 18:48:24 kichiki Exp $
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


/* Classical CG method -- Weiss' Algorithm 2
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
cg (int n, const double *b, double *x,
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

  double *p  = (double *)malloc (sizeof (double) * n);
  double *r  = (double *)malloc (sizeof (double) * n);
  double *ap = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (p,  "cg");
  CHECK_MALLOC (r,  "cg");
  CHECK_MALLOC (ap, "cg");

  double r2;
  double res2 = 0.0;
  double pap;

  double gamma;
  double beta;

  int i;

#ifdef HAVE_CBLAS_H
  /**
   * ATLAS version
   */

  double b2 = cblas_ddot (n, b, 1, b, 1); // (b,b)
  eps2 *= b2;

  atimes (n, x, r, atimes_param); // r = A.x
  cblas_daxpy (n, -1.0, b, 1, r, 1); // r = A.x - b
  
  cblas_dcopy (n, r, 1, p, 1); // p = r

  for (i = 0; i < itmax; i ++)
    {
      r2 = cblas_ddot (n, r, 1, r, 1); // r2 = (r, r)
      
      atimes (n, p, ap, atimes_param); // ap = A.p
      pap = cblas_ddot (n, p, 1, ap, 1); // pap = (p, A.p)

      gamma = - r2 / pap; // gamma = - (r, r) / (p, A.p)
      
      cblas_daxpy (n, gamma, p, 1, x, 1); // x += gamma p
      cblas_daxpy (n, gamma, ap, 1, r, 1); // r += gamma Ap

      // new norm of r
      res2 = cblas_ddot (n, r, 1, r, 1); // (r, r)
      if (it->debug == 2)
	{
	  fprintf (it->out, "libiter-cg %d %e\n", i, res2 / b2);
	}
      if (res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}

      beta = res2 / r2; // beta = (r, r) / (r0, r0)
      
      cblas_dscal (n, beta, p, 1); // p *= beta
      cblas_daxpy (n, 1.0, r, 1, p, 1); // p = r + beta p

      r2 = res2;
    }

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /**
   * BLAS version
   */

  double b2 = ddot_ (&n, b, &i_1, b, &i_1); // (b,b)
  eps2 *= b2;

  atimes (n, x, r, atimes_param); // r = A.x
  daxpy_ (&n, &d_m1, b, &i_1, r, &i_1); // r = A.x - b
  
  dcopy_ (&n, r, &i_1, p, &i_1); // p = r

  for (i = 0; i < itmax; i ++)
    {
      r2 = ddot_ (&n, r, &i_1, r, &i_1); // r2 = (r, r)
      
      atimes (n, p, ap, atimes_param); // ap = A.p
      pap = ddot_ (&n, p, &i_1, ap, &i_1); // pap = (p, A.p)

      gamma = - r2 / pap; // gamma = - (r, r) / (p, A.p)
      
      daxpy_ (&n, &gamma, p, &i_1, x, &i_1); // x += gamma p
      daxpy_ (&n, &gamma, ap, &i_1, r, &i_1); // r += gamma Ap

      // new norm of r
      res2 = ddot_ (&n, r, &i_1, r, &i_1); // (r, r)
      if (it->debug == 2)
	{
	  fprintf (it->out, "libiter-cg %d %e\n", i, res2 / b2);
	}
      if (res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}

      beta = res2 / r2; // beta = (r, r) / (r0, r0)
      
      dscal_ (&n, &beta, p, &i_1); // p *= beta
      daxpy_ (&n, &d_1, r, &i_1, p, &i_1); // p = r + beta p

      r2 = res2;
    }

# else // !HAVE_BLAS_H
  /**
   * local BLAS version
   */

  double b2 = my_ddot (n, b, 1, b, 1); // (b,b)
  eps2 *= b2;

  atimes (n, x, r, atimes_param); // r = A.x
  my_daxpy (n, -1.0, b, 1, r, 1); // r = A.x - b
  
  my_dcopy (n, r, 1, p, 1); // p = r

  for (i = 0; i < itmax; i ++)
    {
      r2 = my_ddot (n, r, 1, r, 1); // r2 = (r, r)
      
      atimes (n, p, ap, atimes_param); // ap = A.p
      pap = my_ddot (n, p, 1, ap, 1); // pap = (p, A.p)

      gamma = - r2 / pap; // gamma = - (r, r) / (p, A.p)
      
      my_daxpy (n, gamma, p, 1, x, 1); // x += gamma p
      my_daxpy (n, gamma, ap, 1, r, 1); // r += gamma Ap

      // new norm of r
      res2 = my_ddot (n, r, 1, r, 1); // (r, r)
      if (it->debug == 2)
	{
	  fprintf (it->out, "libiter-cg %d %e\n", i, res2 / b2);
	}
      if (res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}

      beta = res2 / r2; // beta = (r, r) / (r0, r0)
      
      my_dscal (n, beta, p, 1); // p *= beta
      my_daxpy (n, 1.0, r, 1, p, 1); // p = r + beta p

      r2 = res2;
    }

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H

  free (p);
  free (r);
  free (ap);

  if (it->debug == 1)
    {
      fprintf (it->out, "libiter-cg it= %d res^2= %e\n", i, res2 / b2);
    }

  it->niter = i;
  it->res2  = res2 / b2;
  return (ret);
}
