/* CG method
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: cg_.c,v 1.2 2007/12/01 18:20:20 kichiki Exp $
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



/* conjugate gradient method, another implementation
 */
int
cg_ (int n, const double *b, double *x,
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

  double *ap = (double *)malloc (sizeof (double) * n);
  double *r  = (double *)malloc (sizeof (double) * n);
  double *p  = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (ap, "cg_");
  CHECK_MALLOC (r,  "cg_");
  CHECK_MALLOC (p,  "cg_");

#ifdef HAVE_CBLAS_H
  /**
   * ATLAS version
   */

  double b2 = cblas_ddot (n, b, 1, b, 1); // b2 = (b, b)
  double eps2 = it->eps * it->eps * b2;
  // eps2 is compared with r2 = (r, r)

  int i = 0;

  // set r
  atimes (n, x, r, atimes_param); // r = A.x
  cblas_daxpy (n, -1.0, b, 1, r, 1); // r = A.x - b
  cblas_dscal (n, -1.0, r, 1); // r = b - A.x

  double res0 = cblas_ddot (n, r, 1, r, 1); // r2 = (r, r)
  double res2 = res0;
  if (it->debug == 2)
    {
      fprintf (it->out, "cg_ %d %e\n", i, res2 / b2);
    }
  if (res2 > eps2)
    {
      // set p
      cblas_dcopy (n, r, 1, p, 1); // p = r

      for (i = 0; i < it->max; i ++)
	{
	  atimes (n, p, ap, atimes_param); // ap = A.p
	  double pap = cblas_ddot (n, p, 1, ap, 1); // pap = (p, A.p)

	  double alpha = res0 / pap;

	  // update x and r
	  cblas_daxpy (n,  alpha, p,  1, x, 1); // x += gamma p
	  cblas_daxpy (n, -alpha, ap, 1, r, 1); // r += gamma Ap

	  res2 = cblas_ddot (n, r, 1, r, 1); // r2 = (r, r) for new r
	  if (it->debug == 2)
	    {
	      fprintf (it->out, "cg_ %d %e\n", i, res2 / b2);
	    }
	  if (res2 <= eps2)
	    {
	      ret = 0; // success
	      break;
	    }

	  double beta = res2 / res0; // beta = (r',r') / (r,r)

	  res0 = res2;

	  // set p
	  cblas_dscal (n, beta, p, 1); // p = beta p
	  cblas_daxpy (n, 1.0, r, 1, p, 1); // p = r + beta p
	}
    }

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /**
   * BLAS version
   */

  double b2 = ddot_ (&n, b, &i_1, b, &i_1); // b2 = (b, b)
  double eps2 = it->eps * it->eps * b2;
  // eps2 is compared with r2 = (r, r)

  int i = 0;

  // set r
  atimes (n, x, r, atimes_param); // r = A.x
  daxpy_ (&n, &d_m1, b, &i_1, r, &i_1); // r = A.x - b
  dscal_ (&n, &d_m1, r, &i_1); // r = b - A.x

  double res0 = ddot_ (&n, r, &i_1, r, &i_1); // r2 = (r, r)
  double res2 = res0;
  if (it->debug == 2)
    {
      fprintf (it->out, "cg_ %d %e\n", i, res2 / b2);
    }
  if (res2 > eps2)
    {
      // set p
      dcopy_ (&n, r, &i_1, p, &i_1); // p = r

      for (i = 0; i < it->max; i ++)
	{
	  atimes (n, p, ap, atimes_param); // ap = A.p
	  double pap = ddot_ (&n, p, &i_1, ap, &i_1); // pap = (p, A.p)

	  double alpha = res0 / pap;

	  // update x and r
	  double malpha = -alpha;
	  daxpy_ (&n, &alpha,  p,  &i_1, x, &i_1); // x += gamma p
	  daxpy_ (&n, &malpha, ap, &i_1, r, &i_1); // r += gamma Ap

	  res2 = ddot_ (&n, r, &i_1, r, &i_1); // r2 = (r, r) for new r
	  if (it->debug == 2)
	    {
	      fprintf (it->out, "cg_ %d %e\n", i, res2 / b2);
	    }
	  if (res2 <= eps2)
	    {
	      ret = 0; // success
	      break;
	    }

	  double beta = res2 / res0; // beta = (r',r') / (r,r)

	  res0 = res2;

	  // set p
	  dscal_ (&n, &beta, p, &i_1); // p = beta p
	  daxpy_ (&n, &d_1,  r, &i_1, p, &i_1); // p = r + beta p
	}
    }

# else // !HAVE_BLAS_H
  /**
   * local BLAS version
   */

  double b2 = my_ddot (n, b, 1, b, 1); // b2 = (b, b)
  double eps2 = it->eps * it->eps * b2;
  // eps2 is compared with r2 = (r, r)

  int i = 0;

  // set r
  atimes (n, x, r, atimes_param); // r = A.x
  my_daxpy (n, -1.0, b, 1, r, 1); // r = A.x - b
  my_dscal (n, -1.0, r, 1); // r = b - A.x

  double res0 = my_ddot (n, r, 1, r, 1); // r2 = (r, r)
  double res2 = res0;
  if (it->debug == 2)
    {
      fprintf (it->out, "cg_ %d %e\n", i, res2 / b2);
    }
  if (res2 > eps2)
    {
      // set p
      my_dcopy (n, r, 1, p, 1); // p = r

      for (i = 0; i < it->max; i ++)
	{
	  atimes (n, p, ap, atimes_param); // ap = A.p
	  double pap = my_ddot (n, p, 1, ap, 1); // pap = (p, A.p)

	  double alpha = res0 / pap;

	  // update x and r
	  my_daxpy (n,  alpha, p,  1, x, 1); // x += gamma p
	  my_daxpy (n, -alpha, ap, 1, r, 1); // r += gamma Ap

	  res2 = my_ddot (n, r, 1, r, 1); // r2 = (r, r) for new r
	  if (it->debug == 2)
	    {
	      fprintf (it->out, "cg_ %d %e\n", i, res2 / b2);
	    }
	  if (res2 <= eps2)
	    {
	      ret = 0; // success
	      break;
	    }

	  double beta = res2 / res0; // beta = (r',r') / (r,r)

	  res0 = res2;

	  // set p
	  my_dscal (n, beta, p, 1); // p = beta p
	  my_daxpy (n, 1.0, r, 1, p, 1); // p = r + beta p
	}
    }

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H

  free (ap);
  free (r);
  free (p);

  if (it->debug == 1)
    {
      fprintf (it->out, "cg_ %d %e\n", i, res2 / b2);
    }

  it->niter = i;
  it->res2  = res2 / b2;
  return (ret);
}

/* conjugate gradient method with preconditioner
 */
int
cg_pc (int n, const double *b, double *x,
       void (*atimes) (int, const double *, double *, void *),
       void *atimes_param,
       void (*inv) (int, const double *, double *, void *),
       void *inv_param,
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

  double *ap = (double *)malloc (sizeof (double) * n);
  double *r  = (double *)malloc (sizeof (double) * n);
  double *p  = (double *)malloc (sizeof (double) * n);
  double *Pr = (double *)malloc (sizeof (double) * n);
  CHECK_MALLOC (ap, "cg_pc");
  CHECK_MALLOC (r,  "cg_pc");
  CHECK_MALLOC (p,  "cg_pc");
  CHECK_MALLOC (Pr, "cg_pc");

#ifdef HAVE_CBLAS_H
  /**
   * ATLAS version
   */

  double b2 = cblas_ddot (n, b, 1, b, 1); // b2 = (b, b)
  double eps2 = it->eps * it->eps * b2;
  // eps2 is compared with r2 = (r, r)

  int i = 0;

  // set r
  atimes (n, x, r, atimes_param); // r = A.x
  cblas_daxpy (n, -1.0, b, 1, r, 1); // r = A.x - b
  cblas_dscal (n, -1.0, r, 1); // r = b - A.x

  double res2 = cblas_ddot (n, r, 1, r, 1); // r2 = (r, r)
  if (it->debug == 2)
    {
      fprintf (it->out, "cg_pc %d %e\n", i, res2 / b2);
    }
  if (res2 > eps2)
    {
      // set p
      inv (n, r, Pr, inv_param);
      cblas_dcopy (n, Pr, 1, p, 1); // p = K^{-1}.r

      double Prr = cblas_ddot (n, Pr, 1, r, 1); // Prr = (Pr, r)

      for (i = 0; i < it->max; i ++)
	{
	  atimes (n, p, ap, atimes_param); // ap = A.p
	  double pap = cblas_ddot (n, p, 1, ap, 1); // pap = (p, A.p)

	  double alpha = Prr / pap;

	  // update x and r
	  cblas_daxpy (n,  alpha, p,  1, x, 1); // x += gamma p
	  cblas_daxpy (n, -alpha, ap, 1, r, 1); // r += gamma Ap

	  res2 = cblas_ddot (n, r, 1, r, 1); // r2 = (r, r) for new r
	  if (it->debug == 2)
	    {
	      fprintf (it->out, "cg_pc %d %e\n", i, res2 / b2);
	    }
	  if (res2 <= eps2)
	    {
	      ret = 0; // success
	      break;
	    }

	  // update Pr
	  inv (n, r, Pr, inv_param);
	  double Pres2 = cblas_ddot (n, Pr, 1, r, 1); // Pres2 = (Pr', r')
	  double beta = Pres2 / Prr; // beta = (Pr',r') / (Pr,r)

	  Prr = Pres2;

	  // set p
	  cblas_dscal (n, beta, p, 1); // p = beta p
	  cblas_daxpy (n, 1.0, Pr, 1, p, 1); // p = Pr + beta p
	}
    }

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /**
   * BLAS version
   */

  double b2 = ddot_ (&n, b, &i_1, b, &i_1); // b2 = (b, b)
  double eps2 = it->eps * it->eps * b2;
  // eps2 is compared with r2 = (r, r)

  int i = 0;

  // set r
  atimes (n, x, r, atimes_param); // r = A.x
  daxpy_ (&n, &d_m1, b, &i_1, r, &i_1); // r = A.x - b
  dscal_ (&n, &d_m1, r, &i_1); // r = b - A.x

  double res2 = ddot_ (&n, r, &i_1, r, &i_1); // r2 = (r, r)
  if (it->debug == 2)
    {
      fprintf (it->out, "cg_pc %d %e\n", i, res2 / b2);
    }
  if (res2 > eps2)
    {
      // set p
      inv (n, r, Pr, inv_param);
      dcopy_ (&n, Pr, &i_1, p, &i_1); // p = K^{-1}.r

      double Prr = ddot_ (&n, Pr, &i_1, r, &i_1); // Prr = (Pr, r)

      for (i = 0; i < it->max; i ++)
	{
	  atimes (n, p, ap, atimes_param); // ap = A.p
	  double pap = ddot_ (&n, p, &i_1, ap, &i_1); // pap = (p, A.p)

	  double alpha = Prr / pap;

	  // update x and r
	  double malpha = -alpha;
	  daxpy_ (&n, &alpha,  p,  &i_1, x, &i_1); // x += gamma p
	  daxpy_ (&n, &malpha, ap, &i_1, r, &i_1); // r += gamma Ap

	  res2 = ddot_ (&n, r, &i_1, r, &i_1); // r2 = (r, r) for new r
	  if (it->debug == 2)
	    {
	      fprintf (it->out, "cg_pc %d %e\n", i, res2 / b2);
	    }
	  if (res2 <= eps2)
	    {
	      ret = 0; // success
	      break;
	    }

	  // update Pr
	  inv (n, r, Pr, inv_param);
	  double Pres2 = ddot_ (&n, Pr, &i_1, r, &i_1); // Pres2 = (Pr', r')
	  double beta = Pres2 / Prr; // beta = (Pr',r') / (Pr,r)

	  Prr = Pres2;

	  // set p
	  dscal_ (&n, &beta, p, &i_1); // p = beta p
	  daxpy_ (&n, &d_1,  Pr, &i_1, p, &i_1); // p = Pr + beta p
	}
    }

# else // !HAVE_BLAS_H
  /**
   * local BLAS version
   */

  double b2 = my_ddot (n, b, 1, b, 1); // b2 = (b, b)
  double eps2 = it->eps * it->eps * b2;
  // eps2 is compared with r2 = (r, r)

  int i = 0;

  // set r
  atimes (n, x, r, atimes_param); // r = A.x
  my_daxpy (n, -1.0, b, 1, r, 1); // r = A.x - b
  my_dscal (n, -1.0, r, 1); // r = b - A.x

  double res2 = my_ddot (n, r, 1, r, 1); // r2 = (r, r)
  if (it->debug == 2)
    {
      fprintf (it->out, "cg_pc %d %e\n", i, res2 / b2);
    }
  if (res2 > eps2)
    {
      // set p
      inv (n, r, Pr, inv_param);
      my_dcopy (n, Pr, 1, p, 1); // p = K^{-1}.r

      double Prr = my_ddot (n, Pr, 1, r, 1); // Prr = (Pr, r)

      for (i = 0; i < it->max; i ++)
	{
	  atimes (n, p, ap, atimes_param); // ap = A.p
	  double pap = my_ddot (n, p, 1, ap, 1); // pap = (p, A.p)

	  double alpha = Prr / pap;

	  // update x and r
	  my_daxpy (n,  alpha, p,  1, x, 1); // x += gamma p
	  my_daxpy (n, -alpha, ap, 1, r, 1); // r += gamma Ap

	  res2 = my_ddot (n, r, 1, r, 1); // r2 = (r, r) for new r
	  if (it->debug == 2)
	    {
	      fprintf (it->out, "cg_pc %d %e\n", i, res2 / b2);
	    }
	  if (res2 <= eps2)
	    {
	      ret = 0; // success
	      break;
	    }

	  // update Pr
	  inv (n, r, Pr, inv_param);
	  double Pres2 = my_ddot (n, Pr, 1, r, 1); // Pres2 = (Pr', r')
	  double beta = Pres2 / Prr; // beta = (Pr',r') / (Pr,r)

	  Prr = Pres2;

	  // set p
	  my_dscal (n, beta, p, 1); // p = beta p
	  my_daxpy (n, 1.0, Pr, 1, p, 1); // p = Pr + beta p
	}
    }

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H

  free (ap);
  free (r);
  free (p);
  free (Pr);

  if (it->debug == 1)
    {
      fprintf (it->out, "cg_pc %d %e\n", i, res2 / b2);
    }

  it->niter = i;
  it->res2  = res2 / b2;
  return (ret);
}

