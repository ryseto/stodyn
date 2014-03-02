/* wrapper for iterative solver routines
 * Copyright (C) 1999-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bi-cgstab.c,v 2.9 2007/11/25 18:42:23 kichiki Exp $
 *
 * solver routines are translated into C by K.I. from fortran code
 * originally written by martin h. gutknecht
 *             ===============rog.f==============
 *             problem given by martin h. gutknecht
 *             numerical method: bi-cgsta        method -- sta ()
 *             numerical method: bi-cgsta2       method -- st2 ()
 *             numerical method: gpbi-cg         method -- gpb ()
 *             ver. 1.0 aug. 08 1995  s. l. zhang
 *    at urus.slzhang.fort.iter.real.gpbcg.gutknecht-p(final.f)
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

#include <cstdio> /* fprintf() */
#include <cmath> /* log10() */
#include <cstdlib> /* malloc(), free() */
#include "memory-check.h" // CHECK_MALLOC

#include "libiter.h" // struct iter

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


#include "bi-cgstab.h"


/* bi-cgstab method
 * INPUT
 *   m : dimension of the problem
 *   b [m] : r-h-s vector
 *   atimes (int m, static double *x, double *b, void *param) :
 *        calc matrix-vector product A.x = b.
 *   atimes_param : parameters for atimes().
 *   it : struct iter. following entries are used
 *        it->max = kend : max of iteration
 *        it->eps = eps  : criteria for |r^2|/|b^2|
 * OUTPUT
 *   returned value : 0 == success, otherwise (-1) == failed
 *   x [m] : solution
 *   it->niter : # of iteration
 *   it->res2  : |r^2| / |b^2|
 */
int
sta (int m, const double *b, double *x,
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
  int kend = it->max;
  double eps2 = it->eps * it->eps;

  double *r0 = (double *)calloc (m, sizeof (double));
  double *p  = (double *)calloc (m, sizeof (double));
  double *q  = (double *)calloc (m, sizeof (double));
  double *t  = (double *)calloc (m, sizeof (double));
  double *r  = (double *)calloc (m, sizeof (double));
  CHECK_MALLOC (r0, "sta");
  CHECK_MALLOC (p,  "sta");
  CHECK_MALLOC (q,  "sta");
  CHECK_MALLOC (t,  "sta");
  CHECK_MALLOC (r,  "sta");

  double res2 = 0.0; // for compiler warning.
  double beta = 0.0;
  double zeta = 0.0;

#ifdef HAVE_CBLAS_H
  /**
   * ATLAS version
   */

  double b2 = cblas_ddot (m, b, 1, b, 1); // (b,b)
  eps2 *= b2;

  // r0 = b - A.x
  atimes (m, x, r0, atimes_param);
  cblas_daxpy (m, -1.0, b, 1, r0, 1);
  cblas_dscal (m, -1.0, r0, 1);
  // r  = r0
  cblas_dcopy (m, r0, 1, r, 1);
  // rho0 = (r0, r0)
  double rho0 = cblas_ddot (m, r0, 1, r0, 1);

  int iter;
  for (iter = 0; iter<= kend; iter++)
    {
      // res2 = (r, r)
      res2 = cblas_ddot (m, r, 1, r, 1);
      if (it->debug == 2)
	{
	  fprintf (it->out, "libiter-sta(atlas) %d %e\n", iter, res2 / b2);
	}
      if(res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}

      // p = r + beta * (p - zeta * q)
      cblas_daxpy (m, -zeta, q, 1, p, 1); // p = p - zeta * q
      cblas_dscal (m, beta, p, 1);        // p = beta * (p - zeta * q)
      cblas_daxpy (m, 1.0, r, 1, p, 1);   // p = r + beta * (p - zeta * q)

      // q   = A.p
      atimes (m, p, q, atimes_param);
      // pap = (r0, A.p)
      double pap = cblas_ddot (m, r0, 1, q, 1);

      double alpha = rho0 / pap;
      // r = r - alpha * q
      cblas_daxpy (m, -alpha, q, 1, r, 1);
      // x = x + alpha * p
      cblas_daxpy (m, +alpha, p, 1, x, 1);

      // t = A.r
      atimes (m, r, t, atimes_param);
      // pap = (t, r)
      pap = cblas_ddot (m, t, 1, r, 1);
      // rur = (t, t)
      double rur = cblas_ddot (m, t, 1, t, 1);

      zeta = pap / rur;
      // x = x + zeta * r
      cblas_daxpy (m, zeta, r, 1, x, 1);

      // r = r - zeta * t
      cblas_daxpy (m, -zeta, t, 1, r, 1);
      // rho1 = (r0, r)
      double rho1 = cblas_ddot (m, r0, 1, r, 1);

      beta = (rho1 / rho0) * (alpha / zeta);
      rho0 = rho1;
    }

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /**
   * BLAS version
   */

  double b2 = ddot_ (&m, b, &i_1, b, &i_1); // (b,b)
  eps2 *= b2;

  // r0 = b - A.x
  atimes (m, x, r0, atimes_param);
  daxpy_ (&m, &d_m1, b, &i_1, r0, &i_1);
  dscal_ (&m, &d_m1, r0, &i_1);
  // r  = r0
  dcopy_ (&m, r0, &i_1, r, &i_1);
  // rho0 = (r0, r0)
  double rho0 = ddot_ (&m, r0, &i_1, r0, &i_1);

  double mzeta = -zeta;
  int iter;
  for (iter = 0; iter<= kend; iter++)
    {
      // res2 = (r, r)
      res2 = ddot_ (&m, r, &i_1, r, &i_1);
      if (it->debug == 2)
	{
	  fprintf (it->out, "libiter-sta(blas) %d %e\n", iter, res2 / b2);
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

      // p = r + beta * (p - zeta * q)
      daxpy_ (&m, &mzeta, q, &i_1, p, &i_1); // p = p - zeta * q
      dscal_ (&m, &beta, p, &i_1);           // p = beta * (p - zeta * q)
      daxpy_ (&m, &d_1, r, &i_1, p, &i_1);   // p = r + beta * (p - zeta * q)

      // q   = A.p
      atimes (m, p, q, atimes_param);
      // pap = (r0, A.p)
      double pap = ddot_ (&m, r0, &i_1, q, &i_1);

      double alpha = rho0 / pap;
      // r = r - alpha * q
      double malpha = -alpha;
      daxpy_ (&m, &malpha, q, &i_1, r, &i_1);
      // x = x + alpha * p
      daxpy_ (&m, &alpha, p, &i_1, x, &i_1);

      // t = A.r
      atimes (m, r, t, atimes_param);
      // pap = (t, r)
      pap = ddot_ (&m, t, &i_1, r, &i_1);
      // rur = (t, t)
      double rur = ddot_ (&m, t, &i_1, t, &i_1);

      zeta = pap / rur;
      // x = x + zeta * r
      daxpy_ (&m, &zeta, r, &i_1, x, &i_1);

      // r = r - zeta * t
      mzeta = -zeta;
      daxpy_ (&m, &mzeta, t, &i_1, r, &i_1);
      // rho1 = (r0, r)
      double rho1 = ddot_ (&m, r0, &i_1, r, &i_1);

      beta = (rho1 / rho0) * (alpha / zeta);
      rho0 = rho1;
    }

# else // !HAVE_BLAS_H
  /**
   * local BLAS version
   */

  double b2 = my_ddot (m, b, 1, b, 1); // (b,b)
  eps2 *= b2;

  // r0 = b - A.x
  atimes (m, x, r0, atimes_param);
  my_daxpy (m, -1.0, b, 1, r0, 1);
  my_dscal (m, -1.0, r0, 1);
  // r  = r0
  my_dcopy (m, r0, 1, r, 1);
  // rho0 = (r0, r0)
  double rho0 = my_ddot (m, r0, 1, r0, 1);

  int iter;
  for (iter = 0; iter<= kend; iter++)
    {
      // res2 = (r, r)
      res2 = my_ddot (m, r, 1, r, 1);
      if (it->debug == 2)
	{
	  fprintf (it->out, "libiter-sta(myblas) %d %e\n", iter, res2 / b2);
	}
      if(res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}

      // p = r + beta * (p - zeta * q)
      my_daxpy (m, -zeta, q, 1, p, 1); // p = p - zeta * q
      my_dscal (m, beta, p, 1);         // p = beta * (p - zeta * q)
      my_daxpy (m, 1.0, r, 1, p, 1);    // p = r + beta * (p - zeta * q)

      // q   = A.p
      atimes (m, p, q, atimes_param);
      // pap = (r0, A.p)
      double pap = my_ddot (m, r0, 1, q, 1);

      double alpha = rho0 / pap;
      // r = r - alpha * q
      my_daxpy (m, -alpha, q, 1, r, 1);
      // x = x + alpha * p
      my_daxpy (m, +alpha, p, 1, x, 1);

      // t = A.r
      atimes (m, r, t, atimes_param);
      // pap = (t, r)
      pap = my_ddot (m, t, 1, r, 1);
      // rur = (t, t)
      double rur = my_ddot (m, t, 1, t, 1);

      zeta = pap / rur;
      // x = x + zeta * r
      my_daxpy (m, zeta, r, 1, x, 1);

      // r = r - zeta * t
      my_daxpy (m, -zeta, t, 1, r, 1);
      // rho1 = (r0, r)
      double rho1 = my_ddot (m, r0, 1, r, 1);

      beta = (rho1 / rho0) * (alpha / zeta);
      rho0 = rho1;
    }

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H


  free (r0);
  free (p);
  free (q);
  free (t);
  free (r);

  if (it->debug == 1)
    {
      fprintf (it->out, "libiter-sta %d %e\n", iter, res2);
    }

  it->niter = iter;
  it->res2  = res2 / b2;
  return (ret);
}

/* bi-cgstab method with precondition
 * INPUT
 *   m : dimension of the problem
 *   b [m] : r-h-s vector
 *   atimes (int m, static double *x, double *b, void *param) :
 *        calc matrix-vector product A.x = b.
 *   atimes_param : parameters for atimes().
 *   inv (int m, static double *b, double *x, void *param) :
 *        approx of A^{-1}.b = x for preconditioning.
 *   inv_param : parameters for the preconditioner inv().
 *   it : struct iter. following entries are used
 *        it->max = kend : max of iteration
 *        it->eps = eps : log10 of cutoff
 * OUTPUT
 *   returned value : 0 == success, otherwise (-1) == failed
 *   x [m] : solution
 *   it->niter : # of iteration
 *   it->res2  : |r^2| / |b^2|
 */
int
sta_pc (int m, const double *b, double *x,
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
  int kend = it->max;
  double eps2 = it->eps * it->eps;

  double *r0 = (double *)calloc (m, sizeof (double));
  double *p  = (double *)calloc (m, sizeof (double));
  double *q  = (double *)calloc (m, sizeof (double));
  double *t  = (double *)calloc (m, sizeof (double));
  double *r  = (double *)calloc (m, sizeof (double));
  double *Kp = (double *)calloc (m, sizeof (double));
  CHECK_MALLOC (r0,  "sta_pc");
  CHECK_MALLOC (p,   "sta_pc");
  CHECK_MALLOC (q,   "sta_pc");
  CHECK_MALLOC (t,   "sta_pc");
  CHECK_MALLOC (r,   "sta_pc");
  CHECK_MALLOC (Kp,  "sta_pc");

  double res2 = 0.0; // for compiler warning.
  double beta = 0.0;
  double zeta = 0.0;

#ifdef HAVE_CBLAS_H
  /**
   * ATLAS version
   */

  double b2 = cblas_ddot (m, b, 1, b, 1); // (b,b)
  eps2 *= b2;

  // r0 = b - A.x
  atimes (m, x, r0, atimes_param);
  cblas_daxpy (m, -1.0, b, 1, r0, 1);
  cblas_dscal (m, -1.0, r0, 1);

  // r  = r0
  cblas_dcopy (m, r0, 1, r, 1);

  // rho0 = (r0, r0)
  double rho0 = cblas_ddot (m, r0, 1, r0, 1);

  int iter;
  for (iter = 0; iter<= kend; iter++)
    {
      // res2 = (r, r)
      res2 = cblas_ddot (m, r, 1, r, 1);
      if (it->debug == 2)
	{
	  fprintf (it->out, "libiter-sta_pc(atlas) %d %e\n", iter, res2 / b2);
	}
      if(res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}

      // p = r + beta * (p - zeta * q)
      cblas_daxpy (m, -zeta, q, 1, p, 1); // p = p - zeta * q
      cblas_dscal (m, beta, p, 1);         // p = beta * (p - zeta * q)
      cblas_daxpy (m, 1.0, r, 1, p, 1);    // p = r + beta * (p - zeta * q)

      // Kp = K^{-1}.p
      inv (m, p, Kp, inv_param);
      // q = A.K^{-1}.p
      atimes (m, Kp, q, atimes_param);
      // pap = (r0, A.K^{-1}.p)
      double pap = cblas_ddot (m, r0, 1, q, 1);

      double alpha = rho0 / pap;
      // r = r - alpha * q
      cblas_daxpy (m, -alpha, q, 1, r, 1);
      // x = x + alpha * K^{-1}p
      cblas_daxpy (m, +alpha, Kp, 1, x, 1);

      // use Kp[] for K^{-1}.r
      inv (m, r, Kp, inv_param);
      // t = A.K^{-1}.r
      atimes (m, Kp, t, atimes_param);
      // pap = (t, r)
      pap = cblas_ddot (m, t, 1, r, 1);
      // rur = (t, t)
      double rur = cblas_ddot (m, t, 1, t, 1);

      zeta = pap / rur;
      // x = x + zeta * K^{-1}.r, where K^{-1}.r is in Kp[]
      cblas_daxpy (m, zeta, Kp, 1, x, 1);

      // r = r - zeta * t
      cblas_daxpy (m, -zeta, t, 1, r, 1);
      // rho1 = (r0, r)
      double rho1 = cblas_ddot (m, r0, 1, r, 1);

      beta = (rho1 / rho0) * (alpha / zeta);
      rho0 = rho1;
    }

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /**
   * BLAS version
   */

  double b2 = ddot_ (&m, b, &i_1, b, &i_1); // (b,b)
  eps2 *= b2;

  // r0 = b - A.x
  atimes (m, x, r0, atimes_param);
  daxpy_ (&m, &d_m1, b, &i_1, r0, &i_1);
  dscal_ (&m, &d_m1, r0, &i_1);

  // r  = r0
  dcopy_ (&m, r0, &i_1, r, &i_1);

  // rho0 = (r0, r0)
  double rho0 = ddot_ (&m, r0, &i_1, r0, &i_1);

  double mzeta = -zeta;
  int iter;
  for (iter = 0; iter<= kend; iter++)
    {
      // res2 = (r, r)
      res2 = ddot_ (&m, r, &i_1, r, &i_1);
      if (it->debug == 2)
	{
	  fprintf (it->out, "libiter-sta_pc(blas) %d %e\n", iter, res2 / b2);
	}
      if(res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}

      // p = r + beta * (p - zeta * q)
      daxpy_ (&m, &mzeta, q, &i_1, p, &i_1); // p = p - zeta * q
      dscal_ (&m, &beta, p, &i_1);           // p = beta * (p - zeta * q)
      daxpy_ (&m, &d_1, r, &i_1, p, &i_1);   // p = r + beta * (p - zeta * q)

      // Kp = K^{-1}.p
      inv (m, p, Kp, inv_param);
      // q = A.K^{-1}.p
      atimes (m, Kp, q, atimes_param);
      // pap = (r0, A.K^{-1}.p)
      double pap = ddot_ (&m, r0, &i_1, q, &i_1);

      double alpha = rho0 / pap;
      // r = r - alpha * q
      double malpha = -alpha;
      daxpy_ (&m, &malpha, q, &i_1, r, &i_1);
      // x = x + alpha * K^{-1}p
      daxpy_ (&m, &alpha, Kp, &i_1, x, &i_1);

      // use Kp[] for K^{-1}.r
      inv (m, r, Kp, inv_param);
      // t = A.K^{-1}.r
      atimes (m, Kp, t, atimes_param);
      // pap = (t, r)
      pap = ddot_ (&m, t, &i_1, r, &i_1);
      // rur = (t, t)
      double rur = ddot_ (&m, t, &i_1, t, &i_1);

      zeta = pap / rur;
      // x = x + zeta * K^{-1}.r, where K^{-1}.r is in Kp[]
      daxpy_ (&m, &zeta, Kp, &i_1, x, &i_1);

      // r = r - zeta * t
      mzeta = -zeta;
      daxpy_ (&m, &mzeta, t, &i_1, r, &i_1);
      // rho1 = (r0, r)
      double rho1 = ddot_ (&m, r0, &i_1, r, &i_1);

      beta = (rho1 / rho0) * (alpha / zeta);
      rho0 = rho1;
    }

# else // !HAVE_BLAS_H
  /**
   * local BLAS version
   */

  double b2 = my_ddot (m, b, 1, b, 1); // (b,b)
  eps2 *= b2;

  // r0 = b - A.x
  atimes (m, x, r0, atimes_param);
  my_daxpy (m, -1.0, b, 1, r0, 1);
  my_dscal (m, -1.0, r0, 1);

  // r  = r0
  my_dcopy (m, r0, 1, r, 1);

  // rho0 = (r0, r0)
  double rho0 = my_ddot (m, r0, 1, r0, 1);

  int iter;
  for (iter = 0; iter<= kend; iter++)
    {
      // res2 = (r, r)
      res2 = my_ddot (m, r, 1, r, 1);
      if (it->debug == 2)
	{
	  fprintf (it->out, "libiter-sta_pc(myblas) %d %e\n", iter, res2 / b2);
	}
      if(res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}

      // p = r + beta * (p - zeta * q)
      my_daxpy (m, -zeta, q, 1, p, 1); // p = p - zeta * q
      my_dscal (m, beta, p, 1);         // p = beta * (p - zeta * q)
      my_daxpy (m, 1.0, r, 1, p, 1);    // p = r + beta * (p - zeta * q)

      // Kp = K^{-1}.p
      inv (m, p, Kp, inv_param);
      // q = A.K^{-1}.p
      atimes (m, Kp, q, atimes_param);
      // pap = (r0, A.K^{-1}.p)
      double pap = my_ddot (m, r0, 1, q, 1);

      double alpha = rho0 / pap;
      // r = r - alpha * q
      my_daxpy (m, -alpha, q, 1, r, 1);
      // x = x + alpha * K^{-1}p
      my_daxpy (m, +alpha, Kp, 1, x, 1);

      // use Kp[] for K^{-1}.r
      inv (m, r, Kp, inv_param);
      // t = A.K^{-1}.r
      atimes (m, Kp, t, atimes_param);
      // pap = (t, r)
      pap = my_ddot (m, t, 1, r, 1);
      // rur = (t, t)
      double rur = my_ddot (m, t, 1, t, 1);

      zeta = pap / rur;
      // x = x + zeta * K^{-1}.r, where K^{-1}.r is in Kp[]
      my_daxpy (m, zeta, Kp, 1, x, 1);

      // r = r - zeta * t
      my_daxpy (m, -zeta, t, 1, r, 1);
      // rho1 = (r0, r)
      double rho1 = my_ddot (m, r0, 1, r, 1);

      beta = (rho1 / rho0) * (alpha / zeta);
      rho0 = rho1;
    }

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H


  free (r0);
  free (p);
  free (q);
  free (t);
  free (r);
  free (Kp);

  if (it->debug == 1)
    {
      fprintf (it->out, "libiter-sta_pc %d %e\n", iter, res2);
    }

  it->niter = iter;
  it->res2  = res2 / b2;
  return (ret);
}


/* bi-cgstab2 method
 * INPUT
 *   m : dimension of the problem
 *   b [m] : r-h-s vector
 *   atimes (int m, static double *x, double *b, void *param) :
 *        calc matrix-vector product A.x = b.
 *   atimes_param : parameters for atimes().
 *   it : struct iter. following entries are used
 *        it->max = kend : max of iteration
 *        it->eps = eps : log10 of cutoff
 * OUTPUT
 *   returned value : 0 == success, otherwise (-1) == failed
 *   x [m] : solution
 *   it->niter : # of iteration
 *   it->res2  : |r^2| / |b^2|
 */
int
sta2 (int m, const double *b, double *x,
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
  int kend = it->max;
  double eps2 = it->eps * it->eps;

  double *r0  = (double *)calloc (m, sizeof (double));
  double *w   = (double *)calloc (m, sizeof (double));
  double *q   = (double *)calloc (m, sizeof (double));
  double *u   = (double *)calloc (m, sizeof (double));
  double *z   = (double *)calloc (m, sizeof (double));
  double *y   = (double *)calloc (m, sizeof (double));
  double *r   = (double *)calloc (m, sizeof (double));
  double *p   = (double *)calloc (m, sizeof (double));
  double *tmp = (double *)calloc (m, sizeof (double));
  CHECK_MALLOC (r0,  "sta2");
  CHECK_MALLOC (w,   "sta2");
  CHECK_MALLOC (q,   "sta2");
  CHECK_MALLOC (u,   "sta2");
  CHECK_MALLOC (z,   "sta2");
  CHECK_MALLOC (y,   "sta2");
  CHECK_MALLOC (r,   "sta2");
  CHECK_MALLOC (p,   "sta2");
  CHECK_MALLOC (tmp, "sta2");

  double res2 = 0.0; // for compiler warning.
  double beta = 0.0;

#ifdef HAVE_CBLAS_H
  /**
   * ATLAS version
   */

  double b2 = cblas_ddot (m, b, 1, b, 1); // (b,b)
  eps2 *= b2;

  // r0 = b - A.x
  atimes (m, x, r0, atimes_param);
  cblas_daxpy (m, -1.0, b, 1, r0, 1);
  cblas_dscal (m, -1.0, r0, 1);
  // r = r0
  cblas_dcopy (m, r0, 1, r, 1);
  // y = -r0
  cblas_dcopy (m, r0, 1, y, 1);
  cblas_dscal (m, -1.0, y, 1);
  // rho0 = (r0, r0)
  double rho0 = cblas_ddot (m, r0, 1, r0, 1);

  int iter;
  for (iter = 0; iter<= kend; iter++)
    {
      // res2 = (r, r)
      res2 = cblas_ddot (m, r, 1, r, 1);
      if (it->debug == 2)
	{
	  fprintf (it->out, "libiter-sta2(atlas) %d %e\n", iter, res2 / b2);
	}
      if(res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}

      // p = r + beta * (p - u)
      cblas_daxpy (m, -1.0, u, 1, p, 1); // p = p - u
      cblas_dscal (m, beta, p, 1);       // p = beta * (p - u)
      cblas_daxpy (m, 1.0, r, 1, p, 1);  // p = r + beta * (p - u)
      // u = y + beta * u
      cblas_dscal (m, beta, u, 1);       // u = beta * u
      cblas_daxpy (m, 1.0, y, 1, u, 1);  // u = y + beta * u

      // q = A.p
      atimes (m, p, q, atimes_param);
      // pap = (r0, q)
      double pap = cblas_ddot (m, r0, 1, q, 1);

      double alpha = rho0 / pap;

      // y = y + alpha * (q - w)
      cblas_dcopy (m, q, 1, tmp, 1);        // tmp = q
      cblas_daxpy (m, -1.0, w, 1, tmp, 1);  // tmp = q - w
      cblas_daxpy (m, alpha, tmp, 1, y, 1); // y += alpha * (q - w)
      // r = r - alpha * q
      cblas_daxpy (m, -alpha, q, 1, r, 1);
      // x = x + alpha * p
      cblas_daxpy (m, alpha, p, 1, x, 1);

      atimes (m, r, w, atimes_param); // w [] -> wn = A tn , tn <- r []
      // AA = (A tn, A tn)
      // yy = (  yn,   yn)
      // Ay = (A tn,   yn)
      // At = (A tn,   tn)
      // yt = (  yn,   tn)
      double AA = cblas_ddot (m, w, 1, w, 1);
      double yy = cblas_ddot (m, y, 1, y, 1);
      double Ay = cblas_ddot (m, w, 1, y, 1);
      double At = cblas_ddot (m, w, 1, r, 1);
      double yt = cblas_ddot (m, y, 1, r, 1);

      double eta;
      double zeta;
      // this part is the only difference from GPBi-CG method (gpb() beclow)
      if ((iter % 2) == 0)
	{
	  zeta = At / AA;
	  eta  = 0.0;
        }
      else
	{
	  zeta = (yy * At - Ay * yt) / (AA * yy - Ay * Ay);
	  eta  = (AA * yt - At * Ay) / (AA * yy - Ay * Ay);
	}

      // z = zeta * r + eta * (z - alpha * u)
      cblas_daxpy (m, -alpha, u, 1, z, 1); // z = z - alpha * u
      cblas_dscal (m, eta, z, 1);          // z = eta (z - alpha * u)
      cblas_daxpy (m, zeta, r, 1, z, 1);   // z = zeta * r + eta (z - alpha * u)
      // x = x + z
      cblas_daxpy (m, 1.0, z, 1, x, 1);

      // y = eta * y + zeta * w
      cblas_dscal (m, eta, y, 1);        // y = eta * y
      cblas_daxpy (m, zeta, w, 1, y, 1); // y = eta * y + zeta * w
      // r = r - y
      cblas_daxpy (m, -1.0, y, 1, r, 1);
      // u = zeta * q + eta * u
      cblas_dscal (m, eta, u, 1);        // u = eta * u
      cblas_daxpy (m, zeta, q, 1, u, 1); // u = eta * u + zeta * q
      // rho1 = (r0, r)
      double rho1 = cblas_ddot (m, r0, 1, r, 1);

      beta = (rho1 / rho0) * (alpha / zeta);
      rho0 = rho1;
      // w = w + beta * q
      cblas_daxpy (m, beta, q, 1, w, 1);
    }

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /**
   * BLAS version
   */

  double b2 = ddot_ (&m, b, &i_1, b, &i_1); // (b,b)
  eps2 *= b2;

  // r0 = b - A.x
  atimes (m, x, r0, atimes_param);
  daxpy_ (&m, &d_m1, b, &i_1, r0, &i_1);
  dscal_ (&m, &d_m1, r0, &i_1);
  // r = r0
  dcopy_ (&m, r0, &i_1, r, &i_1);
  // y = -r0
  dcopy_ (&m, r0, &i_1, y, &i_1);
  dscal_ (&m, &d_m1, y, &i_1);
  // rho0 = (r0, r0)
  double rho0 = ddot_ (&m, r0, &i_1, r0, &i_1);

  int iter;
  for (iter = 0; iter<= kend; iter++)
    {
      // res2 = (r, r)
      res2 = ddot_ (&m, r, &i_1, r, &i_1);
      if (it->debug == 2)
	{
	  fprintf (it->out, "libiter-sta2(blas) %d %e\n", iter, res2 / b2);
	}
      if(res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}

      // p = r + beta * (p - u)
      daxpy_ (&m, &d_m1, u, &i_1, p, &i_1); // p = p - u
      dscal_ (&m, &beta, p, &i_1);          // p = beta * (p - u)
      daxpy_ (&m, &d_1, r, &i_1, p, &i_1);  // p = r + beta * (p - u)
      // u = y + beta * u
      dscal_ (&m, &beta, u, &i_1);       // u = beta * u
      daxpy_ (&m, &d_1, y, &i_1, u, &i_1);  // u = y + beta * u

      // q = A.p
      atimes (m, p, q, atimes_param);
      // pap = (r0, q)
      double pap = ddot_ (&m, r0, &i_1, q, &i_1);

      double alpha = rho0 / pap;

      // y = y + alpha * (q - w)
      dcopy_ (&m, q, &i_1, tmp, &i_1);         // tmp = q
      daxpy_ (&m, &d_m1, w, &i_1, tmp, &i_1);  // tmp = q - w
      daxpy_ (&m, &alpha, tmp, &i_1, y, &i_1); // y += alpha * (q - w)
      // r = r - alpha * q
      double malpha = -alpha;
      daxpy_ (&m, &malpha, q, &i_1, r, &i_1);
      // x = x + alpha * p
      daxpy_ (&m, &alpha, p, &i_1, x, &i_1);

      atimes (m, r, w, atimes_param); // w [] -> wn = A tn , tn <- r []
      // AA = (A tn, A tn)
      // yy = (  yn,   yn)
      // Ay = (A tn,   yn)
      // At = (A tn,   tn)
      // yt = (  yn,   tn)
      double AA = ddot_ (&m, w, &i_1, w, &i_1);
      double yy = ddot_ (&m, y, &i_1, y, &i_1);
      double Ay = ddot_ (&m, w, &i_1, y, &i_1);
      double At = ddot_ (&m, w, &i_1, r, &i_1);
      double yt = ddot_ (&m, y, &i_1, r, &i_1);

      double eta;
      double zeta;
      // this part is the only difference from GPBi-CG method (gpb() beclow)
      if ((iter % 2) == 0)
	{
	  zeta = At / AA;
	  eta  = 0.0;
        }
      else
	{
	  zeta = (yy * At - Ay * yt) / (AA * yy - Ay * Ay);
	  eta  = (AA * yt - At * Ay) / (AA * yy - Ay * Ay);
	}

      // z = zeta * r + eta * (z - alpha * u)
      daxpy_ (&m, &malpha, u, &i_1, z, &i_1);// z = z - alpha * u
      dscal_ (&m, &eta, z, &i_1);          // z = eta (z - alpha * u)
      daxpy_ (&m, &zeta, r, &i_1, z, &i_1);// z = zeta * r + eta (z - alpha * u)
      // x = x + z
      daxpy_ (&m, &d_1, z, &i_1, x, &i_1);

      // y = eta * y + zeta * w
      dscal_ (&m, &eta, y, &i_1);           // y = eta * y
      daxpy_ (&m, &zeta, w, &i_1, y, &i_1); // y = eta * y + zeta * w
      // r = r - y
      daxpy_ (&m, &d_m1, y, &i_1, r, &i_1);
      // u = zeta * q + eta * u
      dscal_ (&m, &eta, u, &i_1);           // u = eta * u
      daxpy_ (&m, &zeta, q, &i_1, u, &i_1); // u = eta * u + zeta * q
      // rho1 = (r0, r)
      double rho1 = ddot_ (&m, r0, &i_1, r, &i_1);

      beta = (rho1 / rho0) * (alpha / zeta);
      rho0 = rho1;
      // w = w + beta * q
      daxpy_ (&m, &beta, q, &i_1, w, &i_1);
    }

# else // !HAVE_BLAS_H
  /**
   * local BLAS version
   */

  double b2 = my_ddot (m, b, 1, b, 1); // (b,b)
  eps2 *= b2;

  // r0 = b - A.x
  atimes (m, x, r0, atimes_param);
  my_daxpy (m, -1.0, b, 1, r0, 1);
  my_dscal (m, -1.0, r0, 1);
  // r = r0
  my_dcopy (m, r0, 1, r, 1);
  // y = -r0
  my_dcopy (m, r0, 1, y, 1);
  my_dscal (m, -1.0, y, 1);
  // rho0 = (r0, r0)
  double rho0 = my_ddot (m, r0, 1, r0, 1);

  int iter;
  for (iter = 0; iter<= kend; iter++)
    {
      // res2 = (r, r)
      res2 = my_ddot (m, r, 1, r, 1);
      if (it->debug == 2)
	{
	  fprintf (it->out, "libiter-sta2(myblas) %d %e\n", iter, res2 / b2);
	}
      if(res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}

      // p = r + beta * (p - u)
      my_daxpy (m, -1.0, u, 1, p, 1); // p = p - u
      my_dscal (m, beta, p, 1);       // p = beta * (p - u)
      my_daxpy (m, 1.0, r, 1, p, 1);  // p = r + beta * (p - u)
      // u = y + beta * u
      my_dscal (m, beta, u, 1);       // u = beta * u
      my_daxpy (m, 1.0, y, 1, u, 1);  // u = y + beta * u

      // q = A.p
      atimes (m, p, q, atimes_param);
      // pap = (r0, q)
      double pap = my_ddot (m, r0, 1, q, 1);

      double alpha = rho0 / pap;

      // y = y + alpha * (q - w)
      my_dcopy (m, q, 1, tmp, 1);        // tmp = q
      my_daxpy (m, -1.0, w, 1, tmp, 1);  // tmp = q - w
      my_daxpy (m, alpha, tmp, 1, y, 1); // y += alpha * (q - w)
      // r = r - alpha * q
      my_daxpy (m, -alpha, q, 1, r, 1);
      // x = x + alpha * p
      my_daxpy (m, alpha, p, 1, x, 1);

      atimes (m, r, w, atimes_param); // w [] -> wn = A tn , tn <- r []
      // AA = (A tn, A tn)
      // yy = (  yn,   yn)
      // Ay = (A tn,   yn)
      // At = (A tn,   tn)
      // yt = (  yn,   tn)
      double AA = my_ddot (m, w, 1, w, 1);
      double yy = my_ddot (m, y, 1, y, 1);
      double Ay = my_ddot (m, w, 1, y, 1);
      double At = my_ddot (m, w, 1, r, 1);
      double yt = my_ddot (m, y, 1, r, 1);

      double eta;
      double zeta;
      // this part is the only difference from GPBi-CG method (gpb() beclow)
      if ((iter % 2) == 0)
	{
	  zeta = At / AA;
	  eta  = 0.0;
        }
      else
	{
	  zeta = (yy * At - Ay * yt) / (AA * yy - Ay * Ay);
	  eta  = (AA * yt - At * Ay) / (AA * yy - Ay * Ay);
	}

      // z = zeta * r + eta * (z - alpha * u)
      my_daxpy (m, -alpha, u, 1, z, 1); // z = z - alpha * u
      my_dscal (m, eta, z, 1);          // z = eta (z - alpha * u)
      my_daxpy (m, zeta, r, 1, z, 1);   // z = zeta * r + eta (z - alpha * u)
      // x = x + z
      my_daxpy (m, 1.0, z, 1, x, 1);

      // y = eta * y + zeta * w
      my_dscal (m, eta, y, 1);        // y = eta * y
      my_daxpy (m, zeta, w, 1, y, 1); // y = eta * y + zeta * w
      // r = r - y
      my_daxpy (m, -1.0, y, 1, r, 1);
      // u = zeta * q + eta * u
      my_dscal (m, eta, u, 1);        // u = eta * u
      my_daxpy (m, zeta, q, 1, u, 1); // u = eta * u + zeta * q
      // rho1 = (r0, r)
      double rho1 = my_ddot (m, r0, 1, r, 1);

      beta = (rho1 / rho0) * (alpha / zeta);
      rho0 = rho1;
      // w = w + beta * q
      my_daxpy (m, beta, q, 1, w, 1);
    }

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H

  free (r0);
  free (w);
  free (q);
  free (u);
  free (z);
  free (y);
  free (r);
  free (p);
  free (tmp);

  if (it->debug == 1)
    {
      fprintf (it->out, "libiter-st2 %d %e\n", iter, res2);
    }

  it->niter = iter;
  it->res2  = res2 / b2;
  return (ret);
}

/* bi-cgstab2 method with precondition
 * INPUT
 *   m : dimension of the problem
 *   b [m] : r-h-s vector
 *   atimes (int m, static double *x, double *b, void *param) :
 *        calc matrix-vector product A.x = b.
 *   atimes_param : parameters for atimes().
 *   inv (int m, static double *b, double *x, void *param) :
 *        approx of A^{-1}.b = x for preconditioning.
 *   inv_param : parameters for the preconditioner inv().
 *   it : struct iter. following entries are used
 *        it->max = kend : max of iteration
 *        it->eps = eps : log10 of cutoff
 * OUTPUT
 *   returned value : 0 == success, otherwise (-1) == failed
 *   x [m] : solution
 *   it->niter : # of iteration
 *   it->res2  : |r^2| / |b^2|
 */
int
sta2_pc (int m, const double *b, double *x,
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
  int kend = it->max;
  double eps2 = it->eps * it->eps;

  double *r0  = (double *)calloc (m, sizeof (double));
  double *w   = (double *)calloc (m, sizeof (double));
  double *q   = (double *)calloc (m, sizeof (double));
  double *u   = (double *)calloc (m, sizeof (double));
  double *z   = (double *)calloc (m, sizeof (double));
  double *y   = (double *)calloc (m, sizeof (double));
  double *r   = (double *)calloc (m, sizeof (double));
  double *p   = (double *)calloc (m, sizeof (double));
  double *tmp = (double *)calloc (m, sizeof (double));
  double *Kr  = (double *)calloc (m, sizeof (double));
  CHECK_MALLOC (r0,  "sta2_pc");
  CHECK_MALLOC (w,   "sta2_pc");
  CHECK_MALLOC (q,   "sta2_pc");
  CHECK_MALLOC (u,   "sta2_pc");
  CHECK_MALLOC (z,   "sta2_pc");
  CHECK_MALLOC (y,   "sta2_pc");
  CHECK_MALLOC (r,   "sta2_pc");
  CHECK_MALLOC (p,   "sta2_pc");
  CHECK_MALLOC (tmp, "sta2_pc");
  CHECK_MALLOC (Kr,  "sta2_pc");

  double res2 = 0.0; // for compiler warning.
  double beta = 0.0;

#ifdef HAVE_CBLAS_H
  /**
   * ATLAS version
   */

  double b2 = cblas_ddot (m, b, 1, b, 1); // (b,b)
  eps2 *= b2;

  // r0 = b - A.x
  atimes (m, x, r0, atimes_param);
  cblas_daxpy (m, -1.0, b, 1, r0, 1);
  cblas_dscal (m, -1.0, r0, 1);
  // r = r0
  cblas_dcopy (m, r0, 1, r, 1);
  // y = -r0
  cblas_dcopy (m, r0, 1, y, 1);
  cblas_dscal (m, -1.0, y, 1);
  // rho0 = (r0, r0)
  double rho0 = cblas_ddot (m, r0, 1, r0, 1);

  int iter;
  for (iter = 0; iter<= kend; iter++)
    {
      // res2 = (r, r)
      res2 = cblas_ddot (m, r, 1, r, 1);
      if (it->debug == 2)
	{
	  fprintf (it->out, "libiter-sta2_pc(atlas) %d %e\n", iter, res2 / b2);
	}
      if(res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}

      // Kr = K^{-1}.r
      inv (m, r, Kr, inv_param);
      // p = K^{-1}.r + beta * (p - u)
      cblas_daxpy (m, -1.0, u, 1, p, 1); // p = p - u
      cblas_dscal (m, beta, p, 1);       // p = beta * (p - u)
      cblas_daxpy (m, 1.0, Kr, 1, p, 1);// p = K^{-1}.r + beta * (p - u)
      // u = K^{-1}.y + beta * u
      inv (m, y, tmp, inv_param);       // tmp = K^{-1}.y
      cblas_dscal (m, beta, u, 1);      // u = beta * u
      cblas_daxpy (m, 1.0, tmp, 1, u, 1); // u = K^{-1}.y + beta * u

      // q = A.p
      atimes (m, p, q, atimes_param);
      // pap = (r0, q)
      double pap = cblas_ddot (m, r0, 1, q, 1);

      double alpha = rho0 / pap;
      // y = y + alpha * (q - w)
      cblas_dcopy (m, q, 1, tmp, 1);        // tmp = q
      cblas_daxpy (m, -1.0, w, 1, tmp, 1);  // tmp = q - w
      cblas_daxpy (m, alpha, tmp, 1, y, 1); // y += alpha * (q - w)
      // r = r - alpha * q, note that at this point this is t_n
      cblas_daxpy (m, -alpha, q, 1, r, 1);
      // x = x + alpha * p
      cblas_daxpy (m, alpha, p, 1, x, 1);

      // w = A.K^{-1}.r (= A.K^{-1}.tn)
      inv (m, r, Kr, inv_param); // Kr = K^{-1}.tn, not equal to Kr above
      atimes (m, Kr, w, atimes_param);// w = A.K^{-1}.tn
      // AA = (A tn, A tn)
      // yy = (  yn,   yn)
      // Ay = (A tn,   yn)
      // At = (A tn,   tn)
      // yt = (  yn,   tn)
      double AA = cblas_ddot (m, w, 1, w, 1);
      double yy = cblas_ddot (m, y, 1, y, 1);
      double Ay = cblas_ddot (m, w, 1, y, 1);
      double At = cblas_ddot (m, w, 1, r, 1);
      double yt = cblas_ddot (m, y, 1, r, 1);

      double eta;
      double zeta;
      // this part is the only difference from GPBi-CG method (gpb() beclow)
      if ((iter % 2) == 0)
	{
	  zeta = At / AA;
	  eta  = 0.0;
        }
      else
	{
	  zeta = (yy * At - Ay * yt) / (AA * yy - Ay * Ay);
	  eta  = (AA * yt - At * Ay) / (AA * yy - Ay * Ay);
	}

      // z = zeta * K^{-1}.r + eta * (z - alpha * u)
      cblas_daxpy (m, -alpha, u, 1, z, 1);// z = z - alpha * u
      cblas_dscal (m, eta, z, 1);         // z = eta (z - alpha * u)
      cblas_daxpy (m, zeta, Kr, 1, z, 1); // z = zeta * Kr + eta (z - alpha * u)
      // x = x + z
      cblas_daxpy (m, 1.0, z, 1, x, 1);

      // y = eta * y + zeta * w
      cblas_dscal (m, eta, y, 1);        // y = eta * y
      cblas_daxpy (m, zeta, w, 1, y, 1); // y = eta * y + zeta * w
      // r = r - y
      cblas_daxpy (m, -1.0, y, 1, r, 1);
      // u = zeta * q + eta * u
      inv (m, q, tmp, inv_param);        // tmp = K^{-1}.q
      cblas_dscal (m, eta, u, 1);        // u = eta * u
      cblas_daxpy (m, zeta, tmp, 1, u, 1); // u = eta * u + zeta * K^{-1}.q
      // rho1 = (r0, r)
      double rho1 = cblas_ddot (m, r0, 1, r, 1);

      beta = (rho1 / rho0) * (alpha / zeta);
      rho0 = rho1;
      // w = w + beta * q
      cblas_daxpy (m, beta, q, 1, w, 1);
    }

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /**
   * BLAS version
   */

  double b2 = ddot_ (&m, b, &i_1, b, &i_1); // (b,b)
  eps2 *= b2;

  // r0 = b - A.x
  atimes (m, x, r0, atimes_param);
  daxpy_ (&m, &d_m1, b, &i_1, r0, &i_1);
  dscal_ (&m, &d_m1, r0, &i_1);
  // r = r0
  dcopy_ (&m, r0, &i_1, r, &i_1);
  // y = -r0
  dcopy_ (&m, r0, &i_1, y, &i_1);
  dscal_ (&m, &d_m1, y, &i_1);
  // rho0 = (r0, r0)
  double rho0 = ddot_ (&m, r0, &i_1, r0, &i_1);

  int iter;
  for (iter = 0; iter<= kend; iter++)
    {
      // res2 = (r, r)
      res2 = ddot_ (&m, r, &i_1, r, &i_1);
      if (it->debug == 2)
	{
	  fprintf (it->out, "libiter-sta2_pc(blas) %d %e\n", iter, res2 / b2);
	}
      if(res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}

      // Kr = K^{-1}.r
      inv (m, r, Kr, inv_param);
      // p = K^{-1}.r + beta * (p - u)
      daxpy_ (&m, &d_m1, u, &i_1, p, &i_1); // p = p - u
      dscal_ (&m, &beta, p, &i_1);         // p = beta * (p - u)
      daxpy_ (&m, &d_1, Kr, &i_1, p, &i_1);// p = K^{-1}.r + beta * (p - u)
      // u = K^{-1}.y + beta * u
      inv (m, y, tmp, inv_param);            // tmp = K^{-1}.y
      dscal_ (&m, &beta, u, &i_1);           // u = beta * u
      daxpy_ (&m, &d_1, tmp, &i_1, u, &i_1); // u = K^{-1}.y + beta * u

      // q = A.p
      atimes (m, p, q, atimes_param);
      // pap = (r0, q)
      double pap = ddot_ (&m, r0, &i_1, q, &i_1);

      double alpha = rho0 / pap;
      // y = y + alpha * (q - w)
      dcopy_ (&m, q, &i_1, tmp, &i_1);        // tmp = q
      daxpy_ (&m, &d_m1, w, &i_1, tmp, &i_1);  // tmp = q - w
      daxpy_ (&m, &alpha, tmp, &i_1, y, &i_1); // y += alpha * (q - w)
      // r = r - alpha * q, note that at this point this is t_n
      double malpha = -alpha;
      daxpy_ (&m, &malpha, q, &i_1, r, &i_1);
      // x = x + alpha * p
      daxpy_ (&m, &alpha, p, &i_1, x, &i_1);

      // w = A.K^{-1}.r (= A.K^{-1}.tn)
      inv (m, r, Kr, inv_param); // Kr = K^{-1}.tn, not equal to Kr above
      atimes (m, Kr, w, atimes_param);// w = A.K^{-1}.tn
      // AA = (A tn, A tn)
      // yy = (  yn,   yn)
      // Ay = (A tn,   yn)
      // At = (A tn,   tn)
      // yt = (  yn,   tn)
      double AA = ddot_ (&m, w, &i_1, w, &i_1);
      double yy = ddot_ (&m, y, &i_1, y, &i_1);
      double Ay = ddot_ (&m, w, &i_1, y, &i_1);
      double At = ddot_ (&m, w, &i_1, r, &i_1);
      double yt = ddot_ (&m, y, &i_1, r, &i_1);

      double eta;
      double zeta;
      // this part is the only difference from GPBi-CG method (gpb() beclow)
      if ((iter % 2) == 0)
	{
	  zeta = At / AA;
	  eta  = 0.0;
        }
      else
	{
	  zeta = (yy * At - Ay * yt) / (AA * yy - Ay * Ay);
	  eta  = (AA * yt - At * Ay) / (AA * yy - Ay * Ay);
	}

      // z = zeta * K^{-1}.r + eta * (z - alpha * u)
      daxpy_ (&m, &malpha, u, &i_1, z, &i_1);// z = z - alpha * u
      dscal_ (&m, &eta, z, &i_1);           // z = eta (z - alpha * u)
      daxpy_ (&m, &zeta, Kr, &i_1, z, &i_1);// z = zeta*Kr + eta(z - alpha*u)
      // x = x + z
      daxpy_ (&m, &d_1, z, &i_1, x, &i_1);

      // y = eta * y + zeta * w
      dscal_ (&m, &eta, y, &i_1);           // y = eta * y
      daxpy_ (&m, &zeta, w, &i_1, y, &i_1); // y = eta * y + zeta * w
      // r = r - y
      daxpy_ (&m, &d_m1, y, &i_1, r, &i_1);
      // u = zeta * q + eta * u
      inv (m, q, tmp, inv_param);        // tmp = K^{-1}.q
      dscal_ (&m, &eta, u, &i_1);        // u = eta * u
      daxpy_ (&m, &zeta, tmp, &i_1, u, &i_1); // u = eta * u + zeta * K^{-1}.q
      // rho1 = (r0, r)
      double rho1 = ddot_ (&m, r0, &i_1, r, &i_1);

      beta = (rho1 / rho0) * (alpha / zeta);
      rho0 = rho1;
      // w = w + beta * q
      daxpy_ (&m, &beta, q, &i_1, w, &i_1);
    }

# else // !HAVE_BLAS_H
  /**
   * local BLAS version
   */

  double b2 = my_ddot (m, b, 1, b, 1); // (b,b)
  eps2 *= b2;

  // r0 = b - A.x
  atimes (m, x, r0, atimes_param);
  my_daxpy (m, -1.0, b, 1, r0, 1);
  my_dscal (m, -1.0, r0, 1);
  // r = r0
  my_dcopy (m, r0, 1, r, 1);
  // y = -r0
  my_dcopy (m, r0, 1, y, 1);
  my_dscal (m, -1.0, y, 1);
  // rho0 = (r0, r0)
  double rho0 = my_ddot (m, r0, 1, r0, 1);

  int iter;
  for (iter = 0; iter<= kend; iter++)
    {
      // res2 = (r, r)
      res2 = my_ddot (m, r, 1, r, 1);
      if (it->debug == 2)
	{
	  fprintf (it->out, "libiter-sta2_pc(myblas) %d %e\n", iter, res2 / b2);
	}
      if(res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}

      // Kr = K^{-1}.r
      inv (m, r, Kr, inv_param);
      // p = K^{-1}.r + beta * (p - u)
      my_daxpy (m, -1.0, u, 1, p, 1); // p = p - u
      my_dscal (m, beta, p, 1);       // p = beta * (p - u)
      my_daxpy (m, 1.0, Kr, 1, p, 1);// p = K^{-1}.r + beta * (p - u)
      // u = K^{-1}.y + beta * u
      inv (m, y, tmp, inv_param);       // tmp = K^{-1}.y
      my_dscal (m, beta, u, 1);      // u = beta * u
      my_daxpy (m, 1.0, tmp, 1, u, 1); // u = K^{-1}.y + beta * u

      // q = A.p
      atimes (m, p, q, atimes_param);
      // pap = (r0, q)
      double pap = my_ddot (m, r0, 1, q, 1);

      double alpha = rho0 / pap;
      // y = y + alpha * (q - w)
      my_dcopy (m, q, 1, tmp, 1);        // tmp = q
      my_daxpy (m, -1.0, w, 1, tmp, 1);  // tmp = q - w
      my_daxpy (m, alpha, tmp, 1, y, 1); // y += alpha * (q - w)
      // r = r - alpha * q, note that at this point this is t_n
      my_daxpy (m, -alpha, q, 1, r, 1);
      // x = x + alpha * p
      my_daxpy (m, alpha, p, 1, x, 1);

      // w = A.K^{-1}.r (= A.K^{-1}.tn)
      inv (m, r, Kr, inv_param); // Kr = K^{-1}.tn, not equal to Kr above
      atimes (m, Kr, w, atimes_param);// w = A.K^{-1}.tn
      // AA = (A tn, A tn)
      // yy = (  yn,   yn)
      // Ay = (A tn,   yn)
      // At = (A tn,   tn)
      // yt = (  yn,   tn)
      double AA = my_ddot (m, w, 1, w, 1);
      double yy = my_ddot (m, y, 1, y, 1);
      double Ay = my_ddot (m, w, 1, y, 1);
      double At = my_ddot (m, w, 1, r, 1);
      double yt = my_ddot (m, y, 1, r, 1);

      double eta;
      double zeta;
      // this part is the only difference from GPBi-CG method (gpb() beclow)
      if ((iter % 2) == 0)
	{
	  zeta = At / AA;
	  eta  = 0.0;
        }
      else
	{
	  zeta = (yy * At - Ay * yt) / (AA * yy - Ay * Ay);
	  eta  = (AA * yt - At * Ay) / (AA * yy - Ay * Ay);
	}

      // z = zeta * K^{-1}.r + eta * (z - alpha * u)
      my_daxpy (m, -alpha, u, 1, z, 1);// z = z - alpha * u
      my_dscal (m, eta, z, 1);         // z = eta (z - alpha * u)
      my_daxpy (m, zeta, Kr, 1, z, 1); // z = zeta * Kr + eta (z - alpha * u)
      // x = x + z
      my_daxpy (m, 1.0, z, 1, x, 1);

      // y = eta * y + zeta * w
      my_dscal (m, eta, y, 1);        // y = eta * y
      my_daxpy (m, zeta, w, 1, y, 1); // y = eta * y + zeta * w
      // r = r - y
      my_daxpy (m, -1.0, y, 1, r, 1);
      // u = zeta * q + eta * u
      inv (m, q, tmp, inv_param);        // tmp = K^{-1}.q
      my_dscal (m, eta, u, 1);        // u = eta * u
      my_daxpy (m, zeta, tmp, 1, u, 1); // u = eta * u + zeta * K^{-1}.q
      // rho1 = (r0, r)
      double rho1 = my_ddot (m, r0, 1, r, 1);

      beta = (rho1 / rho0) * (alpha / zeta);
      rho0 = rho1;
      // w = w + beta * q
      my_daxpy (m, beta, q, 1, w, 1);
    }

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H

  free (r0);
  free (w);
  free (q);
  free (u);
  free (z);
  free (y);
  free (r);
  free (p);
  free (tmp);
  free (Kr);

  if (it->debug == 1)
    {
      fprintf (it->out, "libiter-sta2_pc %d %e\n", iter, res2);
    }

  it->niter = iter;
  it->res2  = res2 / b2;
  return (ret);
}

/* gpbi-cg method
 * ref: Zhang, SIAM J.Sci.Comput. 1997 vol.18 pp.537-551.
 * INPUT
 *   m : dimension of the problem
 *   b [m] : r-h-s vector
 *   atimes (int m, double *x, double *b) : calc matrix-vector product
 *   atimes_param : pointer to be passed to atimes routines
 *   it : struct iter. following entries are used
 *        it->max = kend : max of iteration
 *        it->eps = eps : log10 of cutoff
 * OUTPUT
 *   returned value : 0 == success, otherwise (-1) == failed
 *   x [m] : solution
 *   it->niter : # of iteration
 *   it->res2  : |r^2| / |b^2|
 */
int
gpb (int m, const double *b, double *x,
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
  int kend = it->max;
  double eps2 = it->eps * it->eps;

  double *r0  = (double *)calloc (m, sizeof (double));
  double *w   = (double *)calloc (m, sizeof (double));
  double *q   = (double *)calloc (m, sizeof (double));
  double *u   = (double *)calloc (m, sizeof (double));
  double *z   = (double *)calloc (m, sizeof (double));
  double *y   = (double *)calloc (m, sizeof (double));
  double *r   = (double *)calloc (m, sizeof (double));
  double *p   = (double *)calloc (m, sizeof (double));
  double *tmp = (double *)calloc (m, sizeof (double));
  CHECK_MALLOC (r0,  "gpb");
  CHECK_MALLOC (w,   "gpb");
  CHECK_MALLOC (q,   "gpb");
  CHECK_MALLOC (u,   "gpb");
  CHECK_MALLOC (z,   "gpb");
  CHECK_MALLOC (y,   "gpb");
  CHECK_MALLOC (r,   "gpb");
  CHECK_MALLOC (p,   "gpb");
  CHECK_MALLOC (tmp, "gpb");

  double res2 = 0.0; // for compiler warning.
  double beta = 0.0;

#ifdef HAVE_CBLAS_H
  /**
   * ATLAS version
   */

  double b2 = cblas_ddot (m, b, 1, b, 1); // (b,b)
  eps2 *= b2;

  // r0 = b - A.x
  atimes (m, x, r0, atimes_param);
  cblas_daxpy (m, -1.0, b, 1, r0, 1);
  cblas_dscal (m, -1.0, r0, 1);
  // r = r0
  cblas_dcopy (m, r0, 1, r, 1);
  // y = -r0
  cblas_dcopy (m, r0, 1, y, 1);
  cblas_dscal (m, -1.0, y, 1);
  // rho0 = (r0, r0)
  double rho0 = cblas_ddot (m, r0, 1, r0, 1);

  int iter;
  for (iter = 0; iter<= kend; iter++)
    {
      // res2 = (r, r)
      res2 = cblas_ddot (m, r, 1, r, 1);
      if (it->debug == 2)
	{
	  fprintf (it->out, "libiter-gpb(atlas) %d %e\n", iter, res2 / b2);
	}
      if(res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}

      // p = r + beta * (p - u)
      cblas_daxpy (m, -1.0, u, 1, p, 1); // p = p - u
      cblas_dscal (m, beta, p, 1);       // p = beta * (p - u)
      cblas_daxpy (m, 1.0, r, 1, p, 1);  // p = r + beta * (p - u)
      // u = y + beta * u
      cblas_dscal (m, beta, u, 1);       // u = beta * u
      cblas_daxpy (m, 1.0, y, 1, u, 1);  // u = y + beta * u

      // q = A.p
      atimes (m, p, q, atimes_param);
      // pap = (r0, q)
      double pap = cblas_ddot (m, r0, 1, q, 1);

      double alpha = rho0 / pap;

      // y = y + alpha * (q - w)
      cblas_dcopy (m, q, 1, tmp, 1);        // tmp = q
      cblas_daxpy (m, -1.0, w, 1, tmp, 1);  // tmp = q - w
      cblas_daxpy (m, alpha, tmp, 1, y, 1); // y += alpha * (q - w)
      // r = r - alpha * q
      cblas_daxpy (m, -alpha, q, 1, r, 1);
      // x = x + alpha * p
      cblas_daxpy (m, alpha, p, 1, x, 1);

      atimes (m, r, w, atimes_param); // w [] -> wn = A tn , tn <- r []
      // AA = (A tn, A tn)
      // yy = (  yn,   yn)
      // Ay = (A tn,   yn)
      // At = (A tn,   tn)
      // yt = (  yn,   tn)
      double AA = cblas_ddot (m, w, 1, w, 1);
      double yy = cblas_ddot (m, y, 1, y, 1);
      double Ay = cblas_ddot (m, w, 1, y, 1);
      double At = cblas_ddot (m, w, 1, r, 1);
      double yt = cblas_ddot (m, y, 1, r, 1);

      double eta;
      double zeta;
      // this part is the only difference from Bi-CGSTAB2 (sta2() above)
      if (iter == 0)
	{
	  zeta = At / AA;
	  eta  = 0.0;
        }
      else
	{
	  zeta = (yy * At - Ay * yt) / (AA * yy - Ay * Ay);
	  eta  = (AA * yt - At * Ay) / (AA * yy - Ay * Ay);
	}

      // z = zeta * r + eta * (z - alpha * u)
      cblas_daxpy (m, -alpha, u, 1, z, 1); // z = z - alpha * u
      cblas_dscal (m, eta, z, 1);          // z = eta (z - alpha * u)
      cblas_daxpy (m, zeta, r, 1, z, 1);   // z = zeta * r + eta (z - alpha * u)
      // x = x + z
      cblas_daxpy (m, 1.0, z, 1, x, 1);

      // y = eta * y + zeta * w
      cblas_dscal (m, eta, y, 1);        // y = eta * y
      cblas_daxpy (m, zeta, w, 1, y, 1); // y = eta * y + zeta * w
      // r = r - y
      cblas_daxpy (m, -1.0, y, 1, r, 1);
      // u = zeta * q + eta * u
      cblas_dscal (m, eta, u, 1);        // u = eta * u
      cblas_daxpy (m, zeta, q, 1, u, 1); // u = eta * u + zeta * q
      // rho1 = (r0, r)
      double rho1 = cblas_ddot (m, r0, 1, r, 1);

      beta = (rho1 / rho0) * (alpha / zeta);
      rho0 = rho1;
      // w = w + beta * q
      cblas_daxpy (m, beta, q, 1, w, 1);
    }

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /**
   * BLAS version
   */

  double b2 = ddot_ (&m, b, &i_1, b, &i_1); // (b,b)
  eps2 *= b2;

  // r0 = b - A.x
  atimes (m, x, r0, atimes_param);
  daxpy_ (&m, &d_m1, b, &i_1, r0, &i_1);
  dscal_ (&m, &d_m1, r0, &i_1);
  // r = r0
  dcopy_ (&m, r0, &i_1, r, &i_1);
  // y = -r0
  dcopy_ (&m, r0, &i_1, y, &i_1);
  dscal_ (&m, &d_m1, y, &i_1);
  // rho0 = (r0, r0)
  double rho0 = ddot_ (&m, r0, &i_1, r0, &i_1);

  int iter;
  for (iter = 0; iter<= kend; iter++)
    {
      // res2 = (r, r)
      res2 = ddot_ (&m, r, &i_1, r, &i_1);
      if (it->debug == 2)
	{
	  fprintf (it->out, "libiter-gpb(blas) %d %e\n", iter, res2 / b2);
	}
      if(res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}

      // p = r + beta * (p - u)
      daxpy_ (&m, &d_m1, u, &i_1, p, &i_1); // p = p - u
      dscal_ (&m, &beta, p, &i_1);          // p = beta * (p - u)
      daxpy_ (&m, &d_1, r, &i_1, p, &i_1);  // p = r + beta * (p - u)
      // u = y + beta * u
      dscal_ (&m, &beta, u, &i_1);       // u = beta * u
      daxpy_ (&m, &d_1, y, &i_1, u, &i_1);  // u = y + beta * u

      // q = A.p
      atimes (m, p, q, atimes_param);
      // pap = (r0, q)
      double pap = ddot_ (&m, r0, &i_1, q, &i_1);

      double alpha = rho0 / pap;

      // y = y + alpha * (q - w)
      dcopy_ (&m, q, &i_1, tmp, &i_1);         // tmp = q
      daxpy_ (&m, &d_m1, w, &i_1, tmp, &i_1);  // tmp = q - w
      daxpy_ (&m, &alpha, tmp, &i_1, y, &i_1); // y += alpha * (q - w)
      // r = r - alpha * q
      double malpha = -alpha;
      daxpy_ (&m, &malpha, q, &i_1, r, &i_1);
      // x = x + alpha * p
      daxpy_ (&m, &alpha, p, &i_1, x, &i_1);

      atimes (m, r, w, atimes_param); // w [] -> wn = A tn , tn <- r []
      // AA = (A tn, A tn)
      // yy = (  yn,   yn)
      // Ay = (A tn,   yn)
      // At = (A tn,   tn)
      // yt = (  yn,   tn)
      double AA = ddot_ (&m, w, &i_1, w, &i_1);
      double yy = ddot_ (&m, y, &i_1, y, &i_1);
      double Ay = ddot_ (&m, w, &i_1, y, &i_1);
      double At = ddot_ (&m, w, &i_1, r, &i_1);
      double yt = ddot_ (&m, y, &i_1, r, &i_1);

      double eta;
      double zeta;
      // this part is the only difference from Bi-CGSTAB2 (sta2() above)
      if (iter == 0)
	{
	  zeta = At / AA;
	  eta  = 0.0;
        }
      else
	{
	  zeta = (yy * At - Ay * yt) / (AA * yy - Ay * Ay);
	  eta  = (AA * yt - At * Ay) / (AA * yy - Ay * Ay);
	}

      // z = zeta * r + eta * (z - alpha * u)
      daxpy_ (&m, &malpha, u, &i_1, z, &i_1);// z = z - alpha * u
      dscal_ (&m, &eta, z, &i_1);          // z = eta (z - alpha * u)
      daxpy_ (&m, &zeta, r, &i_1, z, &i_1);// z = zeta * r + eta (z - alpha * u)
      // x = x + z
      daxpy_ (&m, &d_1, z, &i_1, x, &i_1);

      // y = eta * y + zeta * w
      dscal_ (&m, &eta, y, &i_1);           // y = eta * y
      daxpy_ (&m, &zeta, w, &i_1, y, &i_1); // y = eta * y + zeta * w
      // r = r - y
      daxpy_ (&m, &d_m1, y, &i_1, r, &i_1);
      // u = zeta * q + eta * u
      dscal_ (&m, &eta, u, &i_1);           // u = eta * u
      daxpy_ (&m, &zeta, q, &i_1, u, &i_1); // u = eta * u + zeta * q
      // rho1 = (r0, r)
      double rho1 = ddot_ (&m, r0, &i_1, r, &i_1);

      beta = (rho1 / rho0) * (alpha / zeta);
      rho0 = rho1;
      // w = w + beta * q
      daxpy_ (&m, &beta, q, &i_1, w, &i_1);
    }

# else // !HAVE_BLAS_H
  /**
   * local BLAS version
   */

  double b2 = my_ddot (m, b, 1, b, 1); // (b,b)
  eps2 *= b2;

  // r0 = b - A.x
  atimes (m, x, r0, atimes_param);
  my_daxpy (m, -1.0, b, 1, r0, 1);
  my_dscal (m, -1.0, r0, 1);
  // r = r0
  my_dcopy (m, r0, 1, r, 1);
  // y = -r0
  my_dcopy (m, r0, 1, y, 1);
  my_dscal (m, -1.0, y, 1);
  // rho0 = (r0, r0)
  double rho0 = my_ddot (m, r0, 1, r0, 1);

  int iter;
  for (iter = 0; iter<= kend; iter++)
    {
      // res2 = (r, r)
      res2 = my_ddot (m, r, 1, r, 1);
      if (it->debug == 2)
	{
	  fprintf (it->out, "libiter-gpb(myblas) %d %e\n", iter, res2 / b2);
	}
      if(res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}

      // p = r + beta * (p - u)
      my_daxpy (m, -1.0, u, 1, p, 1); // p = p - u
      my_dscal (m, beta, p, 1);       // p = beta * (p - u)
      my_daxpy (m, 1.0, r, 1, p, 1);  // p = r + beta * (p - u)
      // u = y + beta * u
      my_dscal (m, beta, u, 1);       // u = beta * u
      my_daxpy (m, 1.0, y, 1, u, 1);  // u = y + beta * u

      // q = A.p
      atimes (m, p, q, atimes_param);
      // pap = (r0, q)
      double pap = my_ddot (m, r0, 1, q, 1);

      double alpha = rho0 / pap;

      // y = y + alpha * (q - w)
      my_dcopy (m, q, 1, tmp, 1);        // tmp = q
      my_daxpy (m, -1.0, w, 1, tmp, 1);  // tmp = q - w
      my_daxpy (m, alpha, tmp, 1, y, 1); // y += alpha * (q - w)
      // r = r - alpha * q
      my_daxpy (m, -alpha, q, 1, r, 1);
      // x = x + alpha * p
      my_daxpy (m, alpha, p, 1, x, 1);

      atimes (m, r, w, atimes_param); // w [] -> wn = A tn , tn <- r []
      // AA = (A tn, A tn)
      // yy = (  yn,   yn)
      // Ay = (A tn,   yn)
      // At = (A tn,   tn)
      // yt = (  yn,   tn)
      double AA = my_ddot (m, w, 1, w, 1);
      double yy = my_ddot (m, y, 1, y, 1);
      double Ay = my_ddot (m, w, 1, y, 1);
      double At = my_ddot (m, w, 1, r, 1);
      double yt = my_ddot (m, y, 1, r, 1);

      double eta;
      double zeta;
      // this part is the only difference from Bi-CGSTAB2 (sta2() above)
      if (iter == 0)
	{
	  zeta = At / AA;
	  eta  = 0.0;
        }
      else
	{
	  zeta = (yy * At - Ay * yt) / (AA * yy - Ay * Ay);
	  eta  = (AA * yt - At * Ay) / (AA * yy - Ay * Ay);
	}

      // z = zeta * r + eta * (z - alpha * u)
      my_daxpy (m, -alpha, u, 1, z, 1); // z = z - alpha * u
      my_dscal (m, eta, z, 1);          // z = eta (z - alpha * u)
      my_daxpy (m, zeta, r, 1, z, 1);   // z = zeta * r + eta (z - alpha * u)
      // x = x + z
      my_daxpy (m, 1.0, z, 1, x, 1);

      // y = eta * y + zeta * w
      my_dscal (m, eta, y, 1);        // y = eta * y
      my_daxpy (m, zeta, w, 1, y, 1); // y = eta * y + zeta * w
      // r = r - y
      my_daxpy (m, -1.0, y, 1, r, 1);
      // u = zeta * q + eta * u
      my_dscal (m, eta, u, 1);        // u = eta * u
      my_daxpy (m, zeta, q, 1, u, 1); // u = eta * u + zeta * q
      // rho1 = (r0, r)
      double rho1 = my_ddot (m, r0, 1, r, 1);

      beta = (rho1 / rho0) * (alpha / zeta);
      rho0 = rho1;
      // w = w + beta * q
      my_daxpy (m, beta, q, 1, w, 1);
    }

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H

  free (r0);
  free (w);
  free (q);
  free (u);
  free (z);
  free (y);
  free (r);
  free (p);
  free (tmp);

  if (it->debug == 1)
    {
      fprintf (it->out, "libiter-gpb %d %e\n", iter, res2);
    }

  it->niter = iter;
  it->res2  = res2 / b2;
  return (ret);
}

/* gpbi-cg method with precondition
 * ref: Zhang, SIAM J.Sci.Comput. 1997 vol.18 pp.537-551.
 * INPUT
 *   m : dimension of the problem
 *   b [m] : r-h-s vector
 *   atimes (int m, static double *x, double *b, void *param) :
 *        calc matrix-vector product A.x = b.
 *   atimes_param : parameters for atimes().
 *   inv (int m, static double *b, double *x, void *param) :
 *        approx of A^{-1}.b = x for preconditioning.
 *   inv_param : parameters for the preconditioner inv().
 *   it : struct iter. following entries are used
 *        it->max = kend : max of iteration
 *        it->eps = eps : log10 of cutoff
 * OUTPUT
 *   returned value : 0 == success, otherwise (-1) == failed
 *   x [m] : solution
 *   it->niter : # of iteration
 *   it->res2  : |r^2| / |b^2|
 */
int
gpb_pc (int m, const double *b, double *x,
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
  int kend = it->max;
  double eps2 = it->eps * it->eps;

  double *r0  = (double *)calloc (m, sizeof (double));
  double *w   = (double *)calloc (m, sizeof (double));
  double *q   = (double *)calloc (m, sizeof (double));
  double *u   = (double *)calloc (m, sizeof (double));
  double *z   = (double *)calloc (m, sizeof (double));
  double *y   = (double *)calloc (m, sizeof (double));
  double *r   = (double *)calloc (m, sizeof (double));
  double *p   = (double *)calloc (m, sizeof (double));
  double *tmp = (double *)calloc (m, sizeof (double));
  double *Kr  = (double *)calloc (m, sizeof (double));
  CHECK_MALLOC (r0,  "gpb_pc");
  CHECK_MALLOC (w,   "gpb_pc");
  CHECK_MALLOC (q,   "gpb_pc");
  CHECK_MALLOC (u,   "gpb_pc");
  CHECK_MALLOC (z,   "gpb_pc");
  CHECK_MALLOC (y,   "gpb_pc");
  CHECK_MALLOC (r,   "gpb_pc");
  CHECK_MALLOC (p,   "gpb_pc");
  CHECK_MALLOC (tmp, "gpb_pc");
  CHECK_MALLOC (Kr,  "gpb_pc");

  double res2 = 0.0; // for compiler warning.
  double beta = 0.0;

#ifdef HAVE_CBLAS_H
  /**
   * ATLAS version
   */

  double b2 = cblas_ddot (m, b, 1, b, 1); // (b,b)
  eps2 *= b2;

  // r0 = b - A.x
  atimes (m, x, r0, atimes_param);
  cblas_daxpy (m, -1.0, b, 1, r0, 1);
  cblas_dscal (m, -1.0, r0, 1);
  // r = r0
  cblas_dcopy (m, r0, 1, r, 1);
  // y = -r0
  cblas_dcopy (m, r0, 1, y, 1);
  cblas_dscal (m, -1.0, y, 1);
  // rho0 = (r0, r0)
  double rho0 = cblas_ddot (m, r0, 1, r0, 1);

  int iter;
  for (iter = 0; iter<= kend; iter++)
    {
      // res2 = (r, r)
      res2 = cblas_ddot (m, r, 1, r, 1);
      if (it->debug == 2)
	{
	  fprintf (it->out, "libiter-gpb_pc(atlas) %d %e\n", iter, res2 / b2);
	}
      if(res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}

      // Kr = K^{-1}.r
      inv (m, r, Kr, inv_param);
      // p = K^{-1}.r + beta * (p - u)
      cblas_daxpy (m, -1.0, u, 1, p, 1); // p = p - u
      cblas_dscal (m, beta, p, 1);       // p = beta * (p - u)
      cblas_daxpy (m, 1.0, Kr, 1, p, 1);// p = K^{-1}.r + beta * (p - u)
      // u = K^{-1}.y + beta * u
      inv (m, y, tmp, inv_param);       // tmp = K^{-1}.y
      cblas_dscal (m, beta, u, 1);      // u = beta * u
      cblas_daxpy (m, 1.0, tmp, 1, u, 1); // u = K^{-1}.y + beta * u

      // q = A.p
      atimes (m, p, q, atimes_param);
      // pap = (r0, q)
      double pap = cblas_ddot (m, r0, 1, q, 1);

      double alpha = rho0 / pap;
      // y = y + alpha * (q - w)
      cblas_dcopy (m, q, 1, tmp, 1);        // tmp = q
      cblas_daxpy (m, -1.0, w, 1, tmp, 1);  // tmp = q - w
      cblas_daxpy (m, alpha, tmp, 1, y, 1); // y += alpha * (q - w)
      // r = r - alpha * q, note that at this point this is t_n
      cblas_daxpy (m, -alpha, q, 1, r, 1);
      // x = x + alpha * p
      cblas_daxpy (m, alpha, p, 1, x, 1);

      // w = A.K^{-1}.r (= A.K^{-1}.tn)
      inv (m, r, Kr, inv_param); // Kr = K^{-1}.tn, not equal to Kr above
      atimes (m, Kr, w, atimes_param);// w = A.K^{-1}.tn
      // AA = (A tn, A tn)
      // yy = (  yn,   yn)
      // Ay = (A tn,   yn)
      // At = (A tn,   tn)
      // yt = (  yn,   tn)
      double AA = cblas_ddot (m, w, 1, w, 1);
      double yy = cblas_ddot (m, y, 1, y, 1);
      double Ay = cblas_ddot (m, w, 1, y, 1);
      double At = cblas_ddot (m, w, 1, r, 1);
      double yt = cblas_ddot (m, y, 1, r, 1);

      double eta;
      double zeta;
      // this part is the only difference from Bi-CGSTAB2 (sta2() above)
      if (iter == 0)
	{
	  zeta = At / AA;
	  eta  = 0.0;
        }
      else
	{
	  zeta = (yy * At - Ay * yt) / (AA * yy - Ay * Ay);
	  eta  = (AA * yt - At * Ay) / (AA * yy - Ay * Ay);
	}

      // z = zeta * K^{-1}.r + eta * (z - alpha * u)
      cblas_daxpy (m, -alpha, u, 1, z, 1);// z = z - alpha * u
      cblas_dscal (m, eta, z, 1);         // z = eta (z - alpha * u)
      cblas_daxpy (m, zeta, Kr, 1, z, 1); // z = zeta * Kr + eta (z - alpha * u)
      // x = x + z
      cblas_daxpy (m, 1.0, z, 1, x, 1);

      // y = eta * y + zeta * w
      cblas_dscal (m, eta, y, 1);        // y = eta * y
      cblas_daxpy (m, zeta, w, 1, y, 1); // y = eta * y + zeta * w
      // r = r - y
      cblas_daxpy (m, -1.0, y, 1, r, 1);
      // u = zeta * q + eta * u
      inv (m, q, tmp, inv_param);        // tmp = K^{-1}.q
      cblas_dscal (m, eta, u, 1);        // u = eta * u
      cblas_daxpy (m, zeta, tmp, 1, u, 1); // u = eta * u + zeta * K^{-1}.q
      // rho1 = (r0, r)
      double rho1 = cblas_ddot (m, r0, 1, r, 1);

      beta = (rho1 / rho0) * (alpha / zeta);
      rho0 = rho1;
      // w = w + beta * q
      cblas_daxpy (m, beta, q, 1, w, 1);
    }

#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
  /**
   * BLAS version
   */

  double b2 = ddot_ (&m, b, &i_1, b, &i_1); // (b,b)
  eps2 *= b2;

  // r0 = b - A.x
  atimes (m, x, r0, atimes_param);
  daxpy_ (&m, &d_m1, b, &i_1, r0, &i_1);
  dscal_ (&m, &d_m1, r0, &i_1);
  // r = r0
  dcopy_ (&m, r0, &i_1, r, &i_1);
  // y = -r0
  dcopy_ (&m, r0, &i_1, y, &i_1);
  dscal_ (&m, &d_m1, y, &i_1);
  // rho0 = (r0, r0)
  double rho0 = ddot_ (&m, r0, &i_1, r0, &i_1);

  int iter;
  for (iter = 0; iter<= kend; iter++)
    {
      // res2 = (r, r)
      res2 = ddot_ (&m, r, &i_1, r, &i_1);
      if (it->debug == 2)
	{
	  fprintf (it->out, "libiter-gpb_pc(blas) %d %e\n", iter, res2 / b2);
	}
      if(res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}

      // Kr = K^{-1}.r
      inv (m, r, Kr, inv_param);
      // p = K^{-1}.r + beta * (p - u)
      daxpy_ (&m, &d_m1, u, &i_1, p, &i_1); // p = p - u
      dscal_ (&m, &beta, p, &i_1);         // p = beta * (p - u)
      daxpy_ (&m, &d_1, Kr, &i_1, p, &i_1);// p = K^{-1}.r + beta * (p - u)
      // u = K^{-1}.y + beta * u
      inv (m, y, tmp, inv_param);            // tmp = K^{-1}.y
      dscal_ (&m, &beta, u, &i_1);           // u = beta * u
      daxpy_ (&m, &d_1, tmp, &i_1, u, &i_1); // u = K^{-1}.y + beta * u

      // q = A.p
      atimes (m, p, q, atimes_param);
      // pap = (r0, q)
      double pap = ddot_ (&m, r0, &i_1, q, &i_1);

      double alpha = rho0 / pap;
      // y = y + alpha * (q - w)
      dcopy_ (&m, q, &i_1, tmp, &i_1);        // tmp = q
      daxpy_ (&m, &d_m1, w, &i_1, tmp, &i_1);  // tmp = q - w
      daxpy_ (&m, &alpha, tmp, &i_1, y, &i_1); // y += alpha * (q - w)
      // r = r - alpha * q, note that at this point this is t_n
      double malpha = -alpha;
      daxpy_ (&m, &malpha, q, &i_1, r, &i_1);
      // x = x + alpha * p
      daxpy_ (&m, &alpha, p, &i_1, x, &i_1);

      // w = A.K^{-1}.r (= A.K^{-1}.tn)
      inv (m, r, Kr, inv_param); // Kr = K^{-1}.tn, not equal to Kr above
      atimes (m, Kr, w, atimes_param);// w = A.K^{-1}.tn
      // AA = (A tn, A tn)
      // yy = (  yn,   yn)
      // Ay = (A tn,   yn)
      // At = (A tn,   tn)
      // yt = (  yn,   tn)
      double AA = ddot_ (&m, w, &i_1, w, &i_1);
      double yy = ddot_ (&m, y, &i_1, y, &i_1);
      double Ay = ddot_ (&m, w, &i_1, y, &i_1);
      double At = ddot_ (&m, w, &i_1, r, &i_1);
      double yt = ddot_ (&m, y, &i_1, r, &i_1);

      double eta;
      double zeta;
      // this part is the only difference from Bi-CGSTAB2 (sta2() above)
      if (iter == 0)
	{
	  zeta = At / AA;
	  eta  = 0.0;
        }
      else
	{
	  zeta = (yy * At - Ay * yt) / (AA * yy - Ay * Ay);
	  eta  = (AA * yt - At * Ay) / (AA * yy - Ay * Ay);
	}

      // z = zeta * K^{-1}.r + eta * (z - alpha * u)
      daxpy_ (&m, &malpha, u, &i_1, z, &i_1);// z = z - alpha * u
      dscal_ (&m, &eta, z, &i_1);           // z = eta (z - alpha * u)
      daxpy_ (&m, &zeta, Kr, &i_1, z, &i_1);// z = zeta*Kr + eta(z - alpha*u)
      // x = x + z
      daxpy_ (&m, &d_1, z, &i_1, x, &i_1);

      // y = eta * y + zeta * w
      dscal_ (&m, &eta, y, &i_1);           // y = eta * y
      daxpy_ (&m, &zeta, w, &i_1, y, &i_1); // y = eta * y + zeta * w
      // r = r - y
      daxpy_ (&m, &d_m1, y, &i_1, r, &i_1);
      // u = zeta * q + eta * u
      inv (m, q, tmp, inv_param);        // tmp = K^{-1}.q
      dscal_ (&m, &eta, u, &i_1);        // u = eta * u
      daxpy_ (&m, &zeta, tmp, &i_1, u, &i_1); // u = eta * u + zeta * K^{-1}.q
      // rho1 = (r0, r)
      double rho1 = ddot_ (&m, r0, &i_1, r, &i_1);

      beta = (rho1 / rho0) * (alpha / zeta);
      rho0 = rho1;
      // w = w + beta * q
      daxpy_ (&m, &beta, q, &i_1, w, &i_1);
    }

# else // !HAVE_BLAS_H
  /**
   * local BLAS version
   */

  double b2 = my_ddot (m, b, 1, b, 1); // (b,b)
  eps2 *= b2;

  // r0 = b - A.x
  atimes (m, x, r0, atimes_param);
  my_daxpy (m, -1.0, b, 1, r0, 1);
  my_dscal (m, -1.0, r0, 1);
  // r = r0
  my_dcopy (m, r0, 1, r, 1);
  // y = -r0
  my_dcopy (m, r0, 1, y, 1);
  my_dscal (m, -1.0, y, 1);
  // rho0 = (r0, r0)
  double rho0 = my_ddot (m, r0, 1, r0, 1);

  int iter;
  for (iter = 0; iter<= kend; iter++)
    {
      // res2 = (r, r)
      res2 = my_ddot (m, r, 1, r, 1);
      if (it->debug == 2)
	{
	  fprintf (it->out, "libiter-gpb_pc(myblas) %d %e\n", iter, res2 / b2);
	}
      if(res2 <= eps2)
	{
	  ret = 0; // success
	  break;
	}

      // Kr = K^{-1}.r
      inv (m, r, Kr, inv_param);
      // p = K^{-1}.r + beta * (p - u)
      my_daxpy (m, -1.0, u, 1, p, 1); // p = p - u
      my_dscal (m, beta, p, 1);       // p = beta * (p - u)
      my_daxpy (m, 1.0, Kr, 1, p, 1);// p = K^{-1}.r + beta * (p - u)
      // u = K^{-1}.y + beta * u
      inv (m, y, tmp, inv_param);       // tmp = K^{-1}.y
      my_dscal (m, beta, u, 1);      // u = beta * u
      my_daxpy (m, 1.0, tmp, 1, u, 1); // u = K^{-1}.y + beta * u

      // q = A.p
      atimes (m, p, q, atimes_param);
      // pap = (r0, q)
      double pap = my_ddot (m, r0, 1, q, 1);

      double alpha = rho0 / pap;
      // y = y + alpha * (q - w)
      my_dcopy (m, q, 1, tmp, 1);        // tmp = q
      my_daxpy (m, -1.0, w, 1, tmp, 1);  // tmp = q - w
      my_daxpy (m, alpha, tmp, 1, y, 1); // y += alpha * (q - w)
      // r = r - alpha * q, note that at this point this is t_n
      my_daxpy (m, -alpha, q, 1, r, 1);
      // x = x + alpha * p
      my_daxpy (m, alpha, p, 1, x, 1);

      // w = A.K^{-1}.r (= A.K^{-1}.tn)
      inv (m, r, Kr, inv_param); // Kr = K^{-1}.tn, not equal to Kr above
      atimes (m, Kr, w, atimes_param);// w = A.K^{-1}.tn
      // AA = (A tn, A tn)
      // yy = (  yn,   yn)
      // Ay = (A tn,   yn)
      // At = (A tn,   tn)
      // yt = (  yn,   tn)
      double AA = my_ddot (m, w, 1, w, 1);
      double yy = my_ddot (m, y, 1, y, 1);
      double Ay = my_ddot (m, w, 1, y, 1);
      double At = my_ddot (m, w, 1, r, 1);
      double yt = my_ddot (m, y, 1, r, 1);

      double eta;
      double zeta;
      // this part is the only difference from Bi-CGSTAB2 (sta2() above)
      if (iter == 0)
	{
	  zeta = At / AA;
	  eta  = 0.0;
        }
      else
	{
	  zeta = (yy * At - Ay * yt) / (AA * yy - Ay * Ay);
	  eta  = (AA * yt - At * Ay) / (AA * yy - Ay * Ay);
	}

      // z = zeta * K^{-1}.r + eta * (z - alpha * u)
      my_daxpy (m, -alpha, u, 1, z, 1);// z = z - alpha * u
      my_dscal (m, eta, z, 1);         // z = eta (z - alpha * u)
      my_daxpy (m, zeta, Kr, 1, z, 1); // z = zeta * Kr + eta (z - alpha * u)
      // x = x + z
      my_daxpy (m, 1.0, z, 1, x, 1);

      // y = eta * y + zeta * w
      my_dscal (m, eta, y, 1);        // y = eta * y
      my_daxpy (m, zeta, w, 1, y, 1); // y = eta * y + zeta * w
      // r = r - y
      my_daxpy (m, -1.0, y, 1, r, 1);
      // u = zeta * q + eta * u
      inv (m, q, tmp, inv_param);        // tmp = K^{-1}.q
      my_dscal (m, eta, u, 1);        // u = eta * u
      my_daxpy (m, zeta, tmp, 1, u, 1); // u = eta * u + zeta * K^{-1}.q
      // rho1 = (r0, r)
      double rho1 = my_ddot (m, r0, 1, r, 1);

      beta = (rho1 / rho0) * (alpha / zeta);
      rho0 = rho1;
      // w = w + beta * q
      my_daxpy (m, beta, q, 1, w, 1);
    }

# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H

  free (r0);
  free (w);
  free (q);
  free (u);
  free (z);
  free (y);
  free (r);
  free (p);
  free (tmp);
  free (Kr);

  if (it->debug == 1)
    {
      fprintf (it->out, "libiter-gpb_pc %d %e\n", iter, res2);
    }

  it->niter = iter;
  it->res2  = res2 / b2;
  return (ret);
}
