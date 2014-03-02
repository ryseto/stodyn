/* generalized minimum residual method
 * Copyright (C) 1998-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: gmres.c,v 2.14 2007/12/08 20:37:17 kichiki Exp $
 *
 * Reference :
 *   GMRES(m) : Y.Saad & M.H.Schultz, SIAM J.Sci.Stat.Comput.
 *   vol7 (1986) pp.856-869
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

#include <cstdlib> /* malloc (), free() */
#include <cstdio>
#include <cmath>
#include "libiter.h"
#include "memory-check.h" // CHECK_MALLOC

#ifdef HAVE_CBLAS_H
/* use ATLAS' CBLAS routines */


#ifdef OSX
#include <Accelerate/Accelerate.h>
#else
extern "C"{
#include <atlas/cblas.h>
	//#include <atlas/clapack.h>
	
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


#include "gmres.h"


/*
 * m  : number of iteration
 * nn : dimension of matrix r [] (nn * nn)
 */
static void
back_sub (int m, int nn,
		  const double *r, const double *g,
		  double *y)
{
	int i, j, jj;
	
	/*for (j = m - 1;j >= 0;j --){*/
	/* above for-loop fail, because j is unsigned!! */
	for (j = m - 1, jj = 0; jj < m; j --, jj ++)
    {
		y [j] = g [j];
		for (i = j + 1; i < m; i ++)
		{
			y [j] -= r [j * nn + i] * y [i];
		}
		y [j] /= r [j * nn + j];
    }
}

/* GMRES(m)
 * INPUT
 *   n : dimension of the problem
 *   b[n] : r-h-s vector
 *   atimes (int n, static double *x, double *b, void *param) :
 *        calc matrix-vector product A.x = b.
 *   atimes_param : parameters for atimes().
 *   it : struct iter. the following entries are used.
 *        it->max     : max of iteration
 *        it->restart : number of iteration at once
 *        it->eps     : criteria for |r^2|/|b^2|
 * OUTPUT
 *   returned value : 0 == success, otherwise (-1) == failed
 *   x[m] : solution
 *   it->niter : # of iteration
 *   it->res2  : |r^2| / |b^2|
 */
int
gmres_m (int n, const double *f, double *x,
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
	double scale;
	
# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H
	
	int ret = -1;
	int m = it->restart;
	int itmax = it->max;
	double eps = it->eps;
	
	double *v = (double *)malloc (sizeof (double) * (m + 1) * n);
	double *h = (double *)malloc (sizeof (double) * m * m);
	double *g = (double *)malloc (sizeof (double) * (m + 1));
	double *c = (double *)malloc (sizeof (double) * m);
	double *s = (double *)malloc (sizeof (double) * m);
	double *yv= (double *)malloc (sizeof (double) * n);
	CHECK_MALLOC (v, "gmres_m");
	CHECK_MALLOC (h, "gmres_m");
	CHECK_MALLOC (g, "gmres_m");
	CHECK_MALLOC (c, "gmres_m");
	CHECK_MALLOC (s, "gmres_m");
	CHECK_MALLOC (yv,"gmres_m");
	
	double res = 0.0;
	
	int i, j, k;
	double hv;
	double rr, hh;
	double r1, r2;
	double g0;
	
#ifdef HAVE_CBLAS_H
	/**
	 * ATLAS version
	 */
	
	double b2 = cblas_ddot (n, f, 1, f, 1); // (f,f)
	eps *= sqrt (b2);
	
	int iter = 0;
	/* 1. start: */
	/* compute r0 */
	/* compute v1 */
	/* beta */
	atimes (n, x, v + 0, atimes_param); /* use v [0] temporaliry */
	
	// v = f - v
	cblas_dscal (n, -1.0, v, 1); // v = - v
	cblas_daxpy (n, 1.0, f, 1, v, 1); // v = f - v
	
	//g[0] = cblas_dnrm2 (n, v, 1);
	g[0] = sqrt (cblas_ddot (n, v, 1, v, 1));
	cblas_dscal (n, 1.0 / g[0], v, 1);
	
	/* main loop */
	while (iter <= itmax)
    {
		++iter;
		/* 2. iterate: */
		for (j = 0; j < m; j ++)
		{
			/* tmp = A.vj : use v [(j + 1) * n] directly */
			atimes (n, v + j * n, v + (j + 1) * n, atimes_param);
			/* h_i,j (i=1,...,j) */
			for (i = 0; i <= j; i ++)
			{
				h [i * m + j] =
				cblas_ddot (n, v + (j + 1) * n, 1,
							v + i * n, 1);
			}
			/* vv_j+1 */
			for (k = 0; k < n; k ++)
			{
				hv = 0.0;
				for (i = 0; i <= j; i ++)
				{
					hv += h [i * m + j] * v [i * n + k];
				}
				v [(j + 1) * n + k] -= hv;
			}
			/* h_j+1,j */
			/* v_j+1 */
			//hh = cblas_dnrm2 (n, v + (j + 1) * n, 1);
			hh = sqrt (cblas_ddot (n, v + (j + 1) * n, 1, v + (j + 1) * n, 1));
			cblas_dscal (n, 1.0 / hh, v + (j + 1) * n, 1);
			
			/* rotate */
			for (i = 0; i < j; i ++)
			{
				r1 = h [ i      * m + j];
				r2 = h [(i + 1) * m + j];
				h [ i      * m + j] = c [i] * r1 - s [i] * r2;
				h [(i + 1) * m + j] = s [i] * r1 + c [i] * r2;
			}
			rr = h [j * m + j];
			hv = sqrt (rr * rr + hh * hh); /* temporary variable */
			c [j] =  rr / hv;
			s [j] = -hh / hv;
			h [j * m + j] = hv; /* resultant (after rotated) element */
			
			g0 = g [j];
			g [j    ] = c [j] * g0;
			g [j + 1] = s [j] * g0;
		}
		/* 3. form the approximate solution */
		/* solve y_k */
		back_sub (j, m, h, g, c); /* use c[m] for y_k */
		for (i = 0; i < n; i ++)
		{
			yv[i] = 0.0;
		}
		for (k = 0; k < j; k ++)
		{
			cblas_daxpy (n, c[k], v + k*n, 1, yv, 1);
		}
		// x_m = x_0 + (y_k v_k)
		cblas_daxpy (n, 1.0, yv, 1, x, 1);
		
		/* 4. restart */
		res = fabs (g [j/*m*/]); /* residual */
		/*fprintf (stderr, "# iter %d res %e\n", iter, *res);*/
		/* if satisfied, */
		if (it->debug == 2)
		{
			fprintf (it->out, "libiter-gmres(%d) %d %d %e\n",
					 m, iter, j, res*res);
		}
		if (res <= eps)
		{
			ret = 0; // success
			break;
		}
		/* else */
		/* compute r_m */
		atimes (n, x, v + 0, atimes_param);
		/* r_m */
		/* compute v1 */
		
		// v = f - v
		cblas_dscal (n, -1.0, v, 1); // v = - v
		cblas_daxpy (n, 1.0, f, 1, v, 1); // v = f - v
		
		//g [0] = cblas_dnrm2 (n, v, 1);
		g[0] = sqrt (cblas_ddot (n, v, 1, v, 1));
		cblas_dscal (n, 1.0 / g[0], v, 1);
    }
	
#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
	/**
	 * BLAS version
	 */
	
	double b2 = ddot_ (&n, f, &i_1, f, &i_1); // (f,f)
	eps *= sqrt (b2);
	
	int iter = 0;
	/* 1. start: */
	/* compute r0 */
	/* compute v1 */
	/* beta */
	atimes (n, x, v + 0, atimes_param); /* use v [0] temporaliry */
	
	// v = f - v
	dscal_ (&n, &d_m1, v, &i_1); // v = - v
	daxpy_ (&n, &d_1, f, &i_1, v, &i_1); // v = f - v
	
	g[0] = sqrt (ddot_ (&n, v, &i_1, v, &i_1));
	scale = 1.0 / g[0];
	dscal_ (&n, &scale, v, &i_1);
	
	/* main loop */
	while (iter <= itmax)
    {
		++iter;
		/* 2. iterate: */
		for (j = 0; j < m; j ++)
		{
			/* tmp = A.vj : use v [(j + 1) * n] directly */
			atimes (n, v + j * n, v + (j + 1) * n, atimes_param);
			/* h_i,j (i=1,...,j) */
			for (i = 0; i <= j; i ++)
			{
				h [i * m + j] =
				ddot_ (&n, v + (j + 1) * n, &i_1,
					   v + i * n, &i_1);
			}
			/* vv_j+1 */
			for (k = 0; k < n; k ++)
			{
				hv = 0.0;
				for (i = 0; i <= j; i ++)
				{
					hv += h [i * m + j] * v [i * n + k];
				}
				v [(j + 1) * n + k] -= hv;
			}
			/* h_j+1,j */
			/* v_j+1 */
			//hh = dnrm2_ (&n, v + (j + 1) * n, &i_1);
			hh = sqrt (ddot_ (&n, v + (j + 1) * n, &i_1, v + (j + 1) * n, &i_1));
			scale = 1.0 / hh;
			dscal_ (&n, &scale, v + (j + 1) * n, &i_1);
			
			/* rotate */
			for (i = 0; i < j; i ++)
			{
				r1 = h [ i      * m + j];
				r2 = h [(i + 1) * m + j];
				h [ i      * m + j] = c [i] * r1 - s [i] * r2;
				h [(i + 1) * m + j] = s [i] * r1 + c [i] * r2;
			}
			rr = h [j * m + j];
			hv = sqrt (rr * rr + hh * hh); /* temporary variable */
			c [j] =  rr / hv;
			s [j] = -hh / hv;
			h [j * m + j] = hv; /* resultant (after rotated) element */
			
			g0 = g [j];
			g [j    ] = c [j] * g0;
			g [j + 1] = s [j] * g0;
		}
		/* 3. form the approximate solution */
		/* solve y_k */
		back_sub (j, m, h, g, c); /* use c[m] for y_k */
		for (i = 0; i < n; i ++)
		{
			yv[i] = 0.0;
		}
		for (k = 0; k < j; k ++)
		{
			daxpy_ (&n, c + k, v + k*n, &i_1, yv, &i_1);
		}
		// x_m = x_0 + (y_k v_k)
		daxpy_ (&n, &d_1, yv, &i_1, x, &i_1);
		
		/* 4. restart */
		res = fabs (g [j/*m*/]); /* residual */
		/*fprintf (stderr, "# iter %d res %e\n", iter, *res);*/
		/* if satisfied, */
		if (it->debug == 2)
		{
			fprintf (it->out, "libiter-gmres(%d) %d %d %e\n",
					 m, iter, j, res*res);
		}
		if (res <= eps)
		{
			ret = 0; // success
			break;
		}
		/* else */
		/* compute r_m */
		atimes (n, x, v + 0, atimes_param);
		/* r_m */
		/* compute v1 */
		
		// v = f - v
		dscal_ (&n, &d_m1, v, &i_1); // v = - v
		daxpy_ (&n, &d_1, f, &i_1, v, &i_1); // v = f - v
		
		g[0] = sqrt (ddot_ (&n, v, &i_1, v, &i_1));
		scale = 1.0 / g[0];
		dscal_ (&n, &scale, v, &i_1);
    }
	
# else // !HAVE_BLAS_H
	/**
	 * local BLAS version
	 */
	
	double b2 = my_ddot (n, f, 1, f, 1); // (f,f)
	eps *= sqrt (b2);
	
	int iter = 0;
	/* 1. start: */
	/* compute r0 */
	/* compute v1 */
	/* beta */
	atimes (n, x, v + 0, atimes_param); /* use v [0] temporaliry */
	
	my_daxpyz (n, -1.0, v + 0, 1, f, 1, v + 0, 1);
	
	g [0] = my_dnrm2 (n, v + 0, 1);
	my_dscal (n, 1.0 / g [0], v + 0, 1);
	
	/* main loop */
	while (iter <= itmax)
    {
		++iter;
		/* 2. iterate: */
		for (j = 0; j < m; j ++)
		{
			/* tmp = A.vj : use v [(j + 1) * n] directly */
			atimes (n, v + j * n, v + (j + 1) * n, atimes_param);
			/* h_i,j (i=1,...,j) */
			for (i = 0; i <= j; i ++)
			{
				h [i * m + j] =
				my_ddot (n, v + (j + 1) * n, 1,
						 v + i * n, 1);
			}
			/* vv_j+1 */
			for (k = 0; k < n; k ++)
			{
				hv = 0.0;
				for (i = 0; i <= j; i ++)
				{
					hv += h [i * m + j] * v [i * n + k];
				}
				v [(j + 1) * n + k] -= hv;
			}
			/* h_j+1,j */
			/* v_j+1 */
			hh = my_dnrm2 (n, v + (j + 1) * n, 1);
			my_dscal (n, 1.0 / hh, v + (j + 1) * n, 1);
			
			/* rotate */
			for (i = 0; i < j; i ++)
			{
				r1 = h [ i      * m + j];
				r2 = h [(i + 1) * m + j];
				h [ i      * m + j] = c [i] * r1 - s [i] * r2;
				h [(i + 1) * m + j] = s [i] * r1 + c [i] * r2;
			}
			rr = h [j * m + j];
			hv = sqrt (rr * rr + hh * hh); /* temporary variable */
			c [j] =  rr / hv;
			s [j] = -hh / hv;
			h [j * m + j] = hv; /* resultant (after rotated) element */
			
			g0 = g [j];
			g [j    ] = c [j] * g0;
			g [j + 1] = s [j] * g0;
		}
		/* 3. form the approximate solution */
		/* solve y_k */
		back_sub (j, m, h, g, c); /* use c[m] for y_k */
		for (i = 0; i < n; i ++)
		{
			yv[i] = 0.0;
		}
		for (k = 0; k < j; k ++)
		{
			my_daxpy (n, c[k], v + k*n, 1, yv, 1);
		}
		// x_m = x_0 + (y_k v_k)
		my_daxpy (n, 1.0, yv, 1, x, 1);
		
		/* 4. restart */
		res = fabs (g [j/*m*/]); /* residual */
		/*fprintf (stderr, "# iter %d res %e\n", iter, *res);*/
		/* if satisfied, */
		if (it->debug == 2)
		{
			fprintf (it->out, "libiter-gmres(%d) %d %d %e\n",
					 m, iter, j, res*res);
		}
		if (res <= eps)
		{
			ret = 0; // success
			break;
		}
		/* else */
		/* compute r_m */
		atimes (n, x, v + 0, atimes_param);
		/* r_m */
		/* compute v1 */
		
		my_daxpyz (n, -1.0, v + 0, 1, f, 1, v + 0, 1);
		
		g [0] = my_dnrm2 (n, v + 0, 1);
		my_dscal (n, 1.0 / g [0], v + 0, 1);
    }
	
# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H
	
	free (v);
	free (h);
	free (g);
	free (c);
	free (s);
	free (yv);
	
	/* adjust iter */
	iter *= m;
	
	if (it->debug == 1)
    {
		fprintf (it->out, "libiter-gmres(%d) %d %e\n",
				 m, iter, res*res / b2);
    }
	
	it->niter = iter;
	it->res2  = res * res / b2;
	return (ret);
}

/* GMRES(m) with preconditioning
 * INPUT
 *   n : dimension of the problem
 *   b[n] : r-h-s vector
 *   atimes (int n, static double *x, double *b, void *param) :
 *        calc matrix-vector product A.x = b.
 *   atimes_param : parameters for atimes().
 *   inv (int n, static double *b, double *x, void *param) :
 *        approx of A^{-1}.b = x for preconditioning.
 *   inv_param : parameters for the preconditioner inv().
 *   it : struct iter. the following entries are used.
 *        it->max     : max of iteration
 *        it->restart : number of iteration at once
 *        it->eps     : criteria for |r^2|/|b^2|
 * OUTPUT
 *   returned value : 0 == success, otherwise (-1) == failed
 *   x[n] : solution
 *   it->niter : # of iteration
 *   it->res2  : |r^2| / |b^2|
 */
int
gmres_m_pc (int n, const double *f, double *x,
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
	double scale;
	
# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H
	
	int ret = -1;
	int m = it->restart;
	int itmax = it->max;
	double eps = it->eps;
	
	double *v = (double *)malloc (sizeof (double) * (m + 1) * n);
	double *h = (double *)malloc (sizeof (double) * m * m);
	double *g = (double *)malloc (sizeof (double) * (m + 1));
	double *c = (double *)malloc (sizeof (double) * m);
	double *s = (double *)malloc (sizeof (double) * m);
	double *yv= (double *)malloc (sizeof (double) * n);
	double *Kv= (double *)malloc (sizeof (double) * n);
	CHECK_MALLOC (v, "gmres_m_pc");
	CHECK_MALLOC (h, "gmres_m_pc");
	CHECK_MALLOC (g, "gmres_m_pc");
	CHECK_MALLOC (c, "gmres_m_pc");
	CHECK_MALLOC (s, "gmres_m_pc");
	CHECK_MALLOC (yv,"gmres_m_pc");
	CHECK_MALLOC (Kv,"gmres_m_pc");
	
	double res = 0.0;
	
	int i, j, k;
	double hv;
	double rr, hh;
	double r1, r2;
	double g0;
	
#ifdef HAVE_CBLAS_H
	/**
	 * ATLAS version
	 */
	
	double b2 = cblas_ddot (n, f, 1, f, 1); // (f,f)
	eps *= sqrt (b2);
	
	int iter = 0;
	/* 1. start: */
	/* compute r0 */
	/* compute v1 */
	/* beta */
	atimes (n, x, v + 0, atimes_param); /* use v [0] temporaliry */
	
	// v = f - v
	cblas_dscal (n, -1.0, v, 1); // v = - v
	cblas_daxpy (n, 1.0, f, 1, v, 1); // v = f - v
	
	//g[0] = cblas_dnrm2 (n, v, 1);
	g[0] = sqrt (cblas_ddot (n, v, 1, v, 1));
	cblas_dscal (n, 1.0 / g[0], v, 1);
	
	/* main loop */
	while (iter <= itmax)
    {
		++iter;
		/* 2. iterate: */
		for (j = 0; j < m; j ++)
		{
			/* tmp = A.K^{-1}.vj : use v [(j + 1) * n] directly */
			// Kv = K^{-1}.vj
			inv (n, v + j * n, Kv, inv_param);
			atimes (n, Kv, v + (j + 1) * n, atimes_param);
			/* h_i,j (i=1,...,j) */
			for (i = 0; i <= j; i ++)
			{
				h [i * m + j] =
				cblas_ddot (n, v + (j + 1) * n, 1,
							v + i * n, 1);
			}
			/* vv_j+1 */
			for (k = 0; k < n; k ++)
			{
				hv = 0.0;
				for (i = 0; i <= j; i ++)
				{
					hv += h [i * m + j] * v [i * n + k];
				}
				v [(j + 1) * n + k] -= hv;
			}
			/* h_j+1,j */
			/* v_j+1 */
			//hh = cblas_dnrm2 (n, v + (j + 1) * n, 1);
			hh = sqrt (cblas_ddot (n, v + (j + 1) * n, 1, v + (j + 1) * n, 1));
			cblas_dscal (n, 1.0 / hh, v + (j + 1) * n, 1);
			
			/* rotate */
			for (i = 0; i < j; i ++)
			{
				r1 = h [ i      * m + j];
				r2 = h [(i + 1) * m + j];
				h [ i      * m + j] = c [i] * r1 - s [i] * r2;
				h [(i + 1) * m + j] = s [i] * r1 + c [i] * r2;
			}
			rr = h [j * m + j];
			hv = sqrt (rr * rr + hh * hh); /* temporary variable */
			c [j] =  rr / hv;
			s [j] = -hh / hv;
			h [j * m + j] = hv; /* resultant (after rotated) element */
			
			g0 = g [j];
			g [j    ] = c [j] * g0;
			g [j + 1] = s [j] * g0;
		}
		/* 3. form the approximate solution */
		/* solve y_k */
		back_sub (j, m, h, g, c); /* use c[m] for y_k */
		// yv[n] = y_k v_k[]
		for (i = 0; i < n; i ++)
		{
			yv[i] = 0.0;
		}
		for (k = 0; k < j; k ++)
		{
			cblas_daxpy (n, c[k], v + k*n, 1, yv, 1);
		}
		// Kv[] = K^{-1}.(y_k v_k)
		inv (n, yv, Kv, inv_param);
		// x_m = x_0 + K^{-1}.(y_k v_k)
		cblas_daxpy (n, 1.0, Kv, 1, x, 1);
		
		/* 4. restart */
		res = fabs (g [j]); /* residual */
		/* if satisfied, */
		if (it->debug == 2)
		{
			fprintf (it->out, "libiter-gmres(%d) %d %d %e\n",
					 m, iter, j, res*res);
		}
		if (res <= eps)
		{
			ret = 0; // success
			break;
		}
		/* else */
		/* compute r_m */
		atimes (n, x, v + 0, atimes_param);
		/* r_m */
		/* compute v1 */
		
		// v = f - v
		cblas_dscal (n, -1.0, v, 1); // v = - v
		cblas_daxpy (n, 1.0, f, 1, v, 1); // v = f - v
		
		//g [0] = cblas_dnrm2 (n, v, 1);
		g[0] = sqrt (cblas_ddot (n, v, 1, v, 1));
		cblas_dscal (n, 1.0 / g[0], v, 1);
		
    }
	
#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
	/**
	 * BLAS version
	 */
	
	double b2 = ddot_ (&n, f, &i_1, f, &i_1); // (f,f)
	eps *= sqrt (b2);
	
	int iter = 0;
	/* 1. start: */
	/* compute r0 */
	/* compute v1 */
	/* beta */
	atimes (n, x, v + 0, atimes_param); /* use v [0] temporaliry */
	
	// v = f - v
	dscal_ (&n, &d_m1, v, &i_1); // v = - v
	daxpy_ (&n, &d_1, f, &i_1, v, &i_1); // v = f - v
	
	g[0] = sqrt (ddot_ (&n, v, &i_1, v, &i_1));
	scale = 1.0 / g[0];
	dscal_ (&n, &scale, v, &i_1);
	
	/* main loop */
	while (iter <= itmax)
    {
		++iter;
		/* 2. iterate: */
		for (j = 0; j < m; j ++)
		{
			/* tmp = A.K^{-1}.vj : use v [(j + 1) * n] directly */
			// Kv = K^{-1}.vj
			inv (n, v + j * n, Kv, inv_param);
			atimes (n, Kv, v + (j + 1) * n, atimes_param);
			/* h_i,j (i=1,...,j) */
			for (i = 0; i <= j; i ++)
			{
				h [i * m + j] =
				ddot_ (&n, v + (j + 1) * n, &i_1,
					   v + i * n, &i_1);
			}
			/* vv_j+1 */
			for (k = 0; k < n; k ++)
			{
				hv = 0.0;
				for (i = 0; i <= j; i ++)
				{
					hv += h [i * m + j] * v [i * n + k];
				}
				v [(j + 1) * n + k] -= hv;
			}
			/* h_j+1,j */
			/* v_j+1 */
			//hh = dnrm2_ (&n, v + (j + 1) * n, &i_1);
			hh = sqrt (ddot_ (&n, v + (j + 1) * n, &i_1, v + (j + 1) * n, &i_1));
			scale = 1.0 / hh;
			dscal_ (&n, &scale, v + (j + 1) * n, &i_1);
			
			/* rotate */
			for (i = 0; i < j; i ++)
			{
				r1 = h [ i      * m + j];
				r2 = h [(i + 1) * m + j];
				h [ i      * m + j] = c [i] * r1 - s [i] * r2;
				h [(i + 1) * m + j] = s [i] * r1 + c [i] * r2;
			}
			rr = h [j * m + j];
			hv = sqrt (rr * rr + hh * hh); /* temporary variable */
			c [j] =  rr / hv;
			s [j] = -hh / hv;
			h [j * m + j] = hv; /* resultant (after rotated) element */
			
			g0 = g [j];
			g [j    ] = c [j] * g0;
			g [j + 1] = s [j] * g0;
		}
		/* 3. form the approximate solution */
		/* solve y_k */
		back_sub (j, m, h, g, c); /* use c[m] for y_k */
		// yv[n] = y_k v_k[]
		for (i = 0; i < n; i ++)
		{
			yv[i] = 0.0;
		}
		for (k = 0; k < j; k ++)
		{
			daxpy_ (&n, c + k, v + k*n, &i_1, yv, &i_1);
		}
		// Kv[] = K^{-1}.(y_k v_k)
		inv (n, yv, Kv, inv_param);
		// x_m = x_0 + K^{-1}.(y_k v_k)
		daxpy_ (&n, &d_1, Kv, &i_1, x, &i_1);
		
		/* 4. restart */
		res = fabs (g [j]); /* residual */
		/* if satisfied, */
		if (it->debug == 2)
		{
			fprintf (it->out, "libiter-gmres(%d) %d %d %e\n",
					 m, iter, j, res*res);
		}
		if (res <= eps)
		{
			ret = 0; // success
			break;
		}
		/* else */
		/* compute r_m */
		atimes (n, x, v + 0, atimes_param);
		/* r_m */
		/* compute v1 */
		
		// v = f - v
		dscal_ (&n, &d_m1, v, &i_1); // v = - v
		daxpy_ (&n, &d_1, f, &i_1, v, &i_1); // v = f - v
		
		g[0] = sqrt (ddot_ (&n, v, &i_1, v, &i_1));
		scale = 1.0 / g[0];
		dscal_ (&n, &scale, v, &i_1);
    }
	
# else // !HAVE_BLAS_H
	/**
	 * local BLAS version
	 */
	
	double b2 = my_ddot (n, f, 1, f, 1); // (f,f)
	eps *= sqrt (b2);
	
	int iter = 0;
	/* 1. start: */
	/* compute r0 */
	/* compute v1 */
	/* beta */
	atimes (n, x, v + 0, atimes_param); /* use v [0] temporaliry */
	
	my_daxpyz (n, -1.0, v + 0, 1, f, 1, v + 0, 1);
	
	g [0] = my_dnrm2 (n, v + 0, 1);
	my_dscal (n, 1.0 / g [0], v + 0, 1);
	
	/* main loop */
	while (iter <= itmax)
    {
		++iter;
		/* 2. iterate: */
		for (j = 0; j < m; j ++)
		{
			/* tmp = A.K^{-1}.vj : use v [(j + 1) * n] directly */
			// Kv = K^{-1}.vj
			inv (n, v + j * n, Kv, inv_param);
			atimes (n, Kv, v + (j + 1) * n, atimes_param);
			/* h_i,j (i=1,...,j) */
			for (i = 0; i <= j; i ++)
			{
				h [i * m + j] =
				my_ddot (n, v + (j + 1) * n, 1,
						 v + i * n, 1);
			}
			/* vv_j+1 */
			for (k = 0; k < n; k ++)
			{
				hv = 0.0;
				for (i = 0; i <= j; i ++)
				{
					hv += h [i * m + j] * v [i * n + k];
				}
				v [(j + 1) * n + k] -= hv;
			}
			/* h_j+1,j */
			/* v_j+1 */
			hh = my_dnrm2 (n, v + (j + 1) * n, 1);
			my_dscal (n, 1.0 / hh, v + (j + 1) * n, 1);
			
			/* rotate */
			for (i = 0; i < j; i ++)
			{
				r1 = h [ i      * m + j];
				r2 = h [(i + 1) * m + j];
				h [ i      * m + j] = c [i] * r1 - s [i] * r2;
				h [(i + 1) * m + j] = s [i] * r1 + c [i] * r2;
			}
			rr = h [j * m + j];
			hv = sqrt (rr * rr + hh * hh); /* temporary variable */
			c [j] =  rr / hv;
			s [j] = -hh / hv;
			h [j * m + j] = hv; /* resultant (after rotated) element */
			
			g0 = g [j];
			g [j    ] = c [j] * g0;
			g [j + 1] = s [j] * g0;
		}
		/* 3. form the approximate solution */
		/* solve y_k */
		back_sub (j, m, h, g, c); /* use c[m] for y_k */
		// yv[n] = y_k v_k[]
		for (i = 0; i < n; i ++)
		{
			yv[i] = 0.0;
		}
		for (k = 0; k < j; k ++)
		{
			my_daxpy (n, c[k], v + k*n, 1, yv, 1);
		}
		// Kv[] = K^{-1}.(y_k v_k)
		inv (n, yv, Kv, inv_param);
		// x_m = x_0 + K^{-1}.(y_k v_k)
		my_daxpy (n, 1.0, Kv, 1, x, 1);
		
		/* 4. restart */
		res = fabs (g [j]); /* residual */
		/* if satisfied, */
		if (it->debug == 2)
		{
			fprintf (it->out, "libiter-gmres(%d) %d %d %e\n",
					 m, iter, j, res*res);
		}
		if (res <= eps)
		{
			ret = 0; // success
			break;
		}
		/* else */
		/* compute r_m */
		atimes (n, x, v + 0, atimes_param);
		/* r_m */
		/* compute v1 */
		
		my_daxpyz (n, -1.0, v + 0, 1, f, 1, v + 0, 1);
		
		g [0] = my_dnrm2 (n, v + 0, 1);
		my_dscal (n, 1.0 / g [0], v + 0, 1);
    }
	
# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H
	
	free (v);
	free (h);
	free (g);
	free (c);
	free (s);
	free (yv);
	free (Kv);
	
	/* adjust iter */
	iter *= m;
	
	if (it->debug == 1)
    {
		fprintf (it->out, "libiter-gmres(%d) %d %e\n",
				 m, iter, res*res / b2);
    }
	
	it->niter = iter;
	it->res2  = res * res / b2;
	return (ret);
}

/* plain GMRES (not restart version)
 * INPUT
 *   m : dimension of the problem
 *   b[m] : r-h-s vector
 *   atimes (int m, static double *x, double *b, void *param) :
 *        calc matrix-vector product A.x = b.
 *   atimes_param : parameters for atimes().
 *   it : struct iter. following entries are used
 *        it->max : max of iteration
 *        it->eps : criteria for |r^2|/|b^2|
 * OUTPUT
 *   returned value : 0 == success, otherwise (-1) == failed
 *   x[m] : solution
 *   it->niter : # of iteration
 *   it->res2  : |r^2| / |b^2|
 */
int
gmres (int n, const double *f, double *x,
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
	double scale;
	
# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H
	
	int ret = -1;
	int m = it->max;
	double eps = it->eps;
	
	double *v = (double *)malloc (sizeof (double) * (m + 1) * n);
	double *h = (double *)malloc (sizeof (double) * m * m);
	double *g = (double *)malloc (sizeof (double) * (m + 1));
	double *c = (double *)malloc (sizeof (double) * m);
	double *s = (double *)malloc (sizeof (double) * m);
	CHECK_MALLOC (v, "gmres");
	CHECK_MALLOC (h, "gmres");
	CHECK_MALLOC (g, "gmres");
	CHECK_MALLOC (c, "gmres");
	CHECK_MALLOC (s, "gmres");
	
	double res = 0.0;
	
	int i, j, k;
	double hv;
	double rr, hh;
	double r1, r2;
	double g0;
	
	
#ifdef HAVE_CBLAS_H
	/* use ATLAS' CBLAS routines */
	
	double b2 = cblas_ddot (n, f, 1, f, 1); // (f,f)
	eps *= sqrt (b2);
	
#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
	/* use Fortran BLAS routines */
	
	double b2 = ddot_ (&n, f, &i_1, f, &i_1); // (f,f)
	eps *= sqrt (b2);
	
# else // !HAVE_BLAS_H
	/* use local BLAS routines */
	
	double b2 = my_ddot (n, f, 1, f, 1); // (f,f)
	eps *= sqrt (b2);
	
# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H
	
	
	/* 1. start: */
	/* compute r0 */
	/* compute v1 */
	/* beta */
	atimes (n, x, v + 0, atimes_param); /* use v [0] temporaliry */
	
#ifdef HAVE_CBLAS_H
	/* use ATLAS' CBLAS routines */
	
	// v = f - v
	cblas_dscal (n, -1.0, v, 1); // v = - v
	cblas_daxpy (n, 1.0, f, 1, v, 1); // v = f - v
	
	g[0] = sqrt (cblas_ddot (n, v, 1, v, 1));
	cblas_dscal (n, 1.0 / g[0], v, 1);
	
#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
	/* use Fortran BLAS routines */
	
	// v = f - v
	dscal_ (&n, &d_m1, v, &i_1); // v = - v
	daxpy_ (&n, &d_1, f, &i_1, v, &i_1); // v = f - v
	
	g[0] = sqrt (ddot_ (&n, v, &i_1, v, &i_1));
	scale = 1.0 / g[0];
	dscal_ (&n, &scale, v, &i_1);
	
# else // !HAVE_BLAS_H
	/* use local BLAS routines */
	
	my_daxpyz (n, -1.0, v + 0, 1, f, 1, v + 0, 1);
	
	g [0] = my_dnrm2 (n, v + 0, 1);
	my_dscal (n, 1.0 / g [0], v + 0, 1);
	
# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H
	
	/* main loop */
	/* 2. iterate: */
	for (j = 0; j < m; j ++)
    {
		/* tmp = A.vj : use v [(j + 1) * n] directly */
		atimes (n, v + j * n, v + (j + 1) * n, atimes_param);
		/* h_i,j (i=1,...,j) */
		for (i = 0; i <= j; i ++)
		{
#ifdef HAVE_CBLAS_H
			/* use ATLAS' CBLAS routines */
			
			h [i * m + j] =
			cblas_ddot (n, v + (j + 1) * n, 1,
						v + i * n, 1);
			
#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
			/* use Fortran BLAS routines */
			
			h [i * m + j] =
			ddot_ (&n, v + (j + 1) * n, &i_1,
				   v + i * n, &i_1);
			
# else // !HAVE_BLAS_H
			/* use local BLAS routines */
			
			h [i * m + j] =
			my_ddot (n, v + (j + 1) * n, 1,
					 v + i * n, 1);
			
# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H
		}
		/* vv_j+1 */
		for (k = 0; k < n; k ++)
		{
			hv = 0.0;
			for (i = 0; i <= j; i ++)
			{
				hv += h [i * m + j] * v [i * n + k];
			}
			v [(j + 1) * n + k] -= hv;
		}
		/* h_j+1,j */
		/* v_j+1 */
#ifdef HAVE_CBLAS_H
		/* use ATLAS' CBLAS routines */
		
		//hh = cblas_dnrm2 (n, v + (j + 1) * n, 1);
		hh = sqrt (cblas_ddot (n, v + (j + 1) * n, 1, v + (j + 1) * n, 1));
		cblas_dscal (n, 1.0 / hh, v + (j + 1) * n, 1);
		
#else // !HAVE_CBLAS_H
# ifdef HAVE_BLAS_H
		/* use Fortran BLAS routines */
		
		//hh = dnrm2_ (&n, v + (j + 1) * n, &i_1);
		hh = sqrt (ddot_ (&n, v + (j + 1) * n, &i_1, v + (j + 1) * n, &i_1));
		scale = 1.0 / hh;
		dscal_ (&n, &scale, v + (j + 1) * n, &i_1);
		
# else // !HAVE_BLAS_H
		/* use local BLAS routines */
		
		hh = my_dnrm2 (n, v + (j + 1) * n, 1);
		my_dscal (n, 1.0 / hh, v + (j + 1) * n, 1);
		
# endif // !HAVE_BLAS_H
#endif // !HAVE_CBLAS_H
		
		/* rotate */
		for (i = 0; i < j; i ++)
		{
			r1 = h [ i      * m + j];
			r2 = h [(i + 1) * m + j];
			h [ i      * m + j] = c [i] * r1 - s [i] * r2;
			h [(i + 1) * m + j] = s [i] * r1 + c [i] * r2;
		}
		rr = h [j * m + j];
		hv = sqrt (rr * rr + hh * hh); /* temporary variable */
		c [j] =  rr / hv;
		s [j] = -hh / hv;
		h [j * m + j] = hv; /* resultant (after rotated) element */
		
		g0 = g [j];
		g [j    ] = c [j] * g0;
		g [j + 1] = s [j] * g0;
		
		res = fabs (g [j + 1]); /* residual */
		if (it->debug == 2)
		{
			fprintf (it->out, "libiter-gmres %d %e\n",
					 j, res*res);
		}
		/* if satisfied, */
		if (res <= eps)
		{
			j ++; /* this is because ++(*iter) in gmres(m) */
			ret = 0; // success
			break;
		}
    }
	
	/* 3. form the approximate solution */
	/* solve y_k */
	back_sub (j, m, h, g, c); /* use c [] as y_k */
	/* x_m */
	for (i = 0; i < n; i ++)
    {
		for (k = 0; k < j; k ++)
		{
			x [i] += v [k * n + i] * c [k];
		}
    }
	
	free (v);
	free (h);
	free (g);
	free (c);
	free (s);
	
	if (it->debug == 1)
    {
		fprintf (it->out, "libiter-gmres(%d) %d %e\n",
				 m, j, res*res / b2);
    }
	
	it->niter = j;
	it->res2  = res * res / b2;
	return (ret);
}

