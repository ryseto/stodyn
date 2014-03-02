/* overall wrapper for iterative solver routines
 * Copyright (C) 2006-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: libiter.c,v 1.8 2007/11/25 19:07:25 kichiki Exp $
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
#include <cstdio> /* fprintf() */
#include <cstdlib> /* malloc(), free() */
#include <string.h> /* strcmp() */

#include "libiter.h"

#include "steepest.h"
#include "cg.h"
#include "cgs.h"
#include "bicgstab.h"

#include "cg_.h" // cg_pc()

#include "gmres.h"
#include "bi-cgstab.h"
#include "orthomin.h"

#include "memory-check.h" // CHECK_MALLOC


/* initialize parameters
 * INPUT
 *   solver : string indicating the solver
 *            sta, sta2, gpb, otmk, or gmres (default)
 *   max, restart, eps : iteration parameters
 *   n          : dimension of the problem
 *   guess[n]   : initial guess
 *                if NULL is given, set zero for the guess
 *   flag_guess = 0 : don't use the guess[] in solve_iter()
 *              = 1 : use guess[] for the initial guess in solve_iter()
 *   debug = 0 : no debug info
 *         = 1 : iteration numbs and residue
 *         = 2 : residue for each iteration step
 *   out   : FILE * to output debug info.
 * OUTPUT
 *   (struct iter *) returned value :
 */
struct iter *
iter_init (const char *solver,
		   int max, int restart, double eps,
		   int n, const double *guess, int flag_guess,
		   int debug, FILE *out)
{
	struct iter *param = (struct iter *)malloc (sizeof (struct iter));
	CHECK_MALLOC (param, "iter_init");
	
	param->solver = (char *)malloc (sizeof (char) * (strlen (solver) + 1));
	CHECK_MALLOC (param->solver, "iter_init");
	strcpy (param->solver, solver);
	
	param->max = max;
	param->restart = restart;
	param->eps = eps;
	
	/* to keep good initial guess for the next step */
	param->n = n;
	param->guess = (double *)malloc (sizeof (double) * n);
	CHECK_MALLOC (param->solver, "iter_init");
	int i;
	for (i = 0; i < n; i ++)
    {
		if (guess == NULL)
		{
			param->guess[i] = 0.0;
		}
		else
		{
			param->guess[i] = guess[i];
		}
    }
	param->flag_guess = flag_guess;
	
	/* results of the iteration */
	param->niter = 0;
	param->res2  = 0.0;
	
	/* to report the debug informations */
	param->out   = out;
	param->debug = debug;
	
	return (param);
}

void
iter_free (struct iter *param)
{
	if (param == NULL) return;
	
	if (param->solver != NULL) free (param->solver);
	if (param->guess  != NULL) free (param->guess);
	free (param);
}

/* wrapper routine for iterative solvers
 * INPUT
 *   n : size of vectors v[] and f[] -- expected to be np * nelm for red-sym
 *   b [n] : given vector
 *   atimes (int n, static double *x, double *b, void *param) :
 *        calc matrix-vector product A.x = b.
 *   atimes_param : parameters for atimes().
 *   it : parameters for iterative solvers
 *        solver : string indicating the solver
 *          "steepest" : steepest descent method
 *          "cg"       : conjugate gradient
 *          "cgs"      : conjugate gradient squared
 *          "bicgstab" : bi-conjugate gradient stabilized
 *          "sta", "sta2", "gpb", "otmk" :
 *          "gmres"    : generalized minimum residual method  (default)
 *        max, restart, eps
 *        n, guess[n] : the result at the last process
 *        flag_guess : 0 == don't keep the results,
 *                           1 == keep the results for the next.
 * OUTPUT
 *   returned value : 0 == success, otherwise (-1) == failed
 *   x [n] : solution
 */
int
solve_iter (int n, const double *b,
			double *x,
			void (*atimes) (int, const double *, double *, void *),
			void *atimes_param,
			struct iter *it)
{
	int i;
	if (it->n != n)
    {
		fprintf (stderr, "libiter solve_iter : n is different %d != %d\n",
				 it->n, n);
		exit (1);
    }
	
	// set the initial guess
	if (it->flag_guess == 0)
    {
		for (i = 0; i < n; ++i)
		{
			x[i] = 0.0;
		}
    }
	else
    {
		for (i = 0; i < n; ++i)
		{
			x[i] = it->guess[i];
		}
    }
	
	int ret = -1;
	if (strcmp (it->solver, "gmres") == 0)
    {
		ret = gmres_m (n, b, x,
					   atimes, atimes_param,
					   it);
    } else if (strcmp (it->solver, "steepest") == 0)
    {
		ret = steepest (n, b, x,
						atimes, atimes_param,
						it);
    }
	else if (strcmp (it->solver, "cg") == 0)
    {
		ret = cg (n, b, x,
				  atimes, atimes_param,
				  it);
    }
	else if (strcmp (it->solver, "cgs") == 0)
    {
		ret = cgs (n, b, x,
				   atimes, atimes_param,
				   it);
    }
	else if (strcmp (it->solver, "bicgstab") == 0)
    {
		ret = bicgstab (n, b, x,
						atimes, atimes_param,
						it);
    }
	else if (strcmp (it->solver, "sta") == 0)
    {
		ret = sta (n, b, x,
				   atimes, atimes_param,
				   it);
    }
	else if (strcmp (it->solver, "sta2") == 0)
    {
		ret = sta2 (n, b, x,
					atimes, atimes_param,
					it);
    }
	else if (strcmp (it->solver, "gpb") == 0)
    {
		ret = gpb (n, b, x,
				   atimes, atimes_param,
				   it);
    }
	else if (strcmp (it->solver, "otmk") == 0)
    {
		ret = otmk (n, b, x,
					atimes, atimes_param,
					it);
    }
	else
    {
		fprintf (stderr, "libiter-solve_iter : "
				 "invalid solver %s.\n", it->solver);
    }
	
	/* keep the results for the next first guess */
	if (it->flag_guess != 0)
    {
		for (i = 0; i < n; ++i)
		{
			it->guess[i] = x[i];
		}
    }
	return (ret);
}


/* wrapper routine for iterative solvers with preconditioner
 * INPUT
 *   n : size of vectors v[] and f[] -- expected to be np * nelm for red-sym
 *   b [n] : given vector
 *   atimes (int n, static double *x, double *b, void *param) :
 *        calc matrix-vector product A.x = b.
 *   atimes_param : parameters for atimes().
 *   inv (int n, static double *b, double *x, void *param) :
 *        approx of A^{-1}.b = x for preconditioning.
 *   inv_param : parameters for the preconditioner inv().
 *   it : parameters for iterative solvers
 *        solver : string indicating the solver
 *          "steepest" : steepest descent method
 *          "cg"       : conjugate gradient
 *          "cgs"      : conjugate gradient squared
 *          "bicgstab" : bi-conjugate gradient stabilized
 *          "sta", "sta2", "gpb", "otmk" :
 *          "gmres"    : generalized minimum residual method  (default)
 *        max, restart, eps
 *        n, guess[n] : the result at the last process
 *        flag_guess : 0 == don't keep the results,
 *                           1 == keep the results for the next.
 * OUTPUT
 *   returned value : 0 == success, otherwise (-1) == failed
 *   x [n] : solution
 */
int
solve_iter_pc (int n, const double *b,
			   double *x,
			   void (*atimes) (int, const double *, double *, void *),
			   void *atimes_param,
			   void (*inv) (int, const double *, double *, void *),
			   void *inv_param,
			   struct iter *it)
{
	int i;
	if (it->n != n)
    {
		fprintf (stderr, "libiter solve_iter_pc : n is different %d != %d\n",
				 it->n, n);
		exit (1);
    }
	if (inv == NULL)
    {
		// no preconditioning
		return (solve_iter (n, b, x, atimes, atimes_param, it));
    }
	
	// set the initial guess
	if (it->flag_guess == 0)
    {
		for (i = 0; i < n; ++i)
		{
			x[i] = 0.0;
		}
    }
	else
    {
		for (i = 0; i < n; ++i)
		{
			x[i] = it->guess[i];
		}
    }
	
	int ret = -1;
	if (strcmp (it->solver, "steepest") == 0)
    {
		fprintf (stderr, "libiter-solve_iter_pc : "
				 "preconditioned steepest is not implemented.\n");
    }
	else if (strcmp (it->solver, "cg") == 0)
    {
		ret = cg_pc (n, b, x,
					 atimes, atimes_param,
					 inv, inv_param,
					 it);
    }
	else if (strcmp (it->solver, "cgs") == 0)
    {
		fprintf (stderr, "libiter-solve_iter_pc : "
				 "preconditioned cgs is not implemented.\n");
    }
	else if (strcmp (it->solver, "bicgstab") == 0)
    {
		fprintf (stderr, "libiter-solve_iter_pc : "
				 "preconditioned bicgstab is not implemented.\n");
    }
	else if (strcmp (it->solver, "sta") == 0)
    {
		ret = sta_pc (n, b, x,
					  atimes, atimes_param,
					  inv, inv_param,
					  it);
    }
	else if (strcmp (it->solver, "sta2") == 0)
    {
		ret = sta2_pc (n, b, x,
					   atimes, atimes_param,
					   inv, inv_param,
					   it);
    }
	else if (strcmp (it->solver, "gpb") == 0)
    {
		ret = gpb_pc (n, b, x,
					  atimes, atimes_param,
					  inv, inv_param,
					  it);
    }
	else if (strcmp (it->solver, "otmk") == 0)
    {
		ret = otmk_pc (n, b, x,
					   atimes, atimes_param,
					   inv, inv_param,
					   it);
    }
	else if (strcmp (it->solver, "gmres") == 0)
    {
		ret = gmres_m (n, b, x,
					   atimes, atimes_param,
					   it);
    }
	else
    {
		fprintf (stderr, "libiter-solve_iter_pc : "
				 "invalid solver %s.\n", it->solver);
    }
	
	/* keep the results for the next first guess */
	if (it->flag_guess != 0)
    {
		for (i = 0; i < n; ++i)
		{
			it->guess[i] = x[i];
		}
    }
	return (ret);
}
