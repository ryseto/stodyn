/* header file for my-cg.c --
 * CG method
 * Copyright (C) 2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: cg_.h,v 1.1 2007/11/25 18:46:41 kichiki Exp $
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
#ifndef	_CG__H_
#define	_CG__H_


/* conjugate gradient method, another implementation
 */
int
cg_ (int n, const double *b, double *x,
     void (*atimes) (int, const double *, double *, void *),
     void *atimes_param,
     struct iter *it);

/* conjugate gradient method with preconditioner
 */
int
cg_pc (int n, const double *b, double *x,
       void (*atimes) (int, const double *, double *, void *),
       void *atimes_param,
       void (*inv) (int, const double *, double *, void *),
       void *inv_param,
       struct iter *it);


#endif /* !_CG__H_ */
