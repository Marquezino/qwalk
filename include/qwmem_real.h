/* QWalk (qwmem_real.h) 
 * Copyright (C) 2008  Franklin Marquezino
 *                                                                                                                           
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *                                                                                                                           
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
 * 02110-1301, USA
 */

#ifndef _QWMEM_REAL
#define _QWMEM_REAL


/* This function receives two integers as arguments, and returns a pointer
 * to a 2D real (double precision) matrix, with dimensions according to 
 * those integers passed. If the dimensions passed are not valid, or if 
 * for any other reason the matrix cannot be allocated, then the function 
 * returns NULL.
 */
double **allocReal2D(int dim1, int dim2);


/* This function receives one integer as argument, and returns a pointer
 * to a 1D real (double precision) matrix, with dimensions according to 
 * this integer passed. If the dimension passed is not valid, or if for 
 * any other reason the matrix cannot be allocated, then the function 
 * returns NULL.
 */
double *allocReal1D(int dim);


/* This function receives as arguments a pointer to a 2D real (double 
 * precision) matrix and an integer describing the first dimension of 
 * the matrix. The function frees the memory occupied by the matrix.
 *
 * Error numbers
 *   0: success
 *   1: invalid dimensions
 *   2: invalid matrix
 */
int freeReal2D(double **matrix, int dim1);


/* This function receives as arguments a pointer to a 2D real (double 
 * precision) matrix and two integers describing the dimensions of the
 * matrix. It assigns zero to each entry.
 *
 * Error numbers
 *   0: successful
 *   1: invalid dimensions
 *   2: invalid matrix
 */
int cleanReal2D(double **matrix, int dim1, int dim2);


/* This function receives as arguments a pointer to a 1D real (double 
 * precision) vector and one integer describing the dimension of the
 * vector. It assigns zero to each entry of the vector.
 *
 * Error numbers
 *   0: success
 *   1: invalid dimensions
 *   2: invalid matrix
 */
int cleanReal1D(double *vector, int dim);


/* This function receives as arguments two pointers to 2D real (double 
 * precision) matrices and two integers. It copies the second matrix 
 * into the first one. Both matrices must have the same dimensions, which 
 * is indicated by the two integers passed.
 *
 * Error numbers
 *   0: success
 *   1: invalid dimensions
 *   2: invalid matrix
 */
int copyReal2D(double **dest, double **src, int dim1, int dim2);


#endif

