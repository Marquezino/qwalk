/* QWalk (qwmem_int.h) 
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

#ifndef _QWMEM_INT
#define _QWMEM_INT


/* This function receives four integers as arguments, and returns a pointer
 * to a 4D integer matrix, with dimensions according to those integers 
 * passed. If the dimensions passed are not valid, or if for any other
 * reason the matrix cannot be allocated, then the function returns NULL.
 */ 
int ****allocInt4D(int dim1, int dim2, int dim3, int dim4);


/* This function receives two integers as arguments, and returns a pointer
 * to a 2D integer matrix, with dimensions according to those integers passed. 
 * If the dimensions passed are not valid, or if for any other reason the 
 * matrix cannot be allocated, then the function returns NULL.
 */
int **allocInt2D(int dim1, int dim2);


/* This function receives as arguments a pointer to a 4D integer 
 * matrix and three integers describing the first three dimensions 
 * of the matrix. The function frees the memory occupied by the 
 * matrix.
 *
 * Error numbers
 *   0: success
 *   1: invalid dimensions
 *   2: invalid matrix
 */
int freeInt4D(int ****array, int dim1, int dim2, int dim3);


/* This function receives as arguments a pointer to a 2D integer 
 * matrix and an integer describing the first dimension of 
 * the matrix. The function frees the memory occupied by the matrix.
 *
 * Error numbers
 *   0: success
 *   1: invalid dimensions
 *   2: invalid matrix
 */
int freeInt2D(int **matrix, int dim1);


/* This function receives as arguments a pointer to a 4D integer matrix and
 * four integers describing the dimensions of the matrix. It assigns zero
 * to each entry of the matrix.
 *
 * Error numbers
 *   0: success
 *   1: invalid dimensions
 *   2: invalid matrix
 */
int cleanInt4D(int ****matrix, int dim1, int dim2, int dim3, int dim4);


/* This function receives as arguments a pointer to a 2D integer 
 * matrix and two integers describing the dimensions of the
 * matrix. It assigns zero to each entry of the matrix.
 *
 * Error numbers
 *   0: successful
 *   1: invalid dimensions
 *   2: invalid matrix
 */
int cleanInt2D(int **vector, int dim1, int dim2);

#endif

