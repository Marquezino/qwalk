/* QWalk (qwcoin.h) 
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

#ifndef _QWCOIN
#define _QWCOIN

/* This function receives no argument and returns a pointer to a matrix which 
 * represents the Hadamard coin for a two-dimensional lattice. If there is not
 * enough memory or if, for any other reason, the operating system cannot 
 * allocate enough memory for the matrix, then this function returns NULL.
 * Note that this matrix is stored as a 4-dimensional one. However, if we read 
 * the first two indexes as a binary number we may understand them as meaning 
 * the "line" of the matrix, in the usual mathematical notation, and similarly 
 * the two last indexes meaning the "column".
 */
double complex ****createHadamardCoin2D();


/* This function receives no argument and returns a pointer to a matrix which 
 * represents the Fourier coin for a two-dimensional lattice. If there is not
 * enough memory or if, for any other reason, the operating system cannot 
 * allocate enough memory for the matrix, then this function returns NULL.
 * Note that this matrix is stored as a 4-dimensional one. However, if we read 
 * the first two indexes as a binary number we may understand them as meaning 
 * the "line" of the matrix, in the usual mathematical notation, and similarly 
 * the two last indexes meaning the "column".
 */
double complex ****createFourierCoin2D();


/* This function receives no argument and returns a pointer to a matrix which 
 * represents the Grover coin for a two-dimensional lattice. If there is not
 * enough memory or if, for any other reason, the operating system cannot 
 * allocate enough memory for the matrix, then this function returns NULL.
 * Note that this matrix is stored as a 4-dimensional one. However, if we read 
 * the first two indexes as a binary number we may understand them as meaning 
 * the "line" of the matrix, in the usual mathematical notation, and similarly 
 * the two last indexes meaning the "column".
 */
double complex ****createGroverCoin2D();


/* This function receives no argument and returns a pointer to a matrix which 
 * represents the Hadamard coin for a one-dimensional lattice. If there is not
 * enough memory or if, for any other reason, the operating system cannot 
 * allocate enough memory for the matrix, then this function returns NULL.
 */
double complex **createHadamardCoin1D();

#endif

