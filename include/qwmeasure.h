/* QWalk (qwmeasure.h) 
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

#ifndef _QWMEASURE
#define _QWMEASURE

#include "qwoptions_io.h"

/* This function receives the address of a complex matrix with the quantum
 * state, the address of a temporary complex matrix and a structure with
 * the simulation options. The first dimension of the array
 * gives the number of the detector, while the second dimension gives the
 * x and y coordinates. For example, detectors_pts[3][0] is the x coordinate
 * of the third detector, while detectors_pts[3][1] is the y coordinate.
 * If the operation is successful it returns the result of the measurement 
 * (i.e., the number of the detector where the particle was found). 
 * Otherwise, the function returns a negative error number.
 *
 * Error numbers
 *   >=0: success
 *   -1: invalid matrix passed
 *   -2: invalid lattice size
 *   -3: invalid number of detectors
 *   -4: invalid array of detector coordinates
 *   -5: not enough memory
 *   -6: could not clean array
 */
int measureState2D(double complex *****A, double complex *****Atemp, 
		   options2D_t opts);


/* Performs random measurements over a 2D lattice.
 *
 */
int randMeasure2D(double complex *****A, options2D_t opts);



/* Performs random measurements over a 1D lattice.
 *
 */
int randMeasure1D(double complex ***A, options1D_t opts);


/* This function receives the address of a complex matrix with the quantum
 * state, the address of a temporary complex matrix and a structure with
 * the simulation options. 
 * If the operation is successful it returns the result of the measurement 
 * (i.e., the number of the detector where the particle was found). 
 * Otherwise, the function returns a negative error number.
 *
 * Error numbers
 *   >=0: success
 *   -1: invalid matrix passed
 *   -2: invalid lattice size
 *   -3: invalid number of detectors
 *   -4: invalid array of detector coordinates
 *   -5: not enough memory
 *   -6: could not clean array
 */
int measureState1D(double complex ***A, double complex ***Atemp, 
		   options1D_t opts);


#endif

