/* QWalk (qw1d_sub.h) 
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

#ifndef _QW1D_SUB
#define _QW1D_SUB

#include<complex.h>
#include "qwoptions_io.h"
#include "qwextra_io.h"
#include "qwconsts.h"

/* This subroutine sets the coin for a 1D simulation. It receives the 
 * address of the matrix that will be used to store the coin, an 
 * integer describing the type of the coin and the name of the input 
 * file. It sets the matrix according to the type of the coin (if
 * it is CUSTOM then the input file is be read to get the complete
 * description of the matrix).
 */
void setCoin1D(double complex ***C, int coinType, const char *filename);



/* This subroutine sets the initial state for a 1D simulation. It receives
 * the address of the array that will be used to store the state, a 
 * structure options1D_t with the simulation options and 
 * the name of the input file. It sets the array according to the size of
 * lattice and the type of the state (if it is CUSTOM then the input file 
 * is read to get the complete description of the state).
 */
void setState1D(double complex ***A, options1D_t options, 
		const char *filename);



/* This subroutine breaks random links of a lattice previously initialized.
 * The probability of breaking each link is given by blProb. It receives
 * the address of the matrix used to define the broken link topology,
 * a structure options1D_t with the simulation options. Broken links are 
 * defined according to Physical Review A, 74, 012312 (2006).
 */
void randomBrokenLink1D(int ***B, options1D_t options);


/* This subroutine performs one iteration in the quantum walk. It receives
 * the address of the matrix used to store the state, the address of a
 * temporary matrix used in the iteration, the coin matrix, the matrix
 * of broken links, a structure options1D_t with the simulation options
 *  and an integer describing the number of the iteration.
 */
void iterate1D(double complex ***A, double complex ***Atemp, double complex **C, 
	       int **BrokenLinks, options1D_t options, int iteration);



/* This subroutine obtains the statistics for an interation of the walk
 * and saves it into the appropriate file. When many experiments are
 * being carried out the subroutine takes the average of the results.
 * It receives as arguments the array used to store the state, a 
 * structure options1D_t used to store the simulation options, a
 * structure filenames_t used to store the names of files generated by
 * the simulation, an integer describing the number of the iteration,
 * an integer describing the number of experiments carried out.
 */
void doStatistics1D(double complex **A, double *StatProb, options1D_t options, 
		    filenames_t fnames, int iteration, int experiment);



/* This subroutine performs all checks requested by input file. It
 * receives the matrix used to store the state, a structure options1D_t
 * used to store the simulation options and an integer describing the
 * number of the iteration.
 */
void check1D(double complex **A, options1D_t options, int iteration);



/* This function calculates the approximate stationary distribution with
 * a certain number of steps. It should not be used with unitary decoherence
 * generated by random broken links. The function receives a matrix used
 * to store the state, a matrix used to store the coin, the matrix of broken
 * links (only because it will be used internally, to call the function to
 * perform the iterations), and a structure options1D_t used to store the
 * simulation options.
 */
double *getStationary1D(double complex **A, double complex **C, 
			int **BLinks, options1D_t options);

#endif
