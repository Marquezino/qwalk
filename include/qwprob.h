/* QWalk (qwprob.h) 
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

#ifndef _QWPROB
#define _QWPROB


/* This function receives a complex matrix and an integer as arguments. 
 * The complex matrix must be a valid quantum state. It is a 4D-matrix 
 * because the first two indexes represent the coin and the two last 
 * indexes represent the position of the particle in a 2D lattice. The 
 * value in each entry of this matrix is a complex amplitude of the 
 * quantum state. The integer max, which is also passed as argument, 
 * describes the size of the lattice used in the simulation. We consider 
 * that the lattice ranges from -max to max, in the X axis, and from 
 * -max to max, in the Y axis. The function returns a matrix $P$ of real
 * numbers, where the value of $P_{mn}$ is the probability of finding 
 * the particle in the corresponding position of the lattice.
 * If the function cannot allocate enough memory for the matrix, or
 * if an invalid array is passed then it returns NULL.
 */
double **getProbArray2D(double complex ****matrix, options2D_t opts);


/* This function receives a complex matrix and a structure options1D_t. 
 * The complex matrix must be a valid quantum state. It is a 2D-matrix 
 * because the first index represent the coin and the last index 
 * represent the position of the particle in a 1D lattice. The 
 * value in each entry of this matrix is a complex amplitude of the 
 * quantum state. The options1D_t constains the simulation options. 
 * If the function cannot allocate enough memory for the vector, 
 * it returns NULL.
 */
double *getProbArray1D(double complex **matrix, options1D_t opts);


/* This function receives a real (double precision) matrix and an integer 
 * as arguments. The real matrix contains probabilities of finding the 
 * particle in each site, and the integer describes the size of the 
 * lattice. We consider that the lattice coordinates range from -max to
 * max. This function returns 1 if the sum of the probabilities of finding
 * the particle in each site is equal to 1, with a tolerance given by
 * WALK_TOL, and returns 0 otherwise.
 *
 * Returns
 *   0: not unitary
 *   1: unitary
 */
int checkProb2D(double **matrix, int max, unsigned int lattType);


/* This function receives a real (double precision) vector, an integer 
 * and an unsigned int as arguments. The real vector contains 
 * probabilities of finding the particle in each site, 
 * the integer describes the size of the lattice. We consider 
 * that the lattice coordinates range from -max to max in
 * the LINE lattice and from 0 to max-1 in the CYCLE and SEGMENT
 * lattices. The unsigned int contains the type of lattice.
 * This function returns 1 if the sum of the probabilities of 
 * finding the particle in each site is equal to 1, with a tolerance 
 * given by WALK_TOL, and returns 0 otherwise. The constant WALK_TOL 
 * can be modified in the header file walk.h.
 *
 * Returns
 *   0: not unitary
 *   1: unitary
 */
int checkProb1D(double *matrix, int max, unsigned int lattType);


/* This function receives the address of a real (double precision) vector, 
 * a complex matrix and a structure options1D_t. The real vector stores the
 * average probabilities of the experiments. This function must be called
 * always with the same address for the real vector. The complex matrix 
 * contains the quantum state. The options1D_t strucuture contains 
 * simulation options. In each call the probability of finding 
 * the particle in each site is calculated, divided by the number 
 * of experiments and the vector pointed by *aveMatrix is updated. In 
 * last call the vector pointed by *aveMatrix contains the average of 
 * the probabilities obtained in the experiments.
 *
 * Error numbers
 *   0: success
 *   1: invalid number of experiments; or function called more times than 
 *      the number of experiments expected.
 *   2: invalid address passed
 *   3: could not allocate memory for the vector of probabilities
 *   4: invalid matrix passed in subsequent calls
 */
int averageProbFromState1D(double **aveMatrix, double complex **state, 
			   options1D_t options);


/* This function receives the address of a real (double precision) matrix, 
 * a complex matrix and a structure options2D_t. The real matrix stores the
 * average probabilities of the experiments. This function must be called
 * always with the same address for the real matrix. The complex matrix 
 * contains the quantum state. The options2D_t strucuture contains 
 * simulation options. In each call the probability of 
 * finding the particle in each site is calculated, divided 
 * by the number of experiments and the matrix pointed by *aveMatrix 
 * is updated. In last call the matrix pointed by *aveMatrix 
 * contains the average of the probabilities obtained in the experiments.
 *
 * Error numbers
 *   0: success
 *   1: invalid number of experiments; or function called more times than 
 *      the number of experiments expected.
 *   2: invalid address passed
 *   3: could not allocate memory for the vector of probabilities
 *   4: invalid matrix passed in subsequent calls
 */
int averageProbFromState2D(double ***aveMatrix, double complex ****state, 
			   options2D_t options);

#endif

