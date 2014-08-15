/* QWalk (qwstate.h) 
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

#ifndef _QWSTATE
#define _QWSTATE

#include "qwoptions_io.h"

/* This function receives an integer max as argument, where max describes the
 * size of the lattice used in the simulation. We consider that the lattice
 * ranges from -max to max, in the X axis, and from -max to max, in the Y axis.
 * The function returns a matrix with a quantum state that can be used as 
 * initial condition for a simulation, giving maximum dispersion for Grover 
 * coin. If max is lower or equal zero, of if for any reason the operating 
 * system cannot allocate enough memory, then we return NULL. Note that the 
 * resulting matrix is 4D. The first two indexes indicate the coin and the 
 * two last indexes indicate the position of the particle in a 2D lattice.
 * The value in each entry of this matrix is a complex amplitude of the 
 * quantum state.
 */
double complex ****createGroverState2D(int max, unsigned int lattType);


/* This function receives an integer max as argument, where max describes 
 * the size of the lattice used in the simulation. We consider that the 
 * lattice ranges from -max to max, in the X axis, and from -max to max, 
 * in the Y axis. The function returns a matrix with a quantum state that 
 * can be used as initial condition for a simulation, giving maximum 
 * dispersion for Hadamard coin. If max is lower or equal zero, of if for 
 * any reason the operating system cannot allocate enough memory, then we 
 * return NULL. Note that the resulting matrix is 4D. The first two indexes 
 * indicate the coin and the two last indexes indicate the position of 
 * the particle in a 2D lattice. The value in each entry of this matrix 
 * is a complex amplitude of the quantum state.
 */
double complex ****createHadamardState2D(int max, unsigned int lattType);


/* This function receives an integer max as argument, where max describes 
 * the size of the lattice used in the simulation. We consider that the 
 * lattice ranges from -max to max, in the X axis, and from -max to max, 
 * in the Y axis. The function returns a matrix with a quantum state that 
 * can be used as initial condition for a simulation, giving maximum 
 * dispersion for Fourier coin. If max is lower or equal zero, of if for 
 * any reason the operating system cannot allocate enough memory, then 
 * we return NULL. Note that the resulting matrix is 4D. The first two 
 * indexes indicate the coin and the two last indexes indicate the position 
 * of the particle in a 2D lattice. The value in each entry of this matrix 
 * is a complex amplitude of the quantum state.
 */
double complex ****createFourierState2D(int max, unsigned int lattType);


/* This function receives an integer max as argument, where max describes 
 * the size of the lattice used in the simulation, and an unsigned
 * int lattType describing the type of lattice. We consider that the 
 * lattice ranges from -max to max in the LINE lattice, and 
 * from 0 to max-1 in the CYCLE and SEGMENT lattices. 
 * The function returns a matrix with a quantum state that 
 * can be used as initial condition for a simulation, 
 * giving maximum dispersion for Hadamard coin. If max is lower or equal 
 * zero, of if for any reason the operating system cannot allocate enough 
 * memory, then we return NULL.
 */
double complex **createHadamardState1D(int max,  unsigned int lattType);


/* This function receives a complex matrix and an integer as arguments. The 
 * complex matrix must be a valid quantum state. It is a 4D-matrix because 
 * the first two indexes represent the coin and the two last indexes represent 
 * the position of the particle in a 2D lattice. The value in each entry of 
 * this matrix is a complex amplitude of the quantum state. The integer max, 
 * which is also passed as argument, describes the size of the lattice used 
 * in the simulation. We consider that the lattice ranges from -max to max, 
 * in the X axis, and from -max to max, in the Y axis. The function returns 
 * 1 if the state is valid, i.e., if it has norm 1 and returns 0 otherwise.
 *
 * Returns
 *   0: not unitary
 *   1: unitary
 */
int checkState2D(double complex ****matrix, options2D_t options);


/* This function receives a complex matrix and a structure options1D_t. The
 * complex matrix is a quantum state and the structure contains simulation
 * options.  This function returns 1 if the sum of the probabilities of 
 * finding the particle in each site is equal to 1, with a tolerance given 
 * by WALK_TOL, and returns 0 otherwise.
 *
 * Returns
 *   0: not unitary
 *   1: unitary
 */
int checkState1D(double complex **matrix, options1D_t options);


/* This function receives a complex array representing a quantum state
 * and a positive integer describing the size of the lattice. We consider
 * that the lattice ranges from -max to max. This function checks the 
 * symmetry of the probability array around x=0, i.e., it checks if the 
 * probabilities for the x>0 sites are equal to the probabilities for the 
 * x<0 sites. It returns 1 if the array is symmetrical, with a tolerance 
 * of WALK_TOL, and returns 0 if the matrix is not symmetrical.
 *
 * Returns
 *   0: not symmetrical
 *   1: symmetrical
 */
int checkSymmetry1D(double complex **matrix, int max);


/* This function receives a complex matrix representing a quantum state
 * and a positive integer describing the size of the lattice. We consider
 * that the lattice ranges from -max to max. This function checks the 
 * symmetry of the probability matrix around the y=0 axis, i.e., it 
 * checks if the probabilities for the x>0 sites are equal to the 
 * probabilities for the x<0 sites. It returns 1 if the matrix is 
 * symmetrical, with a tolerance of WALK_TOL, and returns 0 if the
 * matrix is not symmetrical.
 *
 * Returns
 *   0: not symmetrical
 *   1: symmetrical
 */
int checkXSymmetry2D(double complex ****matrix, int max);


/* This function receives a complex matrix representing a quantum state
 * and a positive integer describing the size of the lattice. We consider
 * that the lattice ranges from -max to max. This function checks the 
 * symmetry of the probability matrix around the x=0 axis, i.e., it 
 * checks if the probabilities for the y>0 sites are equal to the 
 * probabilities for the y<0 sites. It returns 1 if the matrix is 
 * symmetrical, with a tolerance of WALK_TOL, and returns 0 if the
 * matrix is not symmetrical.
 *
 * Returns
 *   0: not symmetrical
 *   1: symmetrical
 */
int checkYSymmetry2D(double complex ****matrix, int max);

#endif

