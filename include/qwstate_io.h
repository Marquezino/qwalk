/* QWalk (qwstate_io.h) 
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

#ifndef _QWSTATE_IO
#define _QWSTATE_IO

#include "qwoptions_io.h"

/* This function receives as input the name of the file that contains
 * the definition of the state. It also receives a positive integer
 * describing the size of the lattice and an unsigned int
 * describing the type of lattice. We consider that the LINE lattice
 * coordinates ranges from -max to max and the CYCLE and SEGMENT
 * lattices range from 0 to max-1. If the input file is correct,
 * this function returns a complex 2D-matrix corresponding to the
 * state. Otherwise, the function returns NULL.
 */
double complex **readStateFile1D(const char *filename, int max,
				 unsigned int lattType);


/* This function receives as input the name of the file that contains
 * the definition of the state. It also receives a positive integer
 * describing the size of the lattice. We consider that the lattice
 * coordinates range from -max to max. If the input file is correct,
 * this function returns a complex 4D-matrix corresponding to the
 * state. Otherwise, the function returns NULL.
 */
double complex ****readStateFile2D(const char *filename, int max,
				   unsigned int lattType);



/* This function receives as input a string containing the name of the file
 * that will be written. It also receives an array with the amplitudes and a
 * structure containing the simulation options. The function writes the
 * amplitudes in the file, making the appropriate conversions of 
 * coordinates whenever necessary.
 * The function also writes a header in the file, describing the options
 * used in the simulation.
 *
 * Error numbers:
 *   0: success
 *   1: could not open file
 */
int writeState1D(const char *filename, double complex **wave, options1D_t options);



/* This function receives as input a string containing the name of the file
 * that will be written. It also receives an array with the amplitudes and a
 * structure containing the simulation options. The function writes the
 * amplitudes in the file, making the appropriate conversions of 
 * coordinates whenever necessary.
 * The function also writes a header in the file, describing the options
 * used in the simulation.
 *
 * Error numbers:
 *   0: success
 *   1: could not open file
 */
int writeState2D(const char *filename, double complex ****wave, options2D_t options);




#endif

