/* QWalk (qwstatistics.h) 
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

#ifndef _QWSTATISTICS
#define _QWSTATISTICS

#include "qwoptions_io.h"

typedef struct{
  int iteration;
  float variance;
  float tvd;
  float tvdu;
  float meanX;
  float meanY;
}statistics_t;


/* This function receives a complex matrix representing a quantum state, a
 * real (double precision) matrix representing the approximate stationary
 * distribution, a structure options1D_t with simulation 
 * options, and an integer describing the number of the iteration. 
 * The function returns a structure statistics_t with 
 * the statistics for that iteration (variance, mean, etc). If something goes 
 * wrong, a negative error number is returned in the iteration field.
 *
 * Error numbers
 *   >0: success
 *   -1: invalid lattice size
 *   -2: invalid matrix
 *   -3: invalid iteration number
 */
statistics_t getStatisticsFromState1D(double complex **matrix, double *StatProb,
				      options1D_t opts, int iteration);


/* This function receives a complex matrix representing a quantum state, a
 * real (double precision) matrix representing the approximate stationary
 * distribution, a structure options2D_t representing the simulation
 * options, and an integer describing the number of the iteration. 
 * We consider that the lattice ranges from 
 * -max to max. The function returns a structure statistics_t with the
 * statistics for that iteration (variance, mean, etc). If something goes 
 * wrong, a negative error number is returned in the iteration field.
 *
 * Error numbers
 *   >0: success
 *   -1: invalid lattice size
 *   -2: invalid matrix
 *   -3: invalid iteration number
 */
statistics_t getStatisticsFromState2D(double complex ****matrix, double **StatProb,
				      options2D_t opts, int iteration);


/* This function receives a real (double precision) matrix containing the
 * probability of finding the particle in each site of the lattice, an
 * integer describing the size of the 1D-lattice, and an integer describing
 * the number of the iteration. We consider that the lattice ranges from 
 * -max to max. The function returns a structure statistics_t with the 
 * statistics for that iteration (variance, stardard deviation, etc). If
 * something wrong happens, a negative error number is returned in the
 * iteration field.
 *
 * Error numbers
 *   >0: success
 *   -1: invalid lattice size
 *   -2: invalid matrix
 *   -3: invalid iteration number
 */
statistics_t getStatisticsFromProb1D(double *matrix, int max, int iteration);


/* This function receives a real (double precision) matrix containing the
 * probability of finding the particle in each site of the lattice, an
 * integer describing the size of the 2D-lattice, and an integer describing
 * the number of the iteration. We consider that the lattice ranges from 
 * -max to max. The function returns a structure statistics_t with the 
 * statistics for that iteration (variance, stardard deviation, etc). If
 * something wrong happens, a negative error number is returned in the
 * iteration field.
 *
 * Error numbers
 *   >0: success
 *   -1: invalid lattice size
 *   -2: invalid matrix
 *   -3: invalid iteration number
 */
statistics_t getStatisticsFromProb2D(double **matrix, int max, int iteration);

#endif

