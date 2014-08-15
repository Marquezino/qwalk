/* QWalk (qwprob_io.h) 
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

#ifndef _QWPROB_IO
#define _QWPROB_IO

#include "qwoptions_io.h"


/* This function receives as input a string containing the name of the file
 * that will be written. It also receives a vector of probabilities and a
 * structure containing the simulation options. The function writes the
 * probabilities in the file, making the appropriate conversions of 
 * coordinates whenever necessary.
 * The function also writes a header in the file, describing the options
 * used in the simulation.
 *
 * Error numbers:
 *   0: success
 *   1: could not open file
 *   2: invalid probability vector
 */
int writeData1D(const char *filename, double *probs, options1D_t options);


/* This function receives as input a string containing the name of the data
 * file that will be written. It also receives a matrix of probabilities and
 * a structure containing the simulation options. The function writes the
 * probabilities in the file, making the appropriate conversions of 
 * coordinates whenever necessary.
 * The function also writes a header in the file, describing the options
 * used in the simulation.
 *
 * Error numbers:
 *   0: success
 *   1: could not open file
 *   2: invalid probability matrix
 */
int writeData2D(const char *filename, double **probs, options2D_t options);

#endif
