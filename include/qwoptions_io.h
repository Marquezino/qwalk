/* QWalk (qwoptions_io.h) 
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

#ifndef _QWOPTIONS_IO
#define _QWOPTIONS_IO

typedef struct{
  unsigned char error;
  unsigned char coinType;
  unsigned char stateType;
  unsigned char lattType;
  float blProb;
  float dtProb;
  int steps;
  int max;
  int lattextra;
  int numOfExperiments;
  int stepsAfterMeasure;
  int stepsMix;
  unsigned char calcMix;
  unsigned char checkState;
  unsigned char checkSymmetry;
  int detectors;
  int *detector_pts;
  int seed;
}options1D_t;

typedef struct{
  unsigned char error;
  unsigned char coinType;
  unsigned char stateType;
  unsigned char blType;
  unsigned char lattType;
  float blProbA;
  float blProbB;
  float dtProb;
  int steps;
  int max;
  int lattextra;
  int numOfExperiments;
  int stepsAfterMeasure;
  int stepsMix;
  unsigned char calcMix;
  unsigned char checkState;
  unsigned char checkXSymmetry;
  unsigned char checkYSymmetry;
  unsigned char screen;
  int screen_pta[2];
  int screen_ptb[2];
  int detectors;
  int **detector_pts;
  int seed;
}options2D_t;


/* This function receives a string containing the name of the input file.
 * It reads the input, identifying keywords and storing the simulation
 * options in a structure options1D_t. The function returns a structure 
 * options1D_t containing all important options for simulation, such as 
 * coin type, inicial state, number of steps, size of lattice, etc.
 *
 * Error numbers:
 *   0: operation successful (no error)
 *   1: could not open input file
 *   2: invalid coin in input file.
 *   3: invalid state in input file.
 *   4: invalid number of steps in input file.
 *   5: could not understand CHECK option.
 *   6: invalid probability of broken links.
 *   7: invalid number of experiments
 *   8: error in BEGIN-END structure
 *   9: invalid lattice size
 *  10: invalid LATTEXTRA
 *  11: invalid lattice type
 *  12: invalid number of steps in mixing time calculation
 *  13: invalid probability (measuments or broken links)
 */
options1D_t readOptionsFile1D(const char *filename);


/* This function receives a string containing the name of the input file.
 * It reads the input, identifying keywords and storing the simulation
 * options in a structure options2D_t. The function returns the structure 
 * options2D_t containing all important options for simulation, such as 
 * coin type, inicial state, number of steps, lattice size, etc.
 *
 * Error numbers:
 *   0: success
 *   1: could not open input file
 *   2: invalid coin
 *   3: invalid state
 *   4: invalid number of steps
 *   5: invalid check option
 *   6: invalid number of experiments
 *   7: invalid lattice size
 *   8: invalid number of detectors
 *   9: error in BEGIN-END structure
 *  10: invalid LATTEXTRA option
 *  11: invalid probability (measuments or broken links)
 *  12: invalid number of steps in mixing time calculation
 *  13: invalid lattice type
 */
options2D_t readOptionsFile2D(const char *filename);

#endif

