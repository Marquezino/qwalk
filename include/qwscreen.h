/* QWalk (qwscreen.h) 
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

#ifndef _QWSCREEN
#define _QWSCREEN

typedef struct{
  int xa;
  int xb;
  int ya;
  int yb;
  int xvar;
  int yvar;
  int numpts;
  double *values;
}screen_t;


/* This function receives as input the address of a screen_t structure and 
 * four integers. The four integers describe the limits of the screen
 * detector, meaning that the screen goes from point (xlower,ylower) to
 * point (xhigher,yhigher). The detector must be placed in horizontal,
 * vertical or in 45 degrees. This function initializes the structure
 * screen_t according to the points passed.
 *
 * Error numbers
 *   0: operation successful (no error)
 *   1: invalid screen address
 *   2: invalid slope of screen
 *   3: not enough memory
 */
int initScreen(screen_t *screen, int xa, int ya, int xb, int yb);


/* This function receives as input the address of a structure screen_t,
 * a 4D complex matrix and a positive integer. The screen_t structure 
 * must be initialized by initScreen function. The complex matrix 
 * represents the state of a walker in a 2D lattice. The integer max 
 * describes the size of the lattice, considering that the lattice ranges
 * from -max to max.
 *
 * Error numbers
 *   0: successful
 *   1: invalid screen address
 *   2: invalid state matrix
 *   3: invalid lattice size
 */
int updateScreen(screen_t *screen, double complex ****state, int max);


#endif

