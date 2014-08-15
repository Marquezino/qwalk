/* QWalk (qwlinks.h) 
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

#ifndef _QWLINKS
#define _QWLINKS


/* This function receives the addresses of one integer matrix, a 
 * positive integer and the type of lattice. The matrix indicates 
 * the broken links, according to Phys. Rev. A, 74, 012312 (2006). 
 * The positive integer describes the size of the lattice, 
 * which ranges from -max to max in the LINE lattice, and 
 * from 0 to max-1 in the CYCLE and SEGMENT lattices. 
 * The function initializes all the links closed.
 *
 * Error numbers
 *   0: success
 *   1: invalid lattice size
 *   2: invalid broken link matrix
 */
int initBrokenLink1D(int ***L1, int max, unsigned char type);


/* This function receives the addresses of two integer matrices, a 
 * positive integer and the type of lattice. The matrices indicate 
 * the broken links, according to Phys. Rev. A, 74, 012312 (2006). 
 * The positive integer describes the size of the lattice, 
 * which ranges from -max to max in both axes.
 * The function initializes all the links closed.
 *
 * Error numbers
 *   0: success
 *   1: invalid lattice size
 *   2: invalid broken link matrix
 */
int initBrokenLink2D(int *****L1, int *****L2, int max, unsigned char type);

#endif

