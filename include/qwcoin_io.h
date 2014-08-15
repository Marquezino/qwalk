/* QWalk (qwcoin_io.h) 
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

#ifndef _QWCOIN_IO
#define _QWCOIN_IO


/* This function receives as input the name of the file that
 * contains the definition of the coin. If the input file
 * is correct, the function returns a 2D complex matrix
 * corresponding to the coin. Otherwise, the function returns
 * NULL. In this version the function doesn't check if the
 * matrix entered by the user is unitary.
 */
double complex **readCoinFile1D(const char *filename);


/* This function receives as input the name of the file that
 * contains the definition of the coin. If the input file
 * is correct, the function returns a 4D complex matrix
 * corresponding to the coin. Otherwise, the function returns
 * NULL. In this version the function doesn't check if the
 * matrix entered by the user is unitary.
 */
double complex ****readCoinFile2D(const char *filename);

#endif

