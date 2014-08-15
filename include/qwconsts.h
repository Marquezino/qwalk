/* QWalk (qwconsts.h) 
 * Copyright (C) 2007  Franklin Marquezino
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

#ifndef _QWCONSTS
#define _QWCONSTS

#define WALK_TOL 10e-5 /* it was 10e-8 until version 1.0 */
/* Tolerance to approximation errors */

#define MAXIMUM(A,B) ((A>B) ? (A):(B))
#define MINIMUM(A,B) ((A<B) ? (A):(B))

#define DELTA(A,B) ((A==B) ? 1:0) 
/* 1, if A=B; 0, otherwise*/

#define STREQ(A,B) !strcmp((A),(B)) 
/* 1, if strings are equal; 0, otherwise.*/

#define STRNEQ(A,B) strcmp((A),(B)) 
/* 1, if strings are not equal; 0, otherwise.*/

#define CUSTOM_COIN   10
#define FOURIER_COIN  11
#define GROVER_COIN   12
#define HADAMARD_COIN 13

#define CUSTOM_STATE   20
#define FOURIER_STATE  21
#define GROVER_STATE   22 
#define HADAMARD_STATE 23

#define NO_BROKENLINKS 30
#define PERMANENT_BROKENLINKS 31

#define NATURAL_LATT 40
#define DIAG_LATT 41
#define LINE_LATT 42
#define CYCLE_LATT 43
#define SEGMENT_LATT 44

#endif
