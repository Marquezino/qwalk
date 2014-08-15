/* QWalk (qwstatistics_io.h) 
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

#ifndef _QWSTATISTICS_IO
#define _QWSTATISTICS_IO

#include "qwstatistics.h"


/* This function receives as input a string containing the name of the 
 * statistics file that will be written. It also receives a strucuture
 * containing the statistics concerning a certain step of the simulation.
 * In the first call of this function the file is created and the 
 * statistics are written. In subsequent calls of the function, the 
 * statistics are appended in the file. The filename must be the same in 
 * all the calls of the function or inexpected behaviour may occur. If 
 * the user wants to start writting statistics in a different file then 
 * this function must be called once with a blank string in the filename 
 * ("", open and close quotes without space).
 *
 * Error numbers:
 *   0: success (no error)
 *   1: could not open file
 */
int writeStatistics(const char *filename, statistics_t stat);


#endif

