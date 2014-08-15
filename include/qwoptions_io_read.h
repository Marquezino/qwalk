/* QWalk (qwoptions_io_read.h) 
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

#ifndef _QWOPTIONS_IO_READ
#define _QWOPTIONS_IO_READ

#include "qwoptions_io.h"

int readOptions_coin2D(FILE *in, options2D_t *options);
int readOptions_state2D(FILE *in, options2D_t *options);
int readOptions_steps2D(FILE *in, options2D_t *options);
int readOptions_afterm2D(FILE *in, options2D_t *options);
int readOptions_cmix2D(FILE *in, options2D_t *options);
int readOptions_check2D(FILE *in, options2D_t *options);
int readOptions_blprob2D(FILE *in, options2D_t *options);
int readOptions_dtprob2D(FILE *in, options2D_t *options);
int readOptions_exp2D(FILE *in, options2D_t *options);
int readOptions_lsize2D(FILE *in, options2D_t *options);
int readOptions_lextra2D(FILE *in, options2D_t *options);
int readOptions_ltype2D(FILE *in, options2D_t *options);
int readOptions_detec2D(FILE *in, options2D_t *options);
int readOptions_seed2D(FILE *in, options2D_t *options);
int readOptions_screen2D(FILE *in, options2D_t *options);
int readOptions_blperm2D(FILE *in, options2D_t *options);

int readOptions_coin1D(FILE *in, options1D_t *options);
int readOptions_state1D(FILE *in, options1D_t *options);
int readOptions_steps1D(FILE *in, options1D_t *options);
int readOptions_check1D(FILE *in, options1D_t *options);
int readOptions_blprob1D(FILE *in, options1D_t *options);
int readOptions_dtprob1D(FILE *in, options1D_t *options);
int readOptions_exp1D(FILE *in, options1D_t *options);
int readOptions_seed1D(FILE *in, options1D_t *options);
int readOptions_lsize1D(FILE *in, options1D_t *options);
int readOptions_lextra1D(FILE *in, options1D_t *options);
int readOptions_ltype1D(FILE *in, options1D_t *options);
int readOptions_cmix1D(FILE *in, options1D_t *options);
int readOptions_detec1D(FILE *in, options1D_t *options);
int readOptions_afterm1D(FILE *in, options1D_t *options);

#endif
