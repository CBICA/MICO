/**
 * @file  mico.h
 * @brief Multiplicative Intrinsic Component Optimization (MICO).
 *
 * Copyright (c) 2011, 2012 University of Pennsylvania. All rights reserved.<br />
 * See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 *
 * @ingroup PublicAPI
 */

#pragma once
#ifndef _MICO_H_
#define _MICO_H_


#include <mico/io.h>              // I/O functions, e.g., to parse list file
#include <mico/utilities.h>       // utility functions, in particular delete_image()
#include <mico/cross-sectional.h> // segmentation and bias correction of 3D image
#include <mico/longitudinal.h>    // temporally consistent segmentation of longitudinal series


#endif // _MICO_H
