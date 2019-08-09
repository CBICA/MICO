/**
 * @file  longitudinal.h
 * @brief Implements MICO for longitudinal MR brain image series.
 *
 * Copyright (c) 2011, 2012 University of Pennsylvania. All rights reserved.<br />
 * See https://www.med.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <software at cbica.upenn.edu>
 */

#pragma once
#ifndef _MICO_LONGITUDINAL_H_
#define _MICO_LONGITUDINAL_H_


#include <mico/types.h>


namespace mico {


// ===========================================================================
// API
// ===========================================================================

/**
 * @brief Segments a longitudinal image series and ensures temporal consistency.
 *
 * This function performs a 4D segmentation of a longitudinal image series
 * using a temporal consistency constraint. It is an extension of the temporally
 * unconstraint 3D segmentation.
 *
 * @ingroup PublicAPI
 */
bool segment_longitudinal(const ImageVector& images,
                          ImageMapVector&    membership,
                          const Parameters&  params    = Parameters(),
                          int                verbosity = 0);

// ===========================================================================
// temporally consistency
// ===========================================================================

/**
 * @brief Create temporal consistent membership function.
 *
 * This function calculates temporal consistent membership functions from the
 * unconstraint membership functions as obtained by the cross-sectional
 * segmentation.
 *
 * @param [in] images Data of image series.
 * @param [in] result Result of image segmentation (C, b, and M).
 * @param [in] params Parameters. 
 *
 * @returns Temporal consistent membership functions or an empty vector if
 *          not enough (temporary) memory was available.
 */
ImageMapVector tc_update_membership(const ImageVector& images,
                                    float**            result,
                                    const Parameters&  params);


} // namespace mico


#endif // _MICO_LONGITUDINAL_H
