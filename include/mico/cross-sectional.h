/**
 * @file  cross-sectional.h
 * @brief Implements MICO for cross-sectional MR brain studies.
 *
 * Copyright (c) 2011, 2012 University of Pennsylvania. All rights reserved.<br />
 * See https://www.med.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <software at cbica.upenn.edu>
 */

#pragma once
#ifndef _MICO_CROSS_SECTIONAL_H
#define _MICO_CROSS_SECTIONAL_H


#include <mico/types.h>


namespace mico {


// ===========================================================================
// API
// ===========================================================================

/**
 * @brief Segments and/or bias corrects a MR brain image.
 *
 * @param [in]  image       Input image.
 * @param [out] memberships Images of membership functions. Returns
 *                          one membership function for each considered
 *                          tissue class, CSF, GM, and WM, accessible
 *                          by these class names (case sensitive).
 *                          The returned images are allocated by this
 *                          function and have to be released by the
 *                          caller using the delete_image() function.
 *                          If NULL is given as argument, the membership
 *                          functions are not returned.
 * @param [out] bcimage     Bias corrected image with same datatype as
 *                          input image. The image is allocated by this
 *                          function and has to be destroyed by the caller
 *                          using delete_image(). If NULL is given as argument,
 *                          no bias corrected image is returned.
 * @param [in]  params      Parameters.
 * @param [in]  verbosity   Verbosity of output messages.
 *
 * @returns Whether the execution was successful.
 *
 * @see segment()
 * @see bias_correct()
 * @see create_label_map()
 *
 * @ingroup PublicAPI
 */
bool segment_and_bias_correct(Image*            image,
                              ImageMap*         memberships,
                              Image**           bcimage,
                              const Parameters& params    = Parameters(),
                              int               verbosity = 0);

/**
 * @brief Create a label image from the given membership functions.
 *
 * @param [in] image       Input image.
 * @param [in] memberships Images of membership functions as generated
 *                         by either segment_and_bias_correct() or
 *                         segment_longitudinal().
 * @param [in] params      Parameters, which in particular define the
 *                         labels to use for each tissue class.
 * @param [in] verbosity   Verbosity of output messages.
 *
 * @returns Segmentation image or NULL on failure. The returned image has to
 *          be destroyed by the caller using delete_image().
 *
 * @sa segment_and_bias_correct()
 * @sa segment_longitudinal()
 *
 * @ingroup PublicAPI
 */
Image* create_label_map(const Image*      image,
                        const ImageMap&   memberships,
                        const Parameters& params    = Parameters(),
                        int               verbosity = 0);

/**
 * @brief Segments a MR brain image.
 *
 * @param [in]  image     Input image.
 * @param [in]  params    Segmentation parameters.
 * @param [in]  verbosity Verbosity of output messages.
 *
 * @returns Segmentation image or NULL on failure. The returned image has to
 *          be destroyed by the caller using delete_image().
 *
 * @see segment_and_bias_correct()
 * @see create_label_map()
 *
 * @ingroup PublicAPI
 */
Image* segment(Image*            image,
               const Parameters& params    = Parameters(),
               int               verbosity = 0);

/**
 * @brief Bias corrects a MR brain image.
 *
 * @param [in]  image     Input image.
 * @param [in]  params    Bias correction parameters.
 * @param [in]  verbosity Verbosity of output messages.
 *
 * @returns Bias corrected image with same datatype as input image or NULL on
 *          failure. The returned image has to be destroyed by the caller
 *          using delete_image().
 *
 * @see segment_and_bias_correct()
 *
 * @ingroup PublicAPI
 */
Image* bias_correct(Image*            image,
                    const Parameters& params    = Parameters(),
                    int               verbosity = 0);

// ===========================================================================
// orthogonal basis function
// ===========================================================================

/**
 * @brief Get orthogonal basis function.
 *
 * @param [in] width  The width of input 3D image.
 * @param [in] height The height of input 3D image.
 * @param [in] slices The slices number of input 3D image.
 *
 * @returns Pointer to matrix of 20 orthogonal basis functions;
 *          The returned memory has to be released by the caller using delete [].
 */
float* get_basis_order(int width, int height, int slices);

// ===========================================================================
// update functions
// ===========================================================================

/**
 * @brief Iteration of C M b of 3D segmentation.
 *
 * @param [in]     img    Pointer to data of image  
 * @param [in]     ROI    Pointer to ROI forground and background of image
 * @param [in,out] M      Pointer to initialization of M
 * @param [in,out] C      Pointer to initialization of C
 * @param [in,out] b      Pointer to initialization of b
 * @param [in]     bias   Pointer to 20 orthogonal basis functions
 * @param [in]     Iter   Number of iteration
 * @param [in]     total  Number of image voxels.
 * @param [in]     params Parameters.    
 */
void updateMCB(float *img, float *ROI, float *M, float *C, float *b, float *bias, int total, const Parameters &params);

/**
 * @brief Update M:membership function.
 *
 * @param [in]     Img    Pointer to image
 * @param [in]     W      Pointer to ROI forground and background of image
 * @param [in,out] M      Pointer to membership function
 * @param [in]     C      Pointer to a centre of CSF WM GM
 * @param [in]     b      Pointer to bias filed
 * @param [in]     total  Number of image voxels.
 * @param [in]     params Parameters.   
 *  
 * @returns Pointer to matrix of M;
 *          The returned memory has to be released by the caller using delete [].
 */
float *updateM(float *Img, float *ROI, float *M, float *C, float *b, int total, const Parameters &params);

/**
 * @brief Update C:a centre of CSF, WM, and GM.
 *
 * @param [in]     Img    Pointer to data of image
 * @param [in]     ROI    Pointer to ROI forground and background of image   
 * @param [in]     M      Pointer to membership function
 * @param [in,out] C      Pointer to a centre of CSF WM GM. 
 * @param [in]     b      Pointer to bias filed
 * @param [in]     total  Number of image voxels.
 * @param [in]     params Parameters. 
 *  
 * @returns Pointer to matrix of C;
 *          The returned memory has to be released by the caller using delete [].
 */
float *updateC(float *Img, float *ROI, float *M, float *C, float *b, int total, const Parameters &params);

/**
 * @brief Update b:bias filed.
 *
 * @param [in]     img    Pointer to data of img  
 * @param [in]     ROI    Pointer to ROI forground and background of image   
 * @param [in]     M      Pointer to membership function
 * @param [in]     C      Pointer to a centre of CSF WM GM. 
 * @param [in,out] b      Pointer to bias filed
 * @param [in]     bias   Pointer to 20 orthogonal basis functions
 * @param [in]     total  Number of image voxels.
 * @param [in]     params Parameters.
 *  
 * @returns Pointer to matrix of b;
 *          The returned memory has to be released by the caller using delete [].
 */
float *updateB(float *img, float *ROI, float *M, float *C, float *b, float *bias, int total, const Parameters &params);

/**
 * @brief Segment image for update b C M.
 *
 * @param [in]  roi    Image region-of-interest.
 * @param [in]  bias   Pointer to 20 orthogonal basis functions.
 * @param [in]  p      Pointer to image data.
 * @param [in]  width  The width of input 3D image
 * @param [in]  height The height of input 3D image
 * @param [in]  slices The slices number of input 3D image
 * @param [in]  params Parameters. 
 * @param [out] result The resulting values of the functions M, C, and b.
 */
void updateBCM(const ImageRegion &roi, float *bias, float *p, int width, int height, int slices, const Parameters &params, float *result);


} // namespace mico


#endif // _MICO_CROSS_SECTIONAL_H
