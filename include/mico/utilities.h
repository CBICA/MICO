/**
 * @file  utilities.h
 * @brief Utility functions for MICO.
 *
 * Copyright (c) 2011, 2012 University of Pennsylvania. All rights reserved.<br />
 * See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#pragma once
#ifndef _MICO_UTILITIES_H
#define _MICO_UTILITIES_H

#include <iostream>

#include <nifti1_io.h>

#include <mico/types.h>


namespace mico {


// ===========================================================================
// sorting
// ===========================================================================

/**
 * @brief Sort C in-place. 
 *
 * @param [in,out] C Pointer to C--the center of SCF, GW and WM.
 */
void sort(float* C, int line);

// ===========================================================================
// image (de-)allocation
// ===========================================================================

/**
 * @brief Allocate three-dimensional image.
 *
 * @param [in] width    Width of the slices.
 * @param [in] height   Height of the slizes.
 * @param [in] slices   Number of slices.
 * @param [in] datatype Type of image data.
 *
 * @returns Allocated image of given size. Delete the image using delete_image().
 *
 * @sa delete_image()
 */
Image* new_image(int width, int height, int slices, short datatype);

/**
 * @brief Delete image.
 *
 * If the image does not own the image data, i.e.,image.ownsData is zero,
 * the image data is not managed by the given image. In this case, this
 * function does not delete the image data.
 *
 * Example:
 * @code
 * Image* image = new_image(128, 128, 64);
 * delete_image(image);
 * image = NULL;
 * @endcode
 *
 * @param [in] image Image. If NULL is given as argument, nothing is done.
 *
 * @sa newImage()
 */
void delete_image(Image* image);

/**
 * @brief Delete images.
 *
 * This function calls delete_image() for each element of the vector and then
 * clears the vector, i.e., setting its size to zero.
 *
 * @param [in] images Vector of images.
 *
 * @sa delete_image()
 */
void delete_images(ImageVector& images);

/**
 * @brief Delete images.
 *
 * This function calls delete_image() for each element of the map and then
 * clears the map, i.e., setting its size to zero.
 *
 * @param [in] images Map of images.
 *
 * @sa delete_image()
 */
void delete_images(ImageMap& images);

/**
 * @brief Delete images.
 *
 * This function calls delete_images() for each map of the vector and then
 * clears the vector, i.e., setting its size to zero.
 *
 * @param [in] images Vector of maps of images.
 *
 * @sa delete_images()
 */
void delete_images(ImageMapVector& images);

// ===========================================================================
// image conversion
// ===========================================================================

/**
 * @brief Whether a given NIfTI image is supported and can be converted to
 *        the internal Image type or not.
 *
 * @param [in] image NIfTI image.
 *
 * @retval  0 NIfTI image is supported and can be converted to own Image type.
 * @retval -1 Input image is not three-dimensional.
 * @retval -2 Input image is not scalar.
 * @retval -3 Invalid datatype.
 */
int is_supported_image(nifti_image* image);

/**
 * @brief Convert NIfTI image to own Image type.
 *
 * @param [in] image    NIfTI image.
 * @param [in] datatype Datatype of converted image. If negative, the
 *                      datatype of the NIfTI image is preserved.
 * @param [in] delnim   Whether to free the NIfTI image right after the
 *                      conversion. If false, the NIfTI image is not deleted
 *                      by this function.
 *
 * @returns Converted image or NULL if conversion failed.
 */
Image* convert_nifti_image(nifti_image* image, short datatype = -1, bool delnim = true);

/**
 * @brief Convert own image type to NIfTI.
 *
 * @attention If the image data is not copied, the NIfTI image will take
 *            over the ownership of the data. Hence, when the NIfTI image
 *            is deleted via nifti_image_free(), the image data is deleted
 *            as well. Therefore, image->owns_img is set to false by this
 *            function in this case.
 *
 * @param [in] image Image.
 * @param [in] fname Image file name prefix. If NULL,the filename entries
 *                   of the nifti_image are not set and remain NULL.  
 * @param [in] copy  Whether to copy the image data. If false, the NIfTI
 *                   image takes over the ownership of the image data.
 *
 * @returns NIfTI image or NULL on failure.
 */
nifti_image* convert_to_nifti(Image* image, const char* fname, bool copy = true);

/**
 * @brief Cast image to another data type.
 *
 * This function allocates a new image of the requested data type and copies
 * the image data from the given image to the new image, casting the image
 * values as necessary.
 *
 * If @p scale is true, the scaling function defined by the scl_slope and
 * scl_inter values of the NIfTI-1 header is applied to the image intensities
 * before casting the values to the output data type. Additionally, if the
 * output data type is an integral data type, the values are rescaled to the
 * range of the output data type before the cast, i.e., to [0, 255] in this
 * example. This may still result in a loss of image information, but the
 * degradation is only caused by the change of output data type and not the
 * scaling of intensities and may be further neglectable depending on the
 * application. In case of floating point data type as output data type,
 * no rescaling is performed. By default, a rescaling is always performed
 * whenever the scaling function of the NIfTI-1 header is applied.
 *
 * The scl_slope and scl_inter values of the NIfTI-1 header of the output
 * image are recalculated to reflect any scaling or rescaling that has be
 * applied to the intensity values stored in memory. Hence, by applying this
 * scaling function, one can recover the values of the input image which
 * correspond to the intensities which result from scaling the input values
 * according to the scaling function of the NIfTI-1 header of the input image.
 *
 * @param [in] image    Input image.
 * @param [in] datatype NIfTI-1 datatype code of casted image.
 * @param [in] scale    Whether to scale the intensities according to the
 *                      scaling function defined by the scl_slope and scl_inter
 *                      values of the NIfTI-1 header of the input image.
 *                      Furthermore, if this parameter is true, the intensities
 *                      are rescaled to the range of the output datatype in
 *                      case of integral valued types.
 *
 * @returns New image of given data type with optionally rescaled intensities
 *          or NULL if memory allocation failed. The returned image has to be
 *          destroyed by the caller using the delete_image() function.
 */
Image* cast_image(const Image* image, short datatype, bool scale = true);

/**
 * @brief Cast image to another data type.
 *
 * This function casts an image to the given output data type, applies the
 * scaling function as given by the NIfTI-1 header before, and further
 * rescales the intensities to the range specified by [@p min, @p max].
 * If either @p min or @p max are outside the range of values representable
 * by the chosen output data type, the image values which are mapped to
 * values outside the range of the data type are truncated.
 *
 * @param [in] image    Input image.
 * @param [in] datatype NIfTI-1 datatype code of casted image.
 * @param [in] min      Minimum of output intensity range.
 * @param [in] max      Maximum of output intensity range.
 *
 * @returns New image of given data type with a minimum intensity of @p min
 *          or the minimal value representable by the output data type,
 *          respectively, and a maximum intensity of @p max or the maximum
 *          value representable by the output data type, respectively, or
 *          NULL if memory allocation failed. The returned image has to be
 *          destroyed by the caller using delete_image().
 */
Image* cast_image(const Image* image, short datatype, float min, float max);

/**
 * @brief Copy image.
 *
 * @returns New image copy or NULL if memory allocation failed. The returned
 *          image has to be destroyed by the caller using delete_image().
 */
Image* copy_image(const Image* image);

// ===========================================================================
// image statistics
// ===========================================================================

/**
 * @brief Get intensity range in scalar image.
 *
 * @param [in]  image Image.
 * @param [out] min   Minimum intensity.
 * @param [out] max   Maximum intensity.
 * @param [in]  scale Whether to return intensity range of scaled intensities,
 *                    i.e., after applying scl_slope and scl_inter of image header.
 */
void get_intensity_range(const Image* image, float& min, float& max, bool scale = false);

// ===========================================================================
// intensity scaling
// ===========================================================================

/**
 * @brief Scale image intensities.
 *
 * This function applies the scaling function as given by the NIfTI-1 header
 * and, if @p rescale is true, further rescales the intensities to the maximum
 * range of the output datatype in case of non-floating point datatypes.
 *
 * @param [in, out] image   Image whose intensities shall be rescaled.
 * @param [in]      rescale Whether to rescale the intensities to the
 *                          maximum range of the image datatype for
 *                          non-floating point datatypes.
 *
 * @sa cast_image(const Image*, short, bool)
 */
void scale_image(Image* image, bool rescale = true);

/**
 * @brief Scale image intensities.
 *
 * This function applies the scaling function as given by the NIfTI-1 header
 * and further rescales the intensities to the range specified by [@p min, @p max].
 * If either @p min or @p max are outside the range of values representable
 * by the output datatype of the output image, the image values which are mapped
 * to values outside the range of the datatype are truncated.
 *
 * @param [in, out] image Image whose intensities shall be rescaled.
 * @param [in]      min   Minimum of output intensity range.
 * @param [in]      max   Maximum of output intensity range.
 *
 * @sa cast_image(const Image*, short, float, float)
 */
void scale_image(Image* image, float min, float max);

// ===========================================================================
// image regions
// ===========================================================================

/**
 * @brief Determine the foreground region.
 *
 * @param [in] image     Image.
 * @param [in] threshold Upper threshold of background values.
 *
 * @returns Region of foreground.
 */
ImageRegion get_foreground_region(Image* image, float threshold);

/**
 * @brief Join image regions.
 *
 * @param [in] region1 First image region.
 * @param [in] region2 Second image region.
 *
 * @returns The union of both regions.
 */
ImageRegion join_regions(ImageRegion region1, ImageRegion region2);

// ===========================================================================
// resize image
// ===========================================================================

/**
 * @brief Resize the image data.
 *
 * Resizing the image data corresponds to cropping parts of the image and
 * padding other parts of the image with a given value such that the
 * resulting image is within the bounds specified by the header information.
 *
 * @note The header information of the image is unchanged. Only the image
 *       data is reorganized and the image region information updated.
 *
 * @note The image must be the owner of the image data.
 *
 * @note If the size of the image data corresponds already to the specified
 *       image region, nothing is done.
 *
 * @param [in,out] image    Image whose data is resized.
 * @param [in]     region   The desired image region within the image bounds
 *                          specified by the header information.
 * @param [in]     padValue The value used to pad the image if necessary.
 *
 * @retval  0 On success.
 * @retval -1 If the image is not the owner of the image data.
 * @retval -2 If the image datatype is not supported by this function.
 * @retval -3 Memory allocation failed.
 */
int resize(Image* image, ImageRegion region, float padValue);

/**
 * @brief Crop the image data.
 *
 * This function crops a possibly previously already cropped image. Cropping
 * means, that this operation will only shrink the image but not expand it.
 * Therefore, the given image region is interpreted relative to the current
 * region of the image, i.e., region must be a subregion of image->region.
 *
 * @param [in,out] image  The image whose data should be cropped.
 * @param [in]     region The image region of the data that is kept.
 *
 * @retval  0 On success.
 * @retval -1 If the image is not the owner of the image data.
 * @retval -2 If the image datatype is not supported by this function.
 * @retval -3 Memory allocation failed.
 *
 * @sa resize()
 */
int crop(Image* image, ImageRegion region);

/**
 * @brief Undo the cropping of the image data.
 *
 * This function pads the image data such that the resulting image region has
 * the size of the original image as specified by the header information.
 * This function is useful to restore the original image size right before the
 * image data is written to disk.
 *
 * @param [in,out] image The image whose data should be padded.
 * @param [in]     value The value used for padding the image data, i.e.,
 *                       the values of voxels outside the current image
 *                       region are set to this value.
 *
 * @retval  0 On success.
 * @retval -1 If the image is not the owner of the image data.
 * @retval -2 If the image datatype is not supported by this function.
 * @retval -3 Memory allocation failed.
 *
 * @sa resize()
 */
int pad(Image *image, float value);

//////////////////////////////////////////////////////////////////////////////
// matrix operations
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief Invert nxn matrix.
 *
 * @param [in]  matrix  Input matrix.
 * @param [out] inverse Inverse of matrix.
 * @param [in]  n       Matrix size n, i.e., number of rows/columns.
 */
void invert_matrix(float* matrix, float* inverse, int n);

//////////////////////////////////////////////////////////////////////////////
// convolution
//////////////////////////////////////////////////////////////////////////////

/**
 * @brief Convolution of K and r.
 *
 * @param [in] K  First function.
 * @param [in] r  Second function.
 * @param [in] lk Number of samples of first function.
 * @param [in] lr Number of samples of second function.
 *  
 * @returns The result of the convolution of K by r. The number of output
 *          samples corresponds to @p lk. The returned memory has to be
 *          released by the caller using free().
 */
float* convolve(float* K, float* r, int lk, int lr);


} // namespace mico


#endif // _MICO_UTILITIES_H
