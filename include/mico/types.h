/**
 * @file  types.h
 * @brief Definition of data types.
 *
 * Copyright (c) 2011, 2012 University of Pennsylvania. All rights reserved.<br />
 * See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 *
 * @ingroup PublicAPI
 */

#pragma once
#ifndef _MICO_TYPES_H
#define _MICO_TYPES_H


#include <map>
#include <vector>
#include <string>
#include <nifti1.h>


namespace mico {


/**
 * @struct ImageRegion
 * @brief  Represents an image region in three dimensions.
 *
 * An image region is defined by an offset in voxel units and the size of the
 * image region. It is always relative to the entire image volume as specified
 * by the image's header information.
 *
 * Example:
 * @code
 * ImageRegion region;
 * @endcode
 *
 * Constraints on image region:
 * @code
 * 0 <= region.ox <= image.hdr.nx
 * 0 <= region.oy <= image.hdr.ny
 * 0 <= region.oz <= image.hdr.nz
 * @endcode
 *
 * @note If these constraints on the image region are not met,
 *       the image region is considered invalid. Functions
 *       can return an invalid region to indicate an error.
 *
 * Image region of entire image:
 * @code
 * region.ox = 0;
 * region.oy = 0;
 * region.oz = 0;
 * region.nx = image.hdr.nx;
 * region.ny = image.hdr.ny;
 * region.nz = image.hdr.nz;
 * @endcode
 */
struct ImageRegion
{
    int ox; ///< Offset of image region in x direction.
    int oy; ///< Offset of image region in y direction.
    int oz; ///< Offset of image region in z direction.
    int nx; ///< Size of image region in x direction.
    int ny; ///< Size of image region in y direction.
    int nz; ///< Size of image region in z direction.

    ImageRegion() : ox(0), oy(0), oz(0), nx(0), ny(0), nz(0) {}
};

/**
 * @union VoxelData
 * @brief Type used for image data pointer to enable the representation of the
 *        image data in either one of the supported NIfTI-1 data types.
 *
 * @sa Image
 */
union VoxelData
{
    unsigned char* uc; // datatype as used for labels
    short*         ss; // datatype as used for CT and MR images
    float*         fl; // datatype as used for calculations

    VoxelData()                 : uc(NULL) {}
    VoxelData(unsigned char* p) : uc(p) {}
    VoxelData(short*         p) : ss(p) {}
    VoxelData(float*         p) : fl(p) {}
};

/**
 * @struct Image
 * @brief  Data structure used to represent an image.
 *
 * Example:
 * @code
 * Image image;
 * @endcode
 *
 * Voxel positions in the image are given by three indices:
 * @code
 * int c; // column index
 * int r; // row index
 * int s; // slice index
 * @endcode
 *
 * How to access image data of voxel with coordinates [c, r, s]:
 * @code
 * float value = image.img.fl [s * image.hdr.nx * image.hdr.ny + r * image.hdr.nx + c];
 * @endcode
 *
 * How to efficiently iterate over entire image data
 * (note that the image data may have been cropped):
 * @code
 * float *value = image.img.fl;
 * 
 * for (s = 0; s < image.region.nz; ++ s)
 * {
 *     for (r = 0; r < image.region.ny; ++ r)
 *     {
 *         for (c = 0; c < image.region.nx; ++ c, ++ value)
 *         {
 *             float f = *value;
 *         }
 *     }
 * }
 * @endcode
 *
 * @sa newImage()
 * @sa deleteImage()
 *
 * @ingroup PublicAPI
 */
struct Image
{
    std::string     name;       ///< Image name.
    nifti_1_header  hdr;        ///< NIfTI-1 header of image.
    int             nifti_type; ///< Type of image when read from file.
    bool            compress;   ///< Whether to write image data to compressed image file.
    ImageRegion     region;     ///< Region of image data within image extent specified by NIfTI-1 header.
    VoxelData       img;        ///< Image data stored in continuous memory.
    bool            owns_img;   ///< Whether this image owns the image data.

    /**
     * @brief Default constructor.
     */             
    Image()
    {
        nifti_type = 1; // NIFTI_FTYPE_NIFTI1_1;
        compress   = true;
        owns_img   = true;
    }
};

/// @brief Type used to store a vector of image pointers.
typedef std::vector<Image*> ImageVector;
/// @brief Non-const iterator over elements of ImageVector.
typedef std::vector<Image*>::iterator ImageIterator;
/// @brief Const iterator over elements of ImageVector.
typedef std::vector<Image*>::const_iterator ImageConstIterator;
/// @brief Type used to map class names to image pointers.
typedef std::map<const char*, Image*> ImageMap;
/// @brief Type used to store multiple maps of class names to image pointers.
typedef std::vector<ImageMap> ImageMapVector;

/**
 * @struct Parameters.
 * @brief  Structure of parameters used for main functions.
 *
 * @ingroup PublicAPI
 */
struct Parameters
{
    bool           outputSegments;        ///< Whether to output the segmentation images.
    bool           outputBiasCorrections; ///< Whether to output the bias corrected images.
    bool           outputMemberships;     ///< Whether to output the membership images.
    std::string    outputDir;             ///< Output directory.
    std::string    suffix;                ///< Suffix used for output images.
    std::string    suffixBC;              ///< Suffix used for output bias corrected images.
    std::string    suffixMF;              ///< Suffix used for output membership images.
    unsigned char  labelCSF;              ///< Label used for CSF.
    unsigned char  labelGM;               ///< Label used for GM.
    unsigned char  labelWM;               ///< Label used for WM.
    float          q;                     ///< @brief Fuzzifier, must be >= 1.
                                          ///<
                                          ///< When q = 1, the algorithm performs hard segmentation,
                                          ///< resulting in binary membership functions
                                          ///< (with value 0 or 1) as the segmentation result.
                                          ///< When q > 1, the algorithm performs soft segmentation,
                                          ///< resulting in fuzzy membership functions with value
                                          ///< between 0 and 1. For soft segmentation, the default
                                          ///< value of q is 1.5.
    float          lambda;                ///< @brief Degree of bias correction.
                                          ///< 
                                          ///< A small value should be used for images with significant
                                          ///< inhomogenity. A value of about 10 for images
                                          ///< with no inhomogenity.
    int            iterMCB;               ///< Number of iterations of outer loop in updateMCB().
    int            iterMC;                ///< Number of iterations of inner loop in updateMCB().
    int            iterNum;               ///< Number of iterations of updateMCB() in updateBCM().
    float          th_bg;                 ///< Background threshold used to determine foreground ROI.
    float          tissueWeight[3];       ///< Weight of CSF, GM, and WM (in this order).
    unsigned int   seed;                  ///< Seed used for pseudo-random number generator (@sa srand()).

    /**
     * @brief Default constructor.
     * 
     * Sets default values for parameters.
     */                   
    Parameters()
    {
        outputSegments        = true;
        outputBiasCorrections = false;
        outputMemberships     = false;
        labelCSF              = 10;
        labelGM               = 150;
        labelWM               = 250;
        q                     = 1.5;
        lambda                = 1.0e-4;
        iterMCB               = 1;
        iterMC                = 5;
        iterNum               = 10;
        th_bg                 = 0.5;
        tissueWeight[0]       = 0.8;
        tissueWeight[1]       = 1;
        tissueWeight[2]       = 1.2;
        seed                  = 1;
    }
};


} // namespace mico


#endif // _MICO_TYPES_H
