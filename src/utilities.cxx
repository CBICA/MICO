/**
 * @file  utilities.cxx
 * @brief Utility functions for MICO.
 *
 * Copyright (c) 2011, 2012 University of Pennsylvania. All rights reserved.<br />
 * See https://www.med.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <software at cbica.upenn.edu>
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <nifti1_io.h>
#include <limits>
#include <set>

// #include <basis/assert.h> // assert(), ASSERT()
// #include <basis/except.h> // BASIS_THROW(), std::invalid_argument
// #include <basis/os.h>     // hasext()

#include "mico/utilities.h"
#include "cbicaUtilities.h"

// acceptable in .cxx file
using namespace std;
//using namespace basis;


namespace mico {


// ===========================================================================
// private auxiliary functions
// ===========================================================================

/**
 * @brief Round floating-point value.
 */
inline float round(float v)
{
    return (v > 0.0f) ? floor(v + 0.5f) : ceil(v - 0.5f);
}

/**
 * @brief Get minimum and maximum image intensity.
 *
 * @param [in]  in     Image data.
 * @param [in]  x_size Image size in first dimension.
 * @param [in]  y_size Image size in second dimension.
 * @param [in]  z_size Image size in third dimension.
 * @param [out] min    Minimum image intensity.
 * @param [out] max    Maximum image intensity.
 */
template <typename T>
void get_min_max(const T* in, int x_size, int y_size, int z_size, float& min, float& max)
{
    T lb = numeric_limits<T>::max();
    T ub = numeric_limits<T>::min();
    for (int k = 0; k < z_size; k++) {
        for (int j = 0; j < y_size; j++) {
            for (int i = 0; i < x_size; i++, in++) {
                if (*in < lb) lb = *in;
                if (*in > ub) ub = *in;
            }
        }
    }
    if (lb > ub) {
        min = 0.0f;
        max = 0.0f;
    } else {
        min = static_cast<float>(lb);
        max = static_cast<float>(ub);
    }
}

// ---------------------------------------------------------------------------
void get_min_max(const float* in, int x_size, int y_size, int z_size, float& min, float& max)
{
    float lb = + numeric_limits<float>::max();
    float ub = - numeric_limits<float>::max();
    for (int k = 0; k < z_size; k++) {
        for (int j = 0; j < y_size; j++) {
            for (int i = 0; i < x_size; i++, in++) {
                if (*in < lb) lb = *in;
                if (*in > ub) ub = *in;
            }
        }
    }
    if (lb > ub) {
        min = 0.0f;
        max = 0.0f;
    } else {
        min = lb;
        max = ub;
    }
}

/**
 * @brief Adjust rescale slope and intercept.
 *
 * This function adjusts the slope and intercept of the given linear scaling
 * function such that the intensities of the rescaled image are rescaled to
 * the given maximum range of intensities which can be represented by the
 * output datatype.
 *
 * @param [in]      in     Input image data.
 * @param [in]      x_size Image size in first dimension.
 * @param [in]      y_size Image size in second dimension.
 * @param [in]      z_size Image size in third dimension.
 * @param [in]      lmin   Minimum representable output intensity.
 * @param [in]      lmax   Maximum representable output intensity.
 * @param [in, out] slope  Scaling slope.
 * @param [in, out] inter  Scaling intercept.
 */
template <typename T>
void adjust_scaling(T* const in, int x_size, int y_size, int z_size,
                    float lmin, float lmax,
                    float& slope, float& inter)
{
    // get intensity range of image data
    float min, max;
    get_min_max(in, x_size, y_size, z_size, min, max);
    // scaling function to map [min, max] to [lmin, lmax]
    float scl_slope = (lmax - lmin) / (max - min);
    float scl_inter = lmin - min * scl_slope;
    // adjust input scaling function
    slope = slope * scl_slope;
    inter = inter * scl_slope + scl_inter;
}

/**
 * @brief Rescale image intensities.
 *
 * This function rescales the intensities of the input image using the linear
 * function y = slope * x + inter and writes the rescaled intensities to the
 * given output image. If a value y is outside the range of values which can
 * be represented by the output data type, these values are truncated and
 * mapped to the closest representable value.
 *
 * @param [out] out    Rescaled image. Set to @c in to rescale image in-place.
 * @param [in]  in     Input image.
 * @param [in]  x_size Image size in first dimension.
 * @param [in]  y_size Image size in second dimension.
 * @param [in]  z_size Image size in third dimension.
 * @param [in]  slope  Rescaling slope.
 * @param [in]  inter  Rescaling intercept.
 */
template <typename TOut, typename TIn>
void rescale_intensities(TOut*      out,
                         const TIn* in,
                         int        x_size,
                         int        y_size,
                         int        z_size,
                         float      slope,
                         float      inter)
{
    const float min = static_cast<float>(numeric_limits<TOut>::min());
    const float max = static_cast<float>(numeric_limits<TOut>::max());
    float value;
    for (int k = 0; k < z_size; k++) {
        for (int j = 0; j < y_size; j++) {
            for (int i = 0; i < x_size; i++, in++, out++) {
                value = round(static_cast<float>(*in) * slope + inter);
                if      (value < min) value = min;
                else if (value > max) value = max;
                (*out) = static_cast<TOut>(value);
            }
        }
    }
}

// ---------------------------------------------------------------------------
template <typename TIn>
void rescale_intensities(float*     out,
                         const TIn* in,
                         int        x_size,
                         int        y_size,
                         int        z_size,
                         float      slope,
                         float      inter)
{
    for (int k = 0; k < z_size; k++) {
        for (int j = 0; j < y_size; j++) {
            for (int i = 0; i < x_size; i++, in++, out++) {
                (*out) = static_cast<float>(*in) * slope + inter;
            }
        }
    }
}

/**
 * @brief Copy image data and (optionally) scale intensities.
 *
 * @param [in]  image_out Output image.
 * @param [out] p         Output image data pointer.
 * @param [in]  image_in  Input image.
 * @param [in]  scale     Whether to scale the image intensities.
 * @param [in]  lmin      Minimum output intensity.
 * @param [in]  lmax      Maximum output intensity. Set to value less or equal
 *                        to @p lmin to only apply the scaling of the NIfTI-1
 *                        image header of the input image.
 */
template <typename T>
void copy_image_data(Image*             image_out,
                     T*                 p,
                     const Image* const image_in,
                     bool               scale,
                     float              lmin,
                     float              lmax)
{
    // determine rescaling slope and intercept
    float slope = 1.0f;
    float inter = 0.0f;
    if (scale && image_in->hdr.scl_slope != 0.0f) {
        slope = image_in->hdr.scl_slope;
        inter = image_in->hdr.scl_inter;
    }
    if (scale && (lmax > lmin) &&
        // do not rescale intensities to maximum range if datatype is DT_FLOAT and no
        // image range other than the maximum range is specified
            (image_out->hdr.datatype != DT_FLOAT ||
                - numeric_limits<float>::max() < lmin || lmax < numeric_limits<float>::max())) {
        switch (image_in->hdr.datatype) {
            case DT_UNSIGNED_CHAR:
                adjust_scaling(image_in->img.uc,
                               image_in->region.nx, image_in->region.ny, image_in->region.nz,
                               lmin, lmax, slope, inter);
                break;
            case DT_SIGNED_SHORT:
                adjust_scaling(image_in->img.ss,
                               image_in->region.nx, image_in->region.ny, image_in->region.nz,
                               lmin, lmax, slope, inter);
                break;
            case DT_FLOAT:
                adjust_scaling(image_in->img.fl,
                               image_in->region.nx, image_in->region.ny, image_in->region.nz,
                               lmin, lmax, slope, inter);
                break;
            default:
              std::cerr << "Invalid image datatype: " << image_in->hdr.datatype << "\n";
                abort();
        }
    }
    // copy and (optionally) rescale image data
    switch (image_in->hdr.datatype) {
        case DT_UNSIGNED_CHAR:
            rescale_intensities(p, image_in->img.uc,
                                image_in->region.nx, image_in->region.ny, image_in->region.nz,
                                slope, inter);
            break;
        case DT_SIGNED_SHORT:
            rescale_intensities(p, image_in->img.ss,
                                image_in->region.nx, image_in->region.ny, image_in->region.nz,
                                slope, inter);
            break;
        case DT_FLOAT:
            rescale_intensities(p, image_in->img.fl,
                                image_in->region.nx, image_in->region.ny, image_in->region.nz,
                                slope, inter);
            break;
        default:
            std::cerr << "Invalid image datatype: " << image_in->hdr.datatype << "\n";
            abort();
    }
    // set scl_slope and scl_inter of output image
    // Attention: image_in and image_out may reference the same image!
    if (image_in->hdr.scl_slope == 0.0) {
        image_out->hdr.scl_slope = 1.0 / slope;
        image_out->hdr.scl_inter = - inter / slope;
    } else {
        image_out->hdr.scl_slope = image_in->hdr.scl_slope / slope;
        image_out->hdr.scl_inter = image_in->hdr.scl_inter - inter / slope;
    }
    // if slope is 1 and intercept is 0, set both to 0
    if (fabs(image_out->hdr.scl_slope - 1.0) < 1.0e-9 && fabs(image_out->hdr.scl_inter) < 1.0e-9) {
        image_out->hdr.scl_slope = 0.0;
        image_out->hdr.scl_inter = 0.0;
    }
}

// ===========================================================================
// sorting
// ===========================================================================

// ---------------------------------------------------------------------------
void sort(float* C, int line)
{
    float tmp;
 
    for (int j = 1; j < line; j++) {
        for(int i = 0; i < line - j; i++) { 
            if (C [i] > C [i + 1]) {
                tmp       = C [i];
                C [i    ] = C [i + 1];
                C [i + 1] = tmp;
            }
        }
    }
}

// ===========================================================================
// image (de-)allocation
// ===========================================================================

// ---------------------------------------------------------------------------
Image* new_image(int width, int height, int slices, short datatype)
{
    nifti_1_header* hdr   = NULL;
    int*            dims  = NULL;

    Image* image = new Image();
    if (image == NULL) return NULL;

    // allocate and initialize image data
    image->img.fl     = NULL;
    image->owns_img   = true;
    image->nifti_type = NIFTI_FTYPE_NIFTI1_1;
    image->compress   = true;

    switch (datatype) {
        case DT_UNSIGNED_CHAR:
            image->img.uc = reinterpret_cast<unsigned char*>(calloc(width * height * slices, sizeof(unsigned char)));
            break;
        case DT_SIGNED_SHORT:
            image->img.ss = reinterpret_cast<short*>(calloc(width * height * slices, sizeof(short)));
            break;
        case DT_FLOAT:
            image->img.fl = reinterpret_cast<float*>(calloc(width * height * slices, sizeof(float)));
            break;
        default:
            free(image);
            return NULL;
    };

    if (image->img.fl == NULL) {
        free(image);
        return NULL;
    }

    // image header
    dims = reinterpret_cast<int*>(malloc(8 * sizeof(int)));

    dims[0] = 5;
    dims[1] = width;
    dims[2] = height;
    dims[3] = slices;
    dims[4] = 1;
    dims[5] = 1;
    dims[6] = 0;
    dims[7] = 0;

    hdr = nifti_make_new_header(dims, datatype);
 
    free(dims);

    if (hdr == NULL) {
        free(image->img.fl);
        return NULL;
    }

    image->hdr = *hdr;
    free(hdr);

    // image region
    image->region.ox = 0;
    image->region.oy = 0;
    image->region.oz = 0;
    image->region.nx = width;
    image->region.ny = height;
    image->region.nz = slices;

    return image;
}

// ---------------------------------------------------------------------------
void delete_image(Image* image)
{
    if (image != NULL) {
        if (image->owns_img) free(image->img.fl);
        delete image;
    }
}

// ---------------------------------------------------------------------------
void delete_images(ImageVector& images)
{
    for (ImageIterator it = images.begin(); it != images.end(); ++it) {
        delete_image(*it);
    }
    images.clear();
}

// ---------------------------------------------------------------------------
void delete_images(ImageMap& images)
{
    for (ImageMap::iterator it = images.begin(); it != images.end(); ++it) {
        delete_image(it->second);
    }
    images.clear();
}

// ---------------------------------------------------------------------------
void delete_images(ImageMapVector& images)
{
    for (ImageMapVector::iterator it = images.begin(); it != images.end(); ++it) {
        delete_images(*it);
    }
    images.clear();
}

// ===========================================================================
// image conversion
// ===========================================================================

// ---------------------------------------------------------------------------
int is_supported_image(nifti_image* nim)
{
    // image must be three dimensional
    if (nim->nx <= 0 || nim->ny <= 0 || nim->nz <= 0) return -1;
    // image must be scalar
    if (nim->nt != 0 && nim->nt != 1) return -2;
    if (nim->nu != 0 && nim->nu != 1) return -2;
    if (nim->nv != 0 && nim->nv != 1) return -2;
    if (nim->nw != 0 && nim->nw != 1) return -2;
    // image data type must be either one fo the following
    if (nim->datatype != DT_UNSIGNED_CHAR &&
        nim->datatype != DT_SIGNED_SHORT  &&
        nim->datatype != DT_SIGNED_INT    &&
        nim->datatype != DT_FLOAT ) return -3;
    // image is supported, otherwise
    return 0;
}

// ---------------------------------------------------------------------------
Image* convert_nifti_image(nifti_image* nim, short datatype, bool delnim)
{
    Image* image = NULL; // converted image

    // verify that NIfTI image can be converted
    if (is_supported_image(nim) != 0) return NULL;

    // instantiate own image type which wraps the nifti_image data
    // note that input data can be anything else then DT_FLOAT as well
    set<string> niftiexts;
    niftiexts.insert(".nii");
    niftiexts.insert(".nii.gz");
    niftiexts.insert(".hdr");
    niftiexts.insert(".hdr.gz");
    niftiexts.insert(".img");
    niftiexts.insert(".img.gz");

    string fext;
    string fbase;
    string fpath;
    //os::path::splitext(nim->fname, fbase, fext, &niftiexts);
    cbica::splitFileName(nim->fname, fpath, fbase, fext);

    Image image_in;
    image_in.hdr        = nifti_convert_nim2nhdr(nim);
    image_in.region.ox  = 0;
    image_in.region.oy  = 0;
    image_in.region.oz  = 0;
    image_in.region.nx  = nim->nx;
    image_in.region.ny  = nim->ny;
    image_in.region.nz  = nim->nz;
    image_in.img.fl     = reinterpret_cast<float*>(nim->data);
    image_in.owns_img   = false;
    image_in.name       = fbase/* os::path::basename(fbase)*/;
    image_in.compress   = (fext == ".gz");
    image_in.nifti_type = nim->nifti_type;

    // instantiate own image type
    if (datatype < 0) datatype = image_in.hdr.datatype;
    image = new_image(nim->nx, nim->ny, nim->nz, datatype);
    if (image == NULL) return NULL;

    // copy header information
    image->hdr = image_in.hdr;

    // adjust datatype related fields
    image->hdr.datatype = datatype;
    if      (datatype == DT_UNSIGNED_CHAR) image->hdr.bitpix = 8 * sizeof(unsigned char);
    else if (datatype == DT_SIGNED_SHORT)  image->hdr.bitpix = 8 * sizeof(short);
    else if (datatype == DT_FLOAT)         image->hdr.bitpix = 8 * sizeof(float);
    else 
    {
      std::cerr << "Check implementation of new_image()! It should have returned NULL!\n";
      exit(EXIT_FAILURE);
    }

    // remember original nifti_type
    // (i.e., to know if input image was stored in ANALYZE format)
    image->nifti_type = image_in.nifti_type;

    // set name and compress attributes
    image->name     = image_in.name;
    image->compress = image_in.compress;

    // copy image data and apply scaling
    switch (datatype) {
        case DT_UNSIGNED_CHAR:
            copy_image_data(image, image->img.uc, &image_in, true,
                    numeric_limits<unsigned char>::min(),
                    numeric_limits<unsigned char>::max());
            break;
        case DT_SIGNED_SHORT:
            copy_image_data(image, image->img.ss, &image_in, true,
                    numeric_limits<short>::min(),
                    numeric_limits<short>::max());
            break;
        case DT_FLOAT:
            copy_image_data(image, image->img.fl, &image_in, true,
                    - numeric_limits<float>::max(),
                    + numeric_limits<float>::max());
            break;
    };

    // delete converted NIfTI image
    if (delnim) nifti_image_free(nim);

    return image;
}

// ---------------------------------------------------------------------------
nifti_image* convert_to_nifti(Image* image, const char* fname, bool copy)
{
    nifti_image* nim = NULL;
    if (!image) return NULL;
    // convert header to NIfTI image
    if (nifti_is_complete_filename(fname)) {
        nim = nifti_convert_nhdr2nim(image->hdr, const_cast<char*>(fname));
    } else {
        nim = nifti_convert_nhdr2nim(image->hdr, NULL);
        nim->fname = nifti_makehdrname(fname, nim->nifti_type, 0, image->compress);
        nim->iname = nifti_makeimgname(fname, nim->nifti_type, 0, image->compress);
        if (!nim->fname || !nim->iname) return NULL;
    }
    // force output of ANALYZE image if input image was stored in ANALYZE
    // format and fname was set to a valid ANALYZE extension
    set<string> analyzeexts; // note: can also be NIfTI extension, though
    analyzeexts.insert(".hdr");
    analyzeexts.insert(".hdr.gz");
    analyzeexts.insert(".img");
    analyzeexts.insert(".img.gz");
    set<string> niftiexts;
    niftiexts.insert(".nii");
    niftiexts.insert(".nii.gz");
    niftiexts.insert(".nia");
    niftiexts.insert(".nia.gz");
//    auto it = analyzeexts.find(cbica::getFilenameExtension(nim->fname));
//    if (it != analyzeexts.end())
//    {
//      std::cout << "'at' found in set" << std::endl;
//    }
//    if ((it != analyzeexts.end()) && image->nifti_type == NIFTI_FTYPE_ANALYZE) {
//        nim->nifti_type = NIFTI_FTYPE_ANALYZE;
//    }
    // set image data of NIfTI image
    if (copy) {
        size_t nbytes = image->hdr.dim[1] * image->hdr.dim[2] * image->hdr.dim[3] * (image->hdr.bitpix / 8);
        nim->data = malloc(nbytes);
        if (nim->data == NULL) {
            nifti_image_free(nim);
            nim = NULL;
        } else {
            // note: can be any data, also non-float data...
            memcpy(nim->data, reinterpret_cast<void*>(image->img.fl), nbytes);
        }
    } else {
        // note: can be any data, also non-float data...
        nim->data       = reinterpret_cast<void*>(image->img.fl);
        image->owns_img = false;
    }

    return nim;
}

// ---------------------------------------------------------------------------
Image* cast_image(const Image* image, short datatype, bool scale)
{
    // allocate output image
    Image* casted_image = new_image(image->region.nx,
                                    image->region.ny,
                                    image->region.nz, datatype);
    if (casted_image == NULL) return NULL;
    // copy image attributes
    int bitpix = casted_image->hdr.bitpix;
    casted_image->hdr = image->hdr;
    casted_image->hdr.bitpix   = bitpix;
    casted_image->hdr.datatype = datatype;
    casted_image->name         = image->name;
    casted_image->compress     = image->compress;
    casted_image->nifti_type   = image->nifti_type;
    // copy image data and (optionally) apply scaling
    switch (datatype) {
        case DT_UNSIGNED_CHAR:
            copy_image_data(casted_image, casted_image->img.uc, image, scale,
                    numeric_limits<unsigned char>::min(),
                    numeric_limits<unsigned char>::max());
            break;
        case DT_SIGNED_SHORT:
            copy_image_data(casted_image, casted_image->img.ss, image, scale,
                    numeric_limits<short>::min(),
                    numeric_limits<short>::max());
            break;
        case DT_FLOAT:
            copy_image_data(casted_image, casted_image->img.fl, image, scale,
                    - numeric_limits<float>::max(),
                    + numeric_limits<float>::max());
            break;
    };
    return casted_image;
}

// ---------------------------------------------------------------------------
Image* cast_image(const Image* image, short datatype, float min, float max)
{
    // allocate output image
    Image* casted_image = new_image(image->region.nx,
                                    image->region.ny,
                                    image->region.nz, datatype);
    if (casted_image == NULL) return NULL;
    // copy image attributes
    int bitpix = casted_image->hdr.bitpix;
    casted_image->hdr = image->hdr;
    casted_image->hdr.bitpix   = bitpix;
    casted_image->hdr.datatype = datatype;
    casted_image->name         = image->name;
    casted_image->compress     = image->compress;
    casted_image->nifti_type   = image->nifti_type;
    // copy image data and apply scaling
    switch (datatype) {
        case DT_UNSIGNED_CHAR:
            copy_image_data(casted_image, casted_image->img.uc, image, true, min, max);
            break;
        case DT_SIGNED_SHORT:
            copy_image_data(casted_image, casted_image->img.ss, image, true, min, max);
            break;
        case DT_FLOAT:
            copy_image_data(casted_image, casted_image->img.fl, image, true, min, max);
            break;
    }
    return casted_image;
}

// ---------------------------------------------------------------------------
Image* copy_image(const Image* image)
{
    return cast_image(image, image->hdr.datatype, false);
}

// ===========================================================================
// image statistics
// ===========================================================================

// ---------------------------------------------------------------------------
void get_intensity_range(const Image* image, float& min, float& max, bool scale)
{
    switch (image->hdr.datatype) {
        case DT_UNSIGNED_CHAR:
            get_min_max(image->img.uc,
                        image->region.nx,
                        image->region.ny,
                        image->region.nz,
                        min, max);
            break;
        case DT_SIGNED_SHORT:
            get_min_max(image->img.ss,
                        image->region.nx,
                        image->region.ny,
                        image->region.nz,
                        min, max);
            break;
        case DT_FLOAT:
            get_min_max(image->img.fl,
                        image->region.nx,
                        image->region.ny,
                        image->region.nz,
                        min, max);
            break;
        default:
          std::cout << "Invalid image datatype: " << image->hdr.datatype;
    }
    if (scale && image->hdr.scl_slope != 0) {
        min = min * image->hdr.scl_slope + image->hdr.scl_inter;
        max = max * image->hdr.scl_slope + image->hdr.scl_inter;
    }
}

// ===========================================================================
// image conversion
// ===========================================================================

// ---------------------------------------------------------------------------
void scale_image(Image* image, bool rescale)
{
    if (rescale) {
        switch (image->hdr.datatype) {
            case DT_UNSIGNED_CHAR:
                copy_image_data(image, image->img.uc, image, true,
                        numeric_limits<unsigned char>::min(),
                        numeric_limits<unsigned char>::max());
                break;
            case DT_SIGNED_SHORT:
                copy_image_data(image, image->img.ss, image, true,
                        numeric_limits<short>::min(),
                        numeric_limits<short>::max());
                break;
            case DT_FLOAT:
                copy_image_data(image, image->img.fl, image, true,
                        - numeric_limits<float>::max(),
                        + numeric_limits<float>::max());
                break;
        }
    } else {
        switch (image->hdr.datatype) {
            case DT_UNSIGNED_CHAR:
                copy_image_data(image, image->img.uc, image, true, 0, 0);
                break;
            case DT_SIGNED_SHORT:
                copy_image_data(image, image->img.ss, image, true, 0, 0);
                break;
            case DT_FLOAT:
                copy_image_data(image, image->img.fl, image, true, 0, 0);
                break;
        }
    }
}

// ---------------------------------------------------------------------------
void scale_image(Image* image, float min, float max)
{
    switch (image->hdr.datatype) {
        case DT_UNSIGNED_CHAR:
            copy_image_data(image, image->img.uc, image, true, min, max);
            break;
        case DT_SIGNED_SHORT:
            copy_image_data(image, image->img.ss, image, true, min, max);
            break;
        case DT_FLOAT:
            copy_image_data(image, image->img.fl, image, true, min, max);
            break;
    }
}

// ===========================================================================
// image regions
// ===========================================================================

// ---------------------------------------------------------------------------
ImageRegion get_foreground_region(Image* image, float threshold)
{
    ImageRegion region;  // return value
    int         x, y, z; // loop variables
    float       value;   // voxel value

    int minx = 0;
    int miny = 0;
    int minz = 0;
    int maxx = image->region.nx - 1;
    int maxy = image->region.ny - 1;
    int maxz = image->region.nz - 1;

    int sizeOfSlice = image->region.nx * image->region.ny;
    int datatype    = image->hdr.datatype;

    if (datatype != DT_UNSIGNED_CHAR && datatype != DT_SIGNED_SHORT && datatype != DT_FLOAT) {
        return region;
    }

    // minium z offset
    for (z = minz; z <= maxz; ++ z) {
        for (y = miny; y <= maxy; ++ y) {
            for (x = minx; x <= maxx; ++ x) {
                if (datatype == DT_UNSIGNED_CHAR) {
                    value = static_cast<float>(image->img.uc[z * sizeOfSlice + y * image->region.nx + x]);
                } else if (datatype == DT_SIGNED_SHORT) {
                    value = static_cast<float>(image->img.ss[z * sizeOfSlice + y * image->region.nx + x]);
                } else {
                    value = image->img.fl[z * sizeOfSlice + y * image->region.nx + x];
                }
                if (value > threshold) {
                    minz = z;
                    // break out of ALL loops
                    y = maxy + 1;
                    z = maxz + 1;
                    break;
                }
            }
        }
    }

    // maximum z offset
    for (z = maxz; z >= minz; -- z) {
        for (y = miny; y <= maxy; ++ y) {
            for (x = minx; x <= maxx; ++ x) {
                if (datatype == DT_UNSIGNED_CHAR) {
                    value = static_cast<float>(image->img.uc[z * sizeOfSlice + y * image->region.nx + x]);
                } else if (datatype == DT_SIGNED_SHORT) {
                    value = static_cast<float>(image->img.ss[z * sizeOfSlice + y * image->region.nx + x]);
                } else {
                    value = image->img.fl[z * sizeOfSlice + y * image->region.nx + x];
                }
                if (value > threshold) {
                    maxz = z;
                    // break out of ALL loops
                    y = maxy + 1;
                    z = minz - 1;
                    break;
                }
            }
        }
    }
 
    // minium y offset
    for (y = miny; y <= maxy; ++ y) {
        for (z = minz; z <= maxz; ++ z) {
            for (x = minx; x <= maxx; ++ x) {
                if (datatype == DT_UNSIGNED_CHAR) {
                    value = static_cast<float>(image->img.uc[z * sizeOfSlice + y * image->region.nx + x]);
                } else if (datatype == DT_SIGNED_SHORT) {
                    value = static_cast<float>(image->img.ss[z * sizeOfSlice + y * image->region.nx + x]);
                } else {
                    value = image->img.fl[z * sizeOfSlice + y * image->region.nx + x];
                }
                if (value > threshold) {
                    miny = y;
                    // break out of ALL loops
                    z = maxz + 1;
                    y = maxy + 1;
                    break;
                }
            }
        }
    }

    // maximum y offset
    for (y = maxy; y >= miny; -- y) {
        for (z = minz; z <= maxz; ++ z) {
            for (x = minx; x <= maxx; ++ x) {
                if (datatype == DT_UNSIGNED_CHAR) {
                    value = static_cast<float>(image->img.uc[z * sizeOfSlice + y * image->region.nx + x]);
                } else if (datatype == DT_SIGNED_SHORT) {
                    value = static_cast<float>(image->img.ss[z * sizeOfSlice + y * image->region.nx + x]);
                } else {
                    value = image->img.fl[z * sizeOfSlice + y * image->region.nx + x];
                }
                if (value > threshold) {
                    maxy = y;
                    // break out of ALL loops
                    z = maxz + 1;
                    y = miny - 1;
                    break;
                }
            }
        }
    }

    // minium x offset
    for (x = minx; x <= maxx; ++ x) {
        for (z = minz; z <= maxz; ++ z) {
            for (y = miny; y <= maxy; ++ y) {
                if (datatype == DT_UNSIGNED_CHAR) {
                    value = static_cast<float>(image->img.uc[z * sizeOfSlice + y * image->region.nx + x]);
                } else if (datatype == DT_SIGNED_SHORT) {
                    value = static_cast<float>(image->img.ss[z * sizeOfSlice + y * image->region.nx + x]);
                } else {
                    value = image->img.fl[z * sizeOfSlice + y * image->region.nx + x];
                }
                if (value > threshold) {
                    minx = x;
                    // break out of ALL loops
                    z = maxz + 1;
                    x = maxx + 1;
                    break;
                }
            }
        }
    }

    // maximum x offset
    for (x = maxx; x >= minx; -- x) {
        for (z = minz; z <= maxz; ++ z) {
            for (y = miny; y <= maxy; ++ y) {
                if (datatype == DT_UNSIGNED_CHAR) {
                    value = static_cast<float>(image->img.uc[z * sizeOfSlice + y * image->region.nx + x]);
                } else if (datatype == DT_SIGNED_SHORT) {
                    value = static_cast<float>(image->img.ss[z * sizeOfSlice + y * image->region.nx + x]);
                } else {
                    value = image->img.fl[z * sizeOfSlice + y * image->region.nx + x];
                }
                if (value > threshold) {
                    maxx = x;
                    // break out of ALL loops
                    z = maxz + 1;
                    x = minx - 1;
                    break;
                }
            }
        }
    }

    region.ox = image->region.ox + minx;
    region.oy = image->region.oy + miny;
    region.oz = image->region.oz + minz;
    region.nx = maxx - minx + 1;
    region.ny = maxy - miny + 1;
    region.nz = maxz - minz + 1;

    return region;
}

// ---------------------------------------------------------------------------
ImageRegion join_regions(ImageRegion region1, ImageRegion region2)
{
    ImageRegion region = region1;

    int maxx  = region .ox + region .nx;
    int maxy  = region .oy + region .ny;
    int maxz  = region .oz + region .nz;
    int maxx2 = region2.ox + region2.nx;
    int maxy2 = region2.oy + region2.ny;
    int maxz2 = region2.oz + region2.nz;

    if (region2.ox < region.ox) region.ox = region2.ox;
    if (region2.oy < region.oy) region.oy = region2.oy;
    if (region2.oz < region.oz) region.oz = region2.oz;

    if (maxx2 > maxx) region.nx = maxx2 - region.ox + 1;
    if (maxy2 > maxy) region.ny = maxy2 - region.oy + 1;
    if (maxz2 > maxz) region.nz = maxz2 - region.oz + 1;

    return region;
}

// ===========================================================================
// cropping / padding of image
// ===========================================================================

// ---------------------------------------------------------------------------
int resize(Image* image, ImageRegion region, float padValue)
{
    VoxelData newImg; // data of resized image 
    VoxelData out;    // output image "iterator"

    newImg.fl = NULL;
    out   .fl = NULL;

    short datatype = image->hdr.datatype;

    // image must be owner of the data
    if (!image->owns_img) return -1;

    // make region valid if necessary
    if (region.ox < 0) region.ox = 0;
    if (region.oy < 0) region.oy = 0;
    if (region.oz < 0) region.oz = 0;

    if (region.ox + region.nx > image->hdr.dim [1]) region.nx = image->hdr.dim [1] - region.ox;
    if (region.oy + region.ny > image->hdr.dim [2]) region.ny = image->hdr.dim [2] - region.oy;
    if (region.oz + region.nz > image->hdr.dim [3]) region.nz = image->hdr.dim [3] - region.oz;

    // allocate memory for resized image data
    if (datatype == DT_UNSIGNED_CHAR) {
        newImg.uc = reinterpret_cast<unsigned char*>(malloc(region.nx * region.ny * region.nz));
    } else if (datatype == DT_SIGNED_SHORT) {
        newImg.ss = reinterpret_cast<short*>(malloc(region.nx * region.ny * region.nz));
    } else if (datatype == DT_FLOAT) {
        newImg.fl = reinterpret_cast<float*>(malloc(region.nx * region.ny * region.nz));
    } else {
        return -2;
    }

    if (newImg.fl == NULL) return -3;
 
    // copy image data within specified region
    out = newImg;
    for (int z = region.oz; z < region.nz; ++ z) {
        if (image->region.oz <= z && z < (image->region.oz + image->region.nz)) {
            for (int y = region.oy; y < region.ny; ++ y) {
                if (image->region.oy <= y && y < (image->region.oy + image->region.ny)) {
                    for (int x = region.ox; x < region.nx; ++ x) {
                        if (image->region.ox <= x && x < (image->region.ox + image->region.nx)) {
                            if (datatype == DT_UNSIGNED_CHAR) {
                                *out.uc = image->img.uc [  (z - image->region.oz) * image->hdr.dim [1] * image->hdr.dim [2]
                                                         + (y - image->region.oy) * image->hdr.dim [1]
                                                         + (x - image->region.ox)];
                                ++ out.uc;
                            } else if (datatype == DT_SIGNED_SHORT) {
                                *out.ss = image->img.ss [  (z - image->region.oz) * image->hdr.dim [1] * image->hdr.dim [2]
                                                         + (y - image->region.oy) * image->hdr.dim [1]
                                                         + (x - image->region.ox)];
                                ++ out.ss;
                            } else {
                                *out.fl = image->img.fl [  (z - image->region.oz) * image->hdr.dim [1] * image->hdr.dim [2]
                                                         + (y - image->region.oy) * image->hdr.dim [1]
                                                         + (x - image->region.ox)];
                                ++ out.fl;
                            }
                        } else {
                            if      (datatype == DT_UNSIGNED_CHAR) *out.uc = static_cast <unsigned char> (padValue);
                            else if (datatype == DT_SIGNED_SHORT)  *out.ss = static_cast <short>         (padValue);
                            else                                   *out.fl = padValue;
                        }
                    }
                } else {
                    if (datatype == DT_UNSIGNED_CHAR) {
                        for (int x = region.ox; x < region.nx; ++ x) {
                            *out.uc = static_cast <unsigned char> (padValue);
                            ++ out.uc;
                        }
                    } else if (datatype == DT_SIGNED_SHORT) {
                        for (int x = region.ox; x < region.nx; ++ x) {
                            *out.ss = static_cast <short> (padValue);
                            ++ out.ss;
                        }
                    } else {
                        for (int x = region.ox; x < region.nx; ++ x) {
                            *out.fl = padValue;
                            ++ out.fl;
                        }
                    }
                }
            }
        } else {
            if (datatype == DT_UNSIGNED_CHAR) {
                for (int y = region.oy; y < region.ny; ++ y) {
                    for (int x = region.ox; x < region.nx; ++ x) {
                        *out.uc = static_cast <unsigned char> (padValue);
                        ++ out.uc;
                    }
                }
            } else if (datatype == DT_SIGNED_SHORT) {
                for (int y = region.oy; y < region.ny; ++ y) {
                    for (int x = region.ox; x < region.nx; ++ x) {
                        *out.ss = static_cast <short> (padValue);
                        ++ out.ss;
                    }
                }
            } else {
                for (int y = region.oy; y < region.ny; ++ y) {
                    for (int x = region.ox; x < region.nx; ++ x) {
                        *out.fl = padValue;
                        ++ out.fl;
                    }
                }
            }
        }
    }
 
    // replace image data
    free(image->img.fl);
    image->img    = newImg;
    image->region = region;
 
    return 0;
}

// ---------------------------------------------------------------------------
int crop(Image* image, ImageRegion region)
{
    ImageRegion region2 = region;

    int maxx = image->region.ox + image->region.nx;
    int maxy = image->region.oy + image->region.ny;
    int maxz = image->region.oz + image->region.nz;

    region2.ox = image->region.ox + region.ox;
    region2.oy = image->region.oy + region.oy;
    region2.oz = image->region.oz + region.oz;

    if (region2.ox + region2.nx > maxx) region2.nx = maxx - region2.ox + 1;
    if (region2.oy + region2.ny > maxy) region2.ny = maxy - region2.oy + 1;
    if (region2.oz + region2.nz > maxz) region2.nz = maxz - region2.oz + 1;

    return resize(image, region2, 0.0);
}

// ---------------------------------------------------------------------------
int pad(Image* image, float value)
{
    ImageRegion region;

    region.ox = 0;
    region.oy = 0;
    region.oz = 0;
    region.nx = image->hdr.dim[1];
    region.ny = image->hdr.dim[2];
    region.nz = image->hdr.dim[3];

    return resize(image, region, value);
}

// ===========================================================================
// matrix operations
// ===========================================================================

// ---------------------------------------------------------------------------
void invert_matrix(float* matrix, float* inverse, int n) 
{
    // Notes:
    // - the matrix must be invertible
    // - there's no pivoting of rows or columns, hence,
    //   accuracy might not be adequate for your needs.
    //
    // Code is rewritten from c++ template code Mike Dinolfo

    int    i, j, k; // loop variables
    double sum,x;   // sum variables
 
    // copy the input matrix to output matrix
    for(i = 0; i < n*n; i++) {
	    inverse [i] = matrix [i];
    }

    // add small value to diagonal if diagonal is zero
    for(i = 0; i < n; i++) { 
        j = i * n + i;
        if ((inverse [j] < 1e-12) && (inverse [j] > -1e-12)) {
             inverse [j] = 1e-12;
        }
    }
 
    // matrix size must be larger than one
    if (n <= 1) return;
 
    // normalize row 0
    for (i = 1; i < n; i++) {
        inverse [i] /= inverse [0]; 
    }

    for (i = 1; i < n; i++) {
        // do a column of L
        for (j = i; j < n; j++) {
            sum = 0.0;
            for (int k = 0; k < i; k++) {
                sum += inverse [j * n + k] * inverse [k * n + i];
            }
            inverse [j * n + i] -= sum;
        }
        if (i == n - 1) continue;
        // do a row of U
        for (int j = i + 1; j < n; j++) {
            sum = 0.0;
            for (int k = 0; k < i; k++) {
                sum += inverse [i * n + k] * inverse [k * n + j];
            }
            inverse [i * n + j] = (inverse [i * n + j] - sum) / inverse [i * n + i];
        }
    }
    // invert L
    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            x = 1.0;
            if ( i != j ) {
                x = 0.0;
                for (k = i; k < j; k++) {
                    x -= inverse [j * n + k] * inverse [k * n + i];
                }
            }
            inverse [j * n + i] = x / inverse [j * n + j];
        }
    }
    // invert U
    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            if (i == j) continue;
            sum = 0.0;
            for (k = i; k < j; k++) {
                sum += inverse [k * n + j] * ((i == k) ? 1.0 : inverse [i * n + k]);
            }
            inverse [i * n + j] = -sum;
        }
    }
    // final inversion
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            sum = 0.0;
            for (k = ((i > j) ? i : j); k < n; k++) {
                sum += ((j == k) ? 1.0 : inverse [j * n + k]) * inverse [k * n + i];
            }
            inverse [j * n + i] = sum;
        }
    }
}

// ===========================================================================
// convolution
// ===========================================================================

// ---------------------------------------------------------------------------
float* convolve(float* K, float* r, int lk, int lr)
{
    int l,i,j;

    l=lk+lr-1;
    float *out=NULL;
    out=(float *)malloc(lk*sizeof(float));

    float *out1=NULL;
    out1=(float *)malloc(l*sizeof(float));
    memset(out1,0,l*sizeof(float));

    float *r1=NULL;
    r1=(float *)malloc(l*sizeof(float));
    memset(r1,0,l*sizeof(float));

    float *r2=NULL;
    r2=(float *)malloc(l*sizeof(float));
    memset(r2,0,l*sizeof(float));

    for(i=0;i<lk;i++) r1[i]=K[i];
    for(i=0;i<lr;i++) r2[i]=r[i];

    for(i=0;i<l;i++) {
       for(j=0;j<=i;j++) {
          out1[i]+=r1[j]*r2[i-j];
       }
    }

    for(i=0;i<lk;i++) out[i]=out1[i+2];

    free(out1);
    free(r1);
    free(r2);
   
    return out;
}


} // namespace mico
