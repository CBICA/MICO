/**
 * @file  longitudinal.cxx
 * @brief Implements MICO for longitudinal MR brain image series.
 *
 * Copyright (c) 2011, 2012 University of Pennsylvania. All rights reserved.<br />
 * See https://www.med.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <software at cbica.upenn.edu>
 */

#include <time.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>

//#include <basis/assert.h>

#include "mico/utilities.h"
#include "mico/cross-sectional.h"
#include "mico/longitudinal.h"


// acceptable in .cxx file
using namespace std;


namespace mico {


// ===========================================================================
// API
// ===========================================================================

// ---------------------------------------------------------------------------
bool segment_longitudinal(const ImageVector& images,
                          ImageMapVector&    memberships,
                          const Parameters&  params,
                          int                verbosity)
{
    bool ok = true;

    const int number_of_images = static_cast<int>(images.size());

    // -----------------------------------------------------------------------
    // degraded image series
    if (number_of_images < 1) {
        return false;
    } else if (number_of_images < 2) {
        memberships.resize(1);
        return segment_and_bias_correct(images[0], &memberships[0], NULL, params, verbosity);
    }

    // -----------------------------------------------------------------------
    // reset output
    delete_images(memberships);

    // -----------------------------------------------------------------------
    // common image size
    const int width  = images[0]->hdr.dim[1];
    const int height = images[0]->hdr.dim[2];
    const int slices = images[0]->hdr.dim[3];

    const int total = width * height * slices;

    // -----------------------------------------------------------------------
    // determine union of foreground regions
    ImageRegion foreground;

    if (ok) {
        if (verbosity > 0) {
            printf ("\nDetermining common region-of-interest...");
            fflush (stdout);
        }

        for (int i = 0; i < number_of_images; i++) {
            ImageRegion region = get_foreground_region(images[i], params.th_bg);
            if (i > 0) {
                foreground = join_regions(foreground, region);
            } else {
                foreground = region;
            }
        }

        if (ok) {
            if (verbosity > 0) printf(" - done");
            if (verbosity > 1) {
                printf(": offset = [%d, %d, %d], size = [%d, %d, %d]",
                       foreground.ox, foreground.oy, foreground.oz,
                       foreground.nx, foreground.ny, foreground.nz);
            }
            if (verbosity > 0) {
                printf("\n");
                fflush(stdout);
            }

            if (foreground.nx == 0 || foreground.ny == 0 || foreground.nz == 0) {
                fprintf(stderr, "WARNING: Region-of-interest has zero dimension!\n");
                fprintf(stderr, "WARNING: Check your images that they contain values higher than %f.\n", params.th_bg);
                fprintf(stderr, "WARNING: Possibly rescale intensities linearly if applicable.\n");
            }

        }
    }

    // -----------------------------------------------------------------------
    // calculate b, C, and M
    if (ok && verbosity > 0) printf("\n");

    // Note: Calculate the b, C, and M functions for the current image
    //       and store these in the results variable. Once we have
    //       calculated these functions for all images, we will update
    //       them by taking the temporal consistency constraint into
    //       consideration. This is done in the final step.
    float** results = NULL; // b, C, and M

    if (ok) {
        results = (float**)calloc(number_of_images, sizeof(float*));
        if (results == NULL) {
            fprintf(stderr, "Failed to allocate memory!\n");
            ok = false;
        }
    }

    float* bias = NULL;
    if (ok) {
        bias = get_basis_order(foreground.nx, foreground.ny, foreground.nz);
        if (bias == NULL) {
            fprintf(stderr, "Failed to initialize basis functions!\n");
            ok = false;
        }
    }

    if (ok) {
        for (int i = 0; ok && i < number_of_images; i++) {
            Image* const image = images[i];

            if (verbosity > 0) {
                printf("Segmenting brain image %d of %d (%s)...", i + 1, number_of_images, image->name.c_str());
                fflush(stdout);
            }

            // cast image to DT_FLOAT
            Image* flimage = image;
            if (image->hdr.datatype != DT_FLOAT) {
                flimage = cast_image(image, DT_FLOAT);
                if (flimage == NULL) {
                    fprintf(stderr, "Not enough memory available to convert image to datatype DT_FLOAT!\n");
                    ok = false;
                    break;
                }
            }
            float* const img = flimage->img.fl;

            // calculate b, C, and M functions
            results[i] = (float *)calloc (4 * total + 3, sizeof(float));
            if (results[i] == NULL) {
                if (verbosity > 0) {
                    printf(" - failed\n");
                    fflush(stdout);
                }
                fprintf(stderr, "Failed to allocate memory!\n");
                ok = false;
                break;
            }
            updateBCM(foreground, bias, img, width, height, slices, params, results[i]);

            // destroy temporary image
            if (flimage != image) delete_image(flimage);

            // shield background for b
            for (int j = 0; j < total; j++) {
                if (img[j] <= params.th_bg) results[i][j] = 0.0;
            }
            // shield background for all M
            float* p = NULL;
            for (int j = total + 3; j < 4 * total + 3; j++, p++) {
                if (((j - 3) % total) == 0) p = img;
                if (*p <= params.th_bg) results[i][j] = 0.0;
            }

            // rescale
            float r_rescale = 200.0 / results[i][total + 2];

            for (int j = 0;     j  < total;     j++) img       [j] *= r_rescale;
            for (int j = total; j <= total + 2; j++) results[i][j] *= r_rescale;

            if (verbosity > 0) {
                printf(" - done\n");
                fflush(stdout);
            }
        }

        // clean up
        free(bias);
        bias = NULL;
    }

    // -----------------------------------------------------------------------
    // calculate temporal consistent membership functions
    if (ok) {
        if (verbosity > 0) {
            printf("\nEstablishing temporal consistency...");
            fflush(stdout);
        }

        memberships = tc_update_membership(images, results, params);

        if (memberships.empty()) {
            fprintf(stderr, "Not enough memory available!\n");
            ok = false;
        }

        if (verbosity > 0) {
            printf(ok ? " - done\n" : " - failed\n");
            fflush(stdout);
        }
    }

    // -----------------------------------------------------------------------
    // clean up
    if (results) {
        for (int i = 0; i < number_of_images; i++) free(results[i]);
        free(results);
        results = NULL;
    }

    return ok;
}

// ===========================================================================
// temporally consistency
// ===========================================================================

// ---------------------------------------------------------------------------
ImageMapVector tc_update_membership(const ImageVector& images,
                                    float**            result_ser,
                                    const Parameters&  params)
{
    bool ok = true;

    float kernel[5] = {1, 2, 4, 2, 1};

    const int number_of_images = static_cast<int>(images.size());

  if (number_of_images <= 0)
  {
    std::cerr << "Number of images is less than or equal to 0, something went wrong.\n";
    exit(EXIT_FAILURE);
  }
  
    // assert(result_ser != NULL);

    const int width  = images[0]->hdr.dim[1];
    const int height = images[0]->hdr.dim[2];
    const int slices = images[0]->hdr.dim[3];

    const int total = width * height * slices;

    // normalize kernel
    float sum_kernel = 0;
    for (int i = 0; i < 5; i++) sum_kernel += kernel[i];
    for (int i = 0; i < 5; i++) kernel[i]   = kernel[i] / sum_kernel;
    // ...
    float** Rs = (float**)calloc(number_of_images, sizeof(float*));
    if (Rs == NULL) ok = false;
    for (int i = 0; ok && i < number_of_images; i++) {
        Rs[i] = (float*)calloc (3 * total, sizeof(float));
        if (Rs[i] == NULL) ok = false;
    }
    if (ok) {
        float* img;
        float* result;
        for (int i = 0; i < number_of_images; i++) 
        {
          
          if (images[i]->hdr.datatype != DT_FLOAT)
          {
            std::cerr << "Data type is not float.\n";
            exit(EXIT_FAILURE);
          }
            // assert(result_ser[i] != NULL);
            float* Rs_i = Rs[i];
            // CSF
            img    = images[i]->img.fl;
            result = result_ser[i];
            for (int j = 0; j < total; j++) {
                *(Rs_i++) = params.tissueWeight[0] * pow(img[j] - result[total] * result[j], 2.0);
            }
            // GM
            img    = images[i]->img.fl;
            result = result_ser[i];
            for (int j = 0; j < total; j++) {
                *(Rs_i++) = params.tissueWeight[1] * pow(img[j] - result[total + 1] * result[j], 2.0);
            }
            // WM
            img    = images[i]->img.fl;
            result = result_ser[i];
            for (int j = 0; j < total; j++) {
                *(Rs_i++) = params.tissueWeight[2] * pow(img[j] - result[total + 2] * result[j], 2.0);
            }
        }
    }
    float** Rs_out = (float**)calloc(number_of_images, sizeof(float*));
    if (Rs_out == NULL) ok = false;
    for (int i = 0; ok && i < number_of_images; i++) {
        Rs_out[i] = (float*)calloc (3 * total, sizeof(float));
        if (Rs_out[i] == NULL) ok = false;
    }
    // convolute the membership function with the kernel
    for (int k = 0; ok && k < 3 * total; k++) {
        float* KK = (float*)malloc(number_of_images * sizeof(float));
        if (KK == NULL) {
            ok = false;
            break;
        }
        for (int i = 0; i < number_of_images; i++) KK[i] = Rs[i][k];
        float* out = convolve(KK, kernel, number_of_images, 5);
        if (out == NULL) {
            ok = false;
        } else {
            for (int i = 0; i < number_of_images; i++) Rs_out[i][k] = out[i];
        }
        free(KK);
    }
    // initialize membership functions
    ImageMapVector memberships(number_of_images);
    for (int i = 0; ok && i < number_of_images; i++) {
        memberships[i]["CSF"] = new_image(width, height, slices, DT_FLOAT);
        memberships[i]["GM" ] = new_image(width, height, slices, DT_FLOAT);
        memberships[i]["WM" ] = new_image(width, height, slices, DT_FLOAT);
        for (ImageMap::iterator it = memberships[i].begin (); it != memberships[i].end (); ++it) {
            Image* membership = it->second;
            if (membership == NULL) {
                fprintf(stderr, "Failed to allocate memory for membership function of %s!\n", it->first);
                ok = false;
                break;
            }
            membership->hdr          = images[0]->hdr;
            membership->hdr.bitpix   = 8 * sizeof(float);
            membership->hdr.datatype = DT_FLOAT;
            membership->compress     = images[0]->compress;

            string descrip = "fuzzy ";
            descrip += it->first;
            descrip += " membership (MICO)";
            strncpy(membership->hdr.descrip, descrip.c_str(), descrip.size() + 1);
        }
    }
    // calculate temporally consistent membership functions
    for (int i = 0; i < number_of_images; i++) {
        const float epsilon = 0.01;
        for (int j = 0; j < 3 * total; j++) Rs_out[i][j] += epsilon;
        float p = 1.0 / (params.q - 1.0);
        float* f = (float*)calloc(3 * total, sizeof(float));
        if (f == NULL) {
            ok = false;
            break;
        }
        for (int j = 0; j < 3 * total; j++) f[j] = 1.0 / pow(Rs_out[i][j], p);
        float* f_sum = (float*)calloc(total, sizeof(float));
        if (f_sum == NULL) {
            ok = false;
            free(f);
            break;
        }
        for (int j = 0; j < total; j++) {
            for (int k = 0; k < 3; k++) f_sum[j] += f[j + k * total];
        }
        for (int k = 0; k < 3; k++) {
            float* M = NULL;
            if      (k == 0) M = memberships[i]["CSF"]->img.fl;
            else if (k == 1) M = memberships[i]["GM"] ->img.fl;
            else if (k == 2) M = memberships[i]["WM"] ->img.fl;
            // assert (M != NULL);
            for (int j = 0; j < total; j++) {
                if (result_ser[i][j] != 0) {
                    M[j] = f[k * total + j] / f_sum[j];
                } else {
                    M[j] = 0;
                }
            }
        }
	    free(f);
	    free(f_sum);
    }
    // clean up
    for (int i = 0; i < number_of_images; i++) {
        free(Rs[i]);
        free(Rs_out[i]);
    }
    free(Rs);
    free(Rs_out);
    // done
    if (!ok) delete_images(memberships);
    return memberships;
}


} // namespace mico
