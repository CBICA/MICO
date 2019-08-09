/**
 * @file  mico.cxx
 * @brief Command-line tool.
 *
 * This program implements Chunming Li et al's 3D segmentation and in particular
 * 4D segmentation. The latter for longitudinal consistent segmentations of
 * longitudinal MR brain series images.
 *
 * Copyright (c) 2011, 2012 University of Pennsylvania. All rights reserved.<br />
 * See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */


#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <algorithm>
#include <string>

#include <mico/mico.h>

#include "cbicaUtilities.h"

// #include <mico/basis.h> // print_contact(), print_version()


// acceptable in .cxx file
using namespace std;
// using namespace basis;
using namespace mico;


// ===========================================================================
// constants
// ===========================================================================

const char* const cSegmentationImageSuffix  = "_labels";
const char* const cBiasCorrectedImageSuffix = "_bc";
const char* const cMembershipImageSuffix    = "_<label>";

// NIfTI image header description strings need to be less than 80 characters long
const char* const cSegmentationImageDescription  = "Segmentation obtained by MICO <RELEASE>";
const char* const cBiasCorrectedImageDescription = "Bias corrected using MICO <RELEASE>";
const char* const cMembershipImageDescription    = "Fuzzy segmentation of <LABEL> obtained by MICO <RELEASE>";

// ===========================================================================
// version / help / usage
// ===========================================================================

/**
 * @brief Print options.
 *
 * @sa print_help()
 * @sa print_usage()
 */
void print_options()
{
    Parameters dflt; // default parameters
    printf("Required options:\n");
    printf("  --inputlist, -i <file>           Text file listing the input T1 image files.\n");
    printf("                                   Relative image file paths are considered relative to\n");
    printf("                                   the directory of the input list file.\n");
    printf("  <T1_*.nii>                       Instead of a text file listing the input image files,\n");
    printf("                                   the images can be given as arguments on the command-line.\n");
    printf("\n");
    printf("Options:\n");
    printf("  --outputdir, -o <dir>            Output directory for obtained segmentations.\n");
    printf("                                   Defaults to current working directory.\n");
    printf("  --4d, -4                         Whether to segment the images as a longitudinal series.\n");
    printf("                                   By default, an independent 3D segmentation without\n");
    printf("                                   temporal consistency constraint is computed.\n");
    printf("  --iterations, -n <int>           The number of iterations (default: %d).\n", dflt.iterNum);
    printf("  --csf-weight, -c <float>         The weight for CSF (default: %.2f).\n", dflt.tissueWeight[0]);
    printf("  --gm-weight, -g <float>          The weight for GM  (default: %.2f).\n", dflt.tissueWeight[1]);
    printf("  --wm-weight, -w <float>          The weight for WM  (default: %.2f).\n", dflt.tissueWeight[2]);
    printf("  --csf-label <0..255>             Label used for CSF (default: %d).\n",   dflt.labelCSF);
    printf("  --gm-label <0..255>              Label used for GM  (default: %d).\n",   dflt.labelGM);
    printf("  --wm-label <0..255>              Label used for WM  (default: %d).\n",   dflt.labelWM);
    printf("  --lambda, -l <float>             The parameter lambda > 0 (default: 0.0001).\n");
    printf("  --threshold, -t <float>          Foreground threshold. Only intensities exceeding this threshold\n");
    printf("                                   are considered foreground. (default: %.1f)\n", dflt.th_bg);
    printf("  --suffix, -s <suffix>            Specify suffix to use for output segmentation files.\n");
    printf("                                   (default: %s)\n", cSegmentationImageSuffix);
    printf("  --fuzzy, -f                      Whether to output fuzzy segmentations, one for each segment.\n");
    printf("  --fuzzy-suffix <suffix>          File name suffix to use for the fuzzy segmentations, where the\n");
    printf("                                   string '<label>' or '<LABEL>' in the suffix is replaced by the\n");
    printf("                                   name of the respective structure in lower- or uppercase, respectively.\n");
    printf("                                   If neither '<label>' nor '<LABEL>' is found in the suffix string,\n");
    printf("                                   '_<label>' is appended to the base file name.\n");
    printf("                                   (default: %s)\n", cMembershipImageSuffix);
    printf("  --bias-correct, -b               Whether to output bias corrected images.\n");
    printf("  --bias-correct-only              Use instead of the option --bias-correct in order to disable\n");
    printf("                                   the output of the segmentations, but only output the bias\n");
    printf("                                   corrected images.\n");
    printf("  --bias-correct-suffix <suffix>   File name suffix to use for bias corrected images.\n");
    printf("                                   (default: %s)\n", cBiasCorrectedImageSuffix);
    printf("  --n1                             Force output of .nii image files in NIfTI-1 format.\n");
    printf("  --n2                             Force output of .hdr/.img image files in NIfTI-1 format.\n");
    printf("  --zn1                            Force output of .nii.gz image files in NIfTI-1 format.\n");
    printf("  --zn2                            Force output of .hdr.gz/.img.gz image files in NIfTI-1 format.\n");
    printf("  --srand-time                     Use the execution time to initialize the pseudo-random\n");
    printf("                                   number generator. If this option is not given, it is\n");
    printf("                                   initialized using a fixed seed which ensures identical\n");
    printf("                                   results whenever this program is run on the same machine\n");
    printf("                                   with the same input images and parameter settings.\n");
    printf("                                   Otherwise, the results may differ slightly between\n");
    printf("                                   executions. Note that the differences decrease, however,\n");
    printf("                                   with the number of --iterations.\n");
    printf("  --verbose, -v                    Increase verbosity of output messages.\n");
    printf("                                   Can be given multiple times.\n");
    printf("  --help, -h                       Print help and exit.\n");
    printf("  --helpshort                      Print short help and exit.\n");
    printf("  --version                        Print version and exit.\n");
}

/**
 * @brief Print help.
 */
void print_help()
{
    printf("Usage:\n");
    printf("  mico [options] [<T1_1.nii> [<T1_2.nii> ...]]\n");
    printf("\n");
    printf("Description:\n");
    printf("  This command implements the Multiplicative Intrinsic Component Optimization (MICO)\n");
    printf("  algorithm. In particular, it segments a set of MR brain images for two cases:\n");
    printf("\n");
    printf("    1) The images are acquired from different subjects:\n");
    printf("       In this case, the program segments each input volume independently and\n");
    printf("       hence performs a normal 3D segmentation.\n");
    printf("\n");
    printf("    2) The images are acquired from the same subject at different time points\n");
    printf("       (i.e., a longitudinal series, or 4D data) and the -s option is given:\n");
    printf("       In this case, the program takes into account a temporal consistency constraint\n");
    printf("       in the segmentation to achieve a temporally consistent segmentation of\n");
    printf("       the longitudinal series. Hence, a 4D segmentation is performed.\n");
    printf("\n");
    printf("  If a tissue is over or under segmented, the weight of this tissue needs to be\n");
    printf("  increased or decreased, respectively.\n");
    printf("\n");
    printf("  For images with almost no inhomogeneity (like 1.5T MR images), a large value for\n");
    printf("  lambda, such as 10, can be used. For images with significant inhomogeneity\n");
    printf("  (like 3T MR images), a very small value for lambda, such as 0.0001, should be used\n");
    printf("  or simply lambda be set to 0.\n");
    printf("\n");
    printf("Image formats:\n");
    printf("  This program reads both ANALYZE 7.5 and NIfTI-1 with a scalar datatype of either\n");
    printf("  DT_UNISGNED_CHAR, DT_SIGNED_SHORT, DT_SIGNED_INT, or DT_FLOAT. If no specific file\n");
    printf("  name extension is specified for the output files using one of the --*suffix options,\n");
    printf("  this program will write any outputs related to a particular input in the same format\n");
    printf("  as used for the storage of the input image. Otherwise, the specified file name extension\n");
    printf("  is considered when determining the output file format. Additionally, one of the options\n");
    printf("  --n1, --n2, --zn1, and --zn2 can be used to choose either one of the NIfTI-1 file formats\n");
    printf("  as default for any output for which no explicitly specified file name extension is given.\n");
    printf("\n");
    print_options();
    printf("\n");
    printf("Examples:\n");
    printf("  mico brain.nii\n");
    printf("    Segments the MR image brain.nii and writes the resulting label image\n");
    printf("    brain%s.nii to the current working directory.\n", cSegmentationImageSuffix);
    printf("\n");
    printf("  mico --suffix _3Dseg.nii.gz brain.hdr\n");
    printf("    Segments the MR image brain.hdr/brain.img and writes the resulting label image\n");
    printf("    to the compressed file brain_3Dseg.nii.gz in the current working directory.\n");
    printf("\n");
    printf("  mico --fuzzy --fuzzy-suffix _<label>.nii.gz brain.hdr\n");
    printf("    Segments the MR image brain.hdr/brain.img and writes the resulting fuzzy\n");
    printf("    segmentations to the compressed files brain_csf.nii.gz, brain_gm.nii.gz,\n");
    printf("    and brain_wm.nii.gz in the current working directory.\n");
    printf("\n");
    printf("  mico --4d brain_1.nii brain_2.nii brain_3.nii brain_4.nii brain_5.nii\n");
    printf("    Segments the given longitudinal series of brain images by taking a temporal\n");
    printf("    consistency constraint into account. The resulting segmentations are named\n");
    printf("    after the input images with the suffix \"%s\" and extension \".nii\".\n", cSegmentationImageSuffix);
    printf("\n");
    printf("  mico --4d --inputlist images.lst --outputdir segmentations --lambda 10\n");
    printf("    Segments the longitudinal image series specified by the input list file\n");
    printf("    images.lst (a plain text file) with almost no inhomogeneity (therefore, the\n");
    printf("    choice of lambda = 10) and writes the resulting segmentations into the\n");
    printf("    subdirectory \"segmentations\" under the current working directory.\n");
    printf("\n");
    printf("  mico --bias-correct brain.hdr\n");
    printf("    Segments and bias corrects the MR brain image brain.hdr/brain.img and writes the\n");
    printf("    resulting label image brain_labels.hdr/brain_segments.img and the bias corrrected\n");
    printf("    image brain_biascorrected.hdr/brain_bc.img.\n");
    printf("\n");
    printf("  mico --bias-correct-only --zn1 brain.hdr\n");
    printf("    Bias corrects the MR brain image brain.hdr/brain.img and writes the resulting\n");
    printf("    bias corrected image to the compressed file brain_bc.nii.gz.\n");
    printf("\n");
    std::cout << "Contact software@cbica.upenn.edu";
}

/**
 * @brief Print usage information.
 */
void print_usage()
{
    printf("Usage:\n");
    printf("  mico [options] [<T1_1.nii> [<T1_2.nii> ...]]\n");
    printf("\n");
    print_options();
    printf("\n");
    std::cout << "Contact software@cbica.upenn.edu";
}

// ===========================================================================
// auxiliary functions
// ===========================================================================

// ---------------------------------------------------------------------------
string get_suffix_with_extension(const char* suffix, const char* extension)
{
    string suffix_out(suffix);
    // if an explicit extension is given (--hdr and --nii options) and the
    // suffix itself does not have a valid extension...
    std::set<string> niftiexts;
    niftiexts.insert(".hdr");
    niftiexts.insert(".hdr.gz");
    niftiexts.insert(".img");
    niftiexts.insert(".img.gz");
    niftiexts.insert(".nii");
    niftiexts.insert(".nii.gz");
    niftiexts.insert(".nia");
    niftiexts.insert(".nia.gz");
    // append extension to suffix
    if (extension /*&& !os::path::hasext(suffix, &niftiexts)*/) suffix_out += extension;
    // use .hdr instead of .img
    if (suffix_out.length() >= 4 && suffix_out.substr(suffix_out.length() - 4) == ".img") {
        suffix_out.replace(suffix_out.length() - 4, 4, ".hdr");
    }
    // use .hdr.gz instead of .img.gz
    if (suffix_out.length() >= 7 && suffix_out.substr(suffix_out.length() - 7) == ".img.gz") {
        suffix_out.replace(suffix_out.length() - 7, 7, ".hdr.gz");
    }
    return suffix_out;
}

// ---------------------------------------------------------------------------
string replace_label(const char* str, const char* label)
{
    string text    (str   ? str   : "");
    string label_lc(label ? label : "");
    string label_uc(label ? label : "");
    transform(label_lc.begin(), label_lc.end(), label_lc.begin(), ::tolower);
    transform(label_uc.begin(), label_uc.end(), label_uc.begin(), ::toupper);
    while (true) {
        size_t pos = text.find("<label>");
        if (pos != string::npos) {
            text.replace(pos, 7, label_lc);
            continue;
        }
        pos = text.find("<LABEL>");
        if (pos != string::npos) {
            text.replace(pos, 7, label_uc);
            continue;
        }
        break;
    }
    return text;
}

// ---------------------------------------------------------------------------
string replace_release(const char* str)
{
    string text(str ? str : "");
    while (true) {
        size_t pos = text.find("<RELEASE>");
        if (pos != string::npos) {
            text.replace(pos, 9, "RELEASE");
            continue;
        }
        break;
    }
    return text;
}

// ===========================================================================
// main of cross-sectional segmentation program
// ===========================================================================

/**
 * @brief Main function of 3D segmentation program.
 *
 * Calculate b, C, and M functions for the each image, bias correct it and
 * create the segmentation, i.e., label image from the membership functions M.
 * Output these images and then proceed with the next image.
 */
int main3D(char** filenames, int number_of_images, const Parameters& params, int verbosity)
{
    bool ok = true;

    for (int i = 0; i < number_of_images; i++) {

        // -------------------------------------------------------------------
        // read image
        if (verbosity > 0) {
            printf("\nProcessing image %d of %d\n\n", i + 1, number_of_images);
            printf("Reading image %s...", filenames[i]);
            fflush(stdout);
        }

        nifti_image* nii = nifti_image_read(filenames[i], 1);
        if (nii == NULL) {
            if (verbosity > 0) printf(" - failed\n");
            fprintf(stderr, "Failed to read image %s!\n", filenames[i]);
            ok = false;
            continue;
        }

        if (verbosity > 0) {
            printf(" - done\n");
            fflush(stdout);
        }

        int datatype = nii->datatype; // remember original datatype of input image
                                      // in particular, the bias corrected image
                                      // is converted to this datatype before it
                                      // is written to disk

        Image* image = convert_nifti_image(nii, DT_FLOAT, true);

        if (image == NULL) {
            fprintf(stderr, "Image %s cannot be processed by this program!\n", filenames[i]);
            nifti_image_free(nii); nii = NULL;
            ok = false;
            continue;
        }

        nii = NULL; // NIfTI image data structure no longer valid

        // -------------------------------------------------------------------
        // perform bias correction and/or segmentation
        ImageMap memberships;
        Image* bcimage = NULL;
        if (params.outputBiasCorrections && (params.outputSegments || params.outputMemberships)) {
            ok = segment_and_bias_correct(image, &memberships, &bcimage, params, verbosity);
        } else if (params.outputBiasCorrections) {
            ok = segment_and_bias_correct(image, NULL, &bcimage, params, verbosity);
        } else {
            ok = segment_and_bias_correct(image, &memberships, NULL, params, verbosity);
        }
        // convert bias corrected image
        if (bcimage && bcimage->hdr.datatype != datatype) {
            Image* tmp = cast_image(bcimage, datatype, true);
            if (tmp == NULL) {
                fprintf(stderr, "Not enough memory available to convert bias corrected image to datatype of input image!");
                ok = false;
            }
            delete_image(bcimage);
            bcimage = tmp;
        }
        // create label map
        Image* segmentation = NULL;
        if (!memberships.empty()) {
            segmentation = create_label_map(image, memberships, params, verbosity);
            if (segmentation == NULL) {
                fprintf(stderr, "Not enough memory available to create label map for image %s!\n", image->name.c_str());
                ok = false;
            }
        }
        // skip image on error
        if (!ok) {
            delete_images(memberships);
            if (bcimage) delete_image(bcimage);
            continue;
        }
        // destroy membership functions if output of fuzzy segmentation not requested
        if (!params.outputMemberships) delete_images(memberships);

        // -------------------------------------------------------------------
        // write results

        // build common path prefix for output images
        string prefix  = params.outputDir.empty() ? cbica::constCharToChar(cbica::getCWD()) : params.outputDir;
        prefix        += '/';
        prefix        += image->name;
        // write segmentation image
        if (segmentation) {
            if (verbosity > 0) {
                printf("Writing label image...");
                fflush(stdout);
            }
            string fname = prefix + params.suffix;
            nii = convert_to_nifti(segmentation, fname.c_str(), 0);
            nii->intent_code = NIFTI_INTENT_LABEL;
            string descrip = replace_release(cSegmentationImageDescription);
            strncpy(nii->descrip, descrip.c_str(), 80);
            nifti_image_write(nii);
            if (verbosity > 0) {
                printf(" - done");
                if (verbosity > 1) printf(": %s", nii->fname);
                printf("\n");
                fflush(stdout);
            }
            nifti_image_free(nii);
            delete_image(segmentation);
        }
        // write membership images
        for (ImageMap::iterator it = memberships.begin(); it != memberships.end(); ++it) {
            if (verbosity > 0) {
                printf("Writing segmentation of %s...", it->first);
                fflush(stdout);
            }
            string suffix = replace_label(params.suffixMF.c_str(), it->first);
            if (suffix.empty()) {
                suffix  = '_';
                suffix += it->first;
            }
            string fname = prefix + suffix;
            nii = convert_to_nifti(it->second, fname.c_str(), 0);
            nii->intent_code = NIFTI_FIRST_STATCODE;
            string descrip = replace_label(cMembershipImageDescription, it->first);
            descrip = replace_release(descrip.c_str());
            strncpy(nii->descrip, descrip.c_str(), 80);
            nifti_image_write(nii);
            if (verbosity > 0) {
                printf(" - done");
                if (verbosity > 1) printf(": %s", nii->fname);
                printf("\n");
                fflush(stdout);
            }
            nifti_image_free(nii);
            delete_image(it->second);
            it->second = NULL;
        }
        memberships.clear();
        // write bias corrected image
        if (bcimage) {
            if (verbosity > 0) {
                printf("Writing bias corrected image...");
                fflush(stdout);
            }
            string fname = prefix + params.suffixBC;
            nii = convert_to_nifti(bcimage, fname.c_str(), 0);
            string descrip = replace_release(cBiasCorrectedImageDescription);
            strncpy(nii->descrip, descrip.c_str(), 80);
            nifti_image_write(nii);
            if (verbosity > 0) {
                printf(" - done");
                if (verbosity > 1) printf(": %s", nii->fname);
                printf("\n");
                fflush(stdout);
            }
            nifti_image_free(nii);
            delete_image(bcimage);
        }

        // -------------------------------------------------------------------
        // clean up
        delete_image(image);
        image = NULL;
    }

    return ok ? EXIT_SUCCESS :EXIT_FAILURE;
}

// ===========================================================================
// main of longitudinal segmentation program
// ===========================================================================

/**
 * @brief Main function of 4D segmentation program.
 */
int main4D(char** filenames, int number_of_images, const Parameters& params, int verbosity)
{
    bool ok = true;

    // -----------------------------------------------------------------------
    // read all images into memory
    ImageVector images(number_of_images, NULL);
    int         width  = 0;    // common width of co-registered images
    int         height = 0;    // common height of co-registered images
    int         slices = 0;    // common slices of co-registered images

    if (ok) {
        if (verbosity > 0) printf("\n");
        for (int i = 0; i < number_of_images; i++) {
            if (verbosity > 0) {
                printf("Reading image %d of %d (%s)...", i + 1, number_of_images, filenames[i]);
                fflush(stdout);
            }

            nifti_image* nii = nifti_image_read(filenames[i], 1);
            if (nii == NULL) {
                if (verbosity > 0) printf(" - failed\n");
                fprintf(stderr, "Failed to read image %s!\n", filenames[i]);
                ok = false;
                break;
            }

            if (verbosity > 0) {
                printf(" - done\n");
                fflush(stdout);
            }

            if (i == 0) {
                width  = nii->nx;
                height = nii->ny;
                slices = nii->nz;
            } else {
                if (nii->nx != width || nii->ny != height || nii->nz != slices) {
                    fprintf(stderr, "There is a mismatch in size of the images!\n");
                    ok = false;
                    break;
                }
            }

            images[i] = convert_nifti_image(nii, DT_FLOAT, true);

            if (images[i] == NULL) {
                fprintf(stderr, "Image %s cannot be processed by this program!\n", filenames[i]);
                nifti_image_free(nii);
                ok = false;
                break;
            }
        }
    }

    // -----------------------------------------------------------------------
    // perform temporal consistent segmentation
    ImageMapVector memberships;

    if (ok) {
        ok = segment_longitudinal(images, memberships, params, verbosity);
    }

    // -----------------------------------------------------------------------
    // write results
    if (ok) {
        for (int i = 0; i < number_of_images; i++) {
            // build common path prefix for output images
            string prefix  = params.outputDir.empty() ? cbica::constCharToChar(cbica::getCWD()) : params.outputDir;
            prefix        += '/';
            prefix        += images[i]->name;

            // create label map
            Image* segmentation = NULL;
            if (params.outputSegments) {
                segmentation = create_label_map(images[i], memberships[i], params, verbosity);
                if (segmentation == NULL) {
                    fprintf(stderr, "Not enough memory available to create label map of image %s!", images[i]->name.c_str());
                    ok = false;
                }
            }

            if (!params.outputMemberships) delete_images(memberships[i]);

            // write segmentation image
            if (segmentation) {
                if (verbosity > 0) {
                    printf("Writing label map...");
                    fflush(stdout);
                }
                string fname = prefix + params.suffix;
                nifti_image* nii = convert_to_nifti(segmentation, fname.c_str(), 0);
                nii->intent_code = NIFTI_INTENT_LABEL;
                string descrip = replace_release(cSegmentationImageDescription);
                strncpy(nii->descrip, descrip.c_str(), 80);
                nifti_image_write(nii);
                if (verbosity > 0) {
                    printf(" - done");
                    if (verbosity > 1) printf(": %s", nii->fname);
                    printf("\n");
                    fflush(stdout);
                }
                nifti_image_free(nii);
                delete_image(segmentation);
                segmentation = NULL;
            }

            // write membership images
            for (ImageMap::iterator it = memberships[i].begin(); it != memberships[i].end(); ++it) {
                if (verbosity > 0) {
                    printf("Writing segmentation of %s...", it->first);
                    fflush(stdout);
                }
                string suffix = replace_label(params.suffixMF.c_str(), it->first);
                if (suffix.empty()) {
                    suffix  = '_';
                    suffix += it->first;
                }
                string fname = prefix + suffix;
                nifti_image* nii = convert_to_nifti(it->second, fname.c_str(), 0);
                nii->intent_code = NIFTI_FIRST_STATCODE;
                string descrip = replace_label(cMembershipImageDescription, it->first);
                descrip = replace_release(descrip.c_str());
                strncpy(nii->descrip, descrip.c_str(), 80);
                nifti_image_write(nii);
                if (verbosity > 0) {
                    printf(" - done");
                    if (verbosity > 1) printf(": %s", nii->fname);
                    printf("\n");
                    fflush(stdout);
                }
                nifti_image_free(nii);
                delete_image(it->second);
                it->second = NULL; // avoid 2nd deletion upon clean up
            }
            memberships[i].clear();
        }
    }

    // clean up
    delete_images(images);
    delete_images(memberships);

    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}

// ===========================================================================
// main
// ===========================================================================

// ---------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    nifti_set_debug_level(0); // discard debug message of nifti_is_complete_filename()

    // print help if command run without arguments
    if (argc == 1) {
        print_help();
        exit(EXIT_FAILURE);
    }

    // get current working directory - used to make file paths absolute
    char* cwd    = cbica::constCharToChar(cbica::getCWD());
    int   lencwd = strlen(cwd);

    // -----------------------------------------------------------------------
    // command-line arguments
 
    // default options
    const char* listfile  = NULL;                      // path of image list file
    const char* suffix    = cSegmentationImageSuffix;  // suffix for segmentations
    const char* suffixBC  = cBiasCorrectedImageSuffix; // suffix for bias corrected images
    const char* suffixMF  = cMembershipImageSuffix;    // suffix for membership images
    const char* extension = NULL;                      // common output image file extension
    int         judge     = 3;                         // either 3(D) or 4(D)
    Parameters  params;                                // parameters for main functions
    int         verbosity = 0;                         // verbosity of output messages

    // names of long options
    static struct option long_options[] =
    {
        {"inputlist",           required_argument, NULL, 'i'},
        {"outputdir",           required_argument, NULL, 'o'},
        {"4d",                  no_argument,       NULL, '4'},
        {"iterations",          required_argument, NULL, 'n'},
        {"csf-weight",          required_argument, NULL, 'c'},
        {"gm-weight",           required_argument, NULL, 'g'},
        {"wm-weight",           required_argument, NULL, 'w'},
        {"csf-label",           required_argument, NULL, '/'},
        {"gm-label",            required_argument, NULL, '|'},
        {"wm-label",            required_argument, NULL, '_'},
        {"lambda",              required_argument, NULL, 'l'},
        {"fuzzy",               no_argument,       NULL, 'f'},
        {"bias-correct",        no_argument,       NULL, 'b'},
        {"bias-correct-only",   no_argument,       NULL, 'y'},
        {"suffix",              required_argument, NULL, 's'},
        {"fuzzy-suffix",        required_argument, NULL, 's'},
        {"bias-correct-suffix", required_argument, NULL, 's'},
        {"n2",                  no_argument,       NULL, 'F'},
        {"n1",                  no_argument,       NULL, 'F'},
        {"zn2",                 no_argument,       NULL, 'F'},
        {"zn1",                 no_argument,       NULL, 'F'},
        {"srand-time",          no_argument,       NULL, 'S'},
        {"verbose",             no_argument,       NULL, 'v'},
        {"version",             no_argument,       NULL, 'V'},
        {"help",                no_argument,       NULL, 'h'},
        {"helpshort",           no_argument,       NULL, 'u'},
        {0, 0, 0, 0}
    }; // struct long_options

    // parse command-line options
    int oc = -1;
    int oi = -1;
    while ((oc = getopt_long(argc, argv, "i:o:4n:c:g:w:l:s:fbvh", long_options, &oi)) != -1) {
        switch (oc) {
            case 'i':
                listfile = optarg;
                break;

            case 'o':
                if (optarg [0] == '/') {
                    params.outputDir = optarg;
                } else {
                    params.outputDir  = cwd;
                    params.outputDir += '/';
                    params.outputDir += optarg;
                }
                break;

            case '4':
                judge = 4;
                break;
                
            case 'n':
                params.iterNum = atoi(optarg);
                break;
                
            case 'c':
                params.tissueWeight[0] = atof(optarg);
                break;
                
            case 'g':
                params.tissueWeight[1] = atof(optarg);
                break;
                
            case 'w':
                params.tissueWeight[2] = atof(optarg);
                break;

            case '/':
                params.labelCSF = static_cast<unsigned char>(atoi(optarg));
                break;

            case '|':
                params.labelGM = static_cast<unsigned char>(atoi(optarg));
                break;

            case '_':
                params.labelWM = static_cast<unsigned char>(atoi(optarg));
                break;

            case 'l':
                params.lambda = atof(optarg);
                break;

            case 's':
                if (oi >= 0) {
                    if (strncmp(long_options[oi].name, "fuzzy-suffix", 13) == 0) {
                        suffixMF = optarg;
                    } else if (strncmp(long_options[oi].name, "bias-correct-suffix", 20) == 0) {
                        suffixBC = optarg;
                    } else if (strncmp(long_options[oi].name, "suffix", 20) == 0) {
                        suffix = optarg;
                    } else {
                      std::cerr << "Invalid long option associated with 's': " << long_options[oi].name << "\n";
                      abort();
                    }
                } else {
                    suffix = optarg;
                }
                break;

            case 'f':
                params.outputMemberships = true;
                break;

            case 'b':
            case 'y':
                params.outputBiasCorrections = true;
                if (oc == 'y') params.outputSegments = false;
                break;

            case 'F':
                if      (strncmp(long_options[oi].name, "n2",  5) == 0) extension = ".hdr";
                else if (strncmp(long_options[oi].name, "n1",  5) == 0) extension = ".nii";
                else if (strncmp(long_options[oi].name, "zn2", 6) == 0) extension = ".hdr.gz";
                else if (strncmp(long_options[oi].name, "zn1", 6) == 0) extension = ".nii.gz";
                else {
                  std::cerr << "Invalid long option associated with 'F': " << long_options[oi].name << "\n";
                  abort();
                }
                break;

            case 'S':
                params.seed = static_cast<unsigned int>(time(NULL));
                break;

            case 'v':
                verbosity++;
                break;
                
            case 'V':
                std::cout << "Version: " << std::string(PROJECT_VERSION) << "\n";
                free(cwd);
                exit(EXIT_SUCCESS);
                
            case 'h':
                print_help();
                free(cwd);
                exit(EXIT_SUCCESS);
                break;
                
            case 'u':
                print_usage();
                free(cwd);
                exit(EXIT_SUCCESS);

            case '?':
                // error message already printed
                free(cwd);
                exit(EXIT_FAILURE);

            default:
                // we should never get here
                fprintf(stderr, "See --help for usage information.\n");
                free(cwd);
                exit(EXIT_FAILURE);
        }
    }

    argc -= optind;
    argv += optind;

    // if common extension selected, append it to suffixes
    if (params.outputSegments) {
        params.suffix = get_suffix_with_extension(suffix, extension);
    }
    if (params.outputBiasCorrections) {
        params.suffixBC = get_suffix_with_extension(suffixBC, extension);
    }
    if (params.outputMemberships) {
        params.suffixMF = get_suffix_with_extension(suffixMF, extension);
    }

    // if suffix for memberships does not contain any label placeholder, add one
    if (replace_label(params.suffixMF.c_str(), "") == params.suffixMF) {
        string::size_type pos = params.suffixMF.find('.');
        if (pos == string::npos) params.suffixMF += "_<label>";
        else params.suffixMF.insert(pos, "_<label>");
    }

    // use current working directory as default output directory
    if (params.outputDir.empty()) params.outputDir = cwd;

    // ensure that tissue weights are not negative
    if (params.tissueWeight[0] < 0.0) params.tissueWeight[0] = 0.0;
    if (params.tissueWeight[1] < 0.0) params.tissueWeight[1] = 0.0;
    if (params.tissueWeight[2] < 0.0) params.tissueWeight[2] = 0.0;

    // ensure that lambda is not negative
    if (params.lambda < 0.0) params.lambda = 0.0;

    // -----------------------------------------------------------------------
    // check arguments and read input list file
    if (listfile != NULL && argc > 0) {
        fprintf(stderr, "Either specify an input list file or image files on the command-line. Not both!\n");
        free(cwd);
        exit(EXIT_FAILURE);
    }

    if (judge == 4 && params.outputBiasCorrections) {
        fprintf(stderr, "This program does not yet support the output of the bias corrected images in case of a temporally consistent segmentation.\n");
        free(cwd);
        exit(EXIT_FAILURE);
    }

    // read file names from list file
    char** filenames      = NULL;
    int    number_of_images = 0;

    if (listfile) {
        number_of_images = read_image_file_names(listfile, &filenames);
 
        if (number_of_images < 0) {
            fprintf(stderr, "Failed to read input list file!\n");
            free(cwd);
            exit(EXIT_FAILURE);
        }
    } else {
        number_of_images = argc;
        filenames = reinterpret_cast<char**>(calloc(number_of_images, sizeof(char*)));
        for (int i = 0; i < number_of_images; i++) {
            if (argv[i][0] == '/') {
                filenames[i] = strdup(argv[i]);
            } else {
                int len = strlen (argv[i]);
                filenames [i] = reinterpret_cast<char*>(calloc(lencwd + len + 2, sizeof(char)));
                strncpy(filenames[i], cwd,     lencwd);
                strncat(filenames[i], "/",     1);
                strncat(filenames[i], argv[i], len);
            }
        }
    }

    if (number_of_images == 0) {
        fprintf(stderr, "No input images specified!\n");
        free(cwd);
        exit(EXIT_FAILURE);
    }
 
    // -----------------------------------------------------------------------
    // display arguments of parameters
    if (verbosity > 1) {
        printf("\n");
        printf("Number of images:                %d\n", number_of_images);
        if (listfile) {
            printf("Input list:                      %s\n", listfile);
        }
        printf("Output directory:                %s\n", params.outputDir.c_str());
        printf("Output label map");
        if (number_of_images > 1) printf("s:");
        else                    printf(": ");
        printf("               %s\n", (params.outputSegments ? "yes" : "no"));
        printf("Output fuzzy segmentations:");
        printf("      %s\n", (params.outputMemberships ? "yes" : "no"));
        printf("Output bias corrected image");
        if (number_of_images > 1) printf("s:");
        else                    printf(": ");
        printf("    %s\n", (params.outputBiasCorrections ? "yes" : "no"));
        if (params.outputSegments) {
            printf("Suffix of label map");
            if (number_of_images > 1) printf("s:");
            else                    printf(": ");
            printf("            %s\n", params.suffix.c_str());
        }
        if (params.outputMemberships) {
            printf("Suffix of segmentations:         %s\n", params.suffixMF.c_str());
        }
        if (params.outputBiasCorrections) {
            printf("Suffix of bias corrected image");
            if (number_of_images > 1) printf("s:");
            else                    printf(": ");
            printf(" %s\n", params.suffixBC.c_str());
        }
        printf("Number of iterations:            %d\n", params.iterNum);
        printf("Temporal consistent segments:    %s\n", ((judge == 4) ? "yes" : "no"));
        printf("Weight for CSF:                  %f\n", params.tissueWeight[0]);
        printf("Weight for GM:                   %f\n", params.tissueWeight[1]);
        printf("Weight for WM:                   %f\n", params.tissueWeight[2]);
        printf("Label  for CSF:                  %d\n", static_cast<int>(params.labelCSF));
        printf("Label  for GM:                   %d\n", static_cast<int>(params.labelGM));
        printf("Label  for WM:                   %d\n", static_cast<int>(params.labelWM));
        printf("Lambda:                          %f\n", params.lambda);
        fflush(stdout);
    }

    // -----------------------------------------------------------------------
    // create output directory (if non-existent)
    cbica::createDir(params.outputDir);
 
    // -----------------------------------------------------------------------
    // perform segmentation
    int retval = 0;
 
    if (judge == 3) {
        retval = main3D(filenames, number_of_images, params, verbosity);
    } else {
        retval = main4D(filenames, number_of_images, params, verbosity);
    }

    if (verbosity > 0) printf("\n");

    // -----------------------------------------------------------------------
    // done
    for (int i = 0; i < number_of_images; i++) free(filenames[i]);
    free(filenames);
    free(cwd);

    return retval;
}
