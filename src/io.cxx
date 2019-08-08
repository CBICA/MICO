/**
 * @file  io.cxx
 * @brief I/O functions.
 *
 * Copyright (c) 2011, 2012 University of Pennsylvania. All rights reserved.<br />
 * See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#include <string>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

//#include <basis/os/path.h> // dirname(), normpath()

#include "mico/io.h"


// acceptable in .cxx file
using namespace std;
//using namespace basis;


namespace mico {


// ---------------------------------------------------------------------------
// bool _getline(FILE* fstream, char* line, size_t len)
#if (defined (__APPLE__) || defined (__OSX__))
#  define _getline(fp, line, len) \
    ((line = fgetln(fp, &len)) != NULL)
#else
ssize_t g_dummy_variable_used_by_getline = 0;
#  define _getline(fp, line, len) \
    ((g_dummy_variable_used_by_getline = getline(&line, &len, fp)) != -1 && \
      (len = /*yes, assignment!*/ static_cast<size_t>(g_dummy_variable_used_by_getline)))
#endif


// ---------------------------------------------------------------------------
int length_of_image_name(const char* line, int len)
{
    // line ends in "\r\n" as on Windows => actual length is two less
    if (len > 1 && line[len - 2] == '\r' && line[len - 1] == '\n') return len - 2;
    // line ends in '\n' => actual length is one less
    if (len > 0 && line[len - 1] == '\n') return len - 1;
    // given length is the actual length of the line
    return len;
}

// ---------------------------------------------------------------------------
int get_number_of_images(const char* listfile)
{
    FILE*   fp    = NULL; // file pointer of images list
    char*   line  = NULL; // pointer to start of line in fp
    size_t  len   = 0;    // length of the line read
    int     num   = 0;    // number of images

    fp = fopen(listfile, "r");

    if (fp == NULL) {
        return -1; // indicates list file open error
    }

    while (_getline(fp, line, len))
    {
        if (length_of_image_name(line, len) > 0) num++;
    }

    if (ferror(fp) != 0) {
        num = -2; // indicates error while reading list file
    }

    fclose(fp);
    free(line);

    return num;
}

// ---------------------------------------------------------------------------
int read_image_file_names(const char* listfile, char*** filenames)
{
    char*  listpath = NULL; // non-const copy of input listfile
    FILE*  fp       = NULL; // file pointer of images list
    char*  line     = NULL; // pointer to start of line in fp
    size_t len      = 0;    // length of line
    int    i        = 0;    // loop variable
    int    n        = 0;    // number of image names read
    char*  hname    = NULL; // name of image file
    int    hlen     = 0;    // size of hname buffer - 1

    // get number of images
    int num = get_number_of_images(listfile);
    if (num < 1) return num;

    // initialize output array
    *filenames = reinterpret_cast<char**>(calloc(num, sizeof(char**)));
    if (*filenames == NULL) return -3; // indicates memory allocation error

    // open list file
    fp = fopen(listfile, "r");
    if (fp == NULL) return -1; // indicates list file open error

    // directory of list file
    string root    = os::path::dirname(listfile);
    int    lenroot = static_cast<int>(root.size());

    // read list line by line
    while (_getline(fp, line, len)) {
        // skip empty lines
        len = length_of_image_name(line, len);
        if (len == 0) continue;
        // allocate filename buffer
        if (line[0] != '/') {
            hlen = lenroot + len + 1;
        } else {
            hlen = len;
        }
        hname = (char*)malloc((hlen + 1) * sizeof(char));
        if (hname == NULL) {
            n = -3; // indicates memory allocation error
            break;
        }
        // copy line into filename buffer and terminate string
        if (line[0] != '/') {
            strncpy(hname, root.c_str(), lenroot + 1);
            strncat(hname, "/", 1);
            strncat(hname, line, len);
        } else {
            strncpy(hname, line, len);
        }
        hname[hlen] = '\0';
        // normalize file path
        string tmp = os::path::normpath(hname);
        strncpy(hname, tmp.c_str(), hlen);
        // add file name to output list
        (*filenames)[n] = hname;
        hname = NULL;
        // increment number of image names read
        n++;
    }

    if (ferror(fp) != 0) {
        n = -2; // indicates list file read error
    }

    // close file
    fclose(fp);

    // clean up
    if (n < 0) {
        for (i = 0; i < num; i++) {
            if ((*filenames)[i]) {
                free((*filenames)[i]);
                (*filenames)[i] = NULL;
            }
        }
        free(*filenames);
        *filenames = NULL;
    }

    free(hname);
    free(line);
    free(listpath);

    return n;
}


} // namespace mico
