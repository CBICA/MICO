/**
 * @file  io.h
 * @brief I/O functions.
 *
 * Copyright (c) 2011, 2012 University of Pennsylvania. All rights reserved.<br />
 * See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#pragma once
#ifndef _MICO_IO_H
#define _MICO_IO_H


namespace mico {


/**
 * @brief Get length of line excluding trailing newline.
 *
 * @param [in] line Pointer to line start as returned by fgetline ().
 * @param [in] len  Line length as returned by fgetline ().
 *
 * @returns Length of line without trailing "\r\n" or "\n".
 */
int length_of_image_name(const char* line, int len);

/**
 * @brief Get number of images listed in images list file.
 *
 * @param [in] listfile Name of the images list file.
 *
 * @returns Number of images listed in the images list file.
 *
 * @retval -1 Failed to open list file.
 * @retval -2 Error while reading list file.
 */
int get_number_of_images(const char* listfile);

/**
 * @brief Read image file names given in list file.
 *
 * Example:
 * @code
 * char** filenames = NULL;
 * int number_of_images = read_image_file_names("images.lst", &filenames);
 * @endcode
 *
 * @param [in]  listfile  Name of the input list file.
 * @param [out] filenames Array of absolute image file names.
 *                        The memory is allocated via malloc() and
 *                        has to be deallocated by the caller using free().
 *
 * @return Number of images or negative error code.
 *
 * @retval -1 Failed to open list file.
 * @retval -2 Error while reading list file.
 * @retval -3 Memory allocation failed.
 */
int read_image_file_names(const char* listfile, char*** filenames);


} // namespace mico


#endif // _MICO_IO_H_
