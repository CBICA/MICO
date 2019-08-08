##############################################################################
# @file  CMakeLists.txt
# @brief Root build configuration file.
#
# Copyright (c) 2011 University of Pennsylvania. All rights reserved.<br />
# See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
#
# Contact: SBIA Group <sbia-software at uphs.upenn.edu>
##############################################################################

# ----------------------------------------------------------------------------
# minimum required CMake version
cmake_minimum_required (VERSION 2.8.4)

PROJECT( MICO )

SET(CMAKE_CONFIGURATION_TYPES "Release" CACHE STRING "Default configuration types" FORCE)

IF(NOT EXISTS "${PROJECT_BINARY_DIR}/ep/share/pkgconfig/zlib.pc")

  PROJECT(${PROJECT_NAME}_Depends)

  IF( CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT )
    SET( CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR}/ep CACHE PATH "Changed the default install path" FORCE )
  ENDIF( )
  INCLUDE(ExternalProject)
  ExternalProject_Add( 
    ZLIB
    UPDATE_COMMAND ""
    PATCH_COMMAND ""
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/zlib
    BINARY_DIR thirdparty/zlib_build
    UPDATE_COMMAND ""
    PATCH_COMMAND ""
    #INSTALL_COMMAND ""
    #INSTALL_COMMAND make install
    # INSTALL_COMMAND make install/strip
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
      ${ep_common_args}   
      -DCMAKE_DEBUG_POSTFIX:STRING=d
      -DBUILD_TESTING:BOOL=OFF # testing the ITK build is not required
      #-DCMAKE_INSTALL_PREFIX:PATH=${PROJECT_BINARY_DIR}/ep
  )

  SET(ZLIB_INCLUDE_DIR "${PROJECT_BINARY_DIR}/ep/include" CACHE STRING "ZLIB incude dir" FORCE)
  IF(MSVC)
    SET(ZLIB_LIBRARY "${PROJECT_BINARY_DIR}/ep/lib/zlibstatic.lib" CACHE STRING "ZLIB location" FORCE)
  ELSE()
    SET(ZLIB_LIBRARY "${PROJECT_BINARY_DIR}/ep/lib/zlibstatic.a" CACHE STRING "ZLIB location" FORCE)
  ENDIF()
  
  ExternalProject_Add( 
    NIFTICLIB
    DEPENDS ZLIB
    UPDATE_COMMAND ""
    PATCH_COMMAND ""
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/nifticlib
    BINARY_DIR thirdparty/nifticlib_build
    UPDATE_COMMAND ""
    PATCH_COMMAND ""
    #INSTALL_COMMAND ""
    #INSTALL_COMMAND make install
    # INSTALL_COMMAND make install/strip
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
      ${ep_common_args}   
      -DCMAKE_DEBUG_POSTFIX:STRING=d
      -DBUILD_TESTING:BOOL=OFF # testing the ITK build is not required
      #-DCMAKE_INSTALL_PREFIX:PATH=${PROJECT_BINARY_DIR}/ep
      -DZLIB_LIBRARY=${ZLIB_LIBRARY}
      -DZLIB_INCLUDE_DIR=${ZLIB_INCLUDE_DIR}
  )

 RETURN()
ENDIF()

# ----------------------------------------------------------------------------
# version of BASIS template used by this project
set (TEMPLATE_VERSION 1.5)

# ----------------------------------------------------------------------------
# include BASIS policies, settings, macros, and functions

# circumvent issue with CMake's find_package() interpreting these variables
# relative to the current binary directory instead of the top-level directory
if (BASIS_DIR AND NOT IS_ABSOLUTE "${BASIS_DIR}")
  set (BASIS_DIR "${CMAKE_BINARY_DIR}/${BASIS_DIR}")
  get_filename_component (BASIS_DIR "${BASIS_DIR}" ABSOLUTE)
endif ()
# moreover, users tend to specify the installation prefix instead of the
# actual directory containing the package configuration file
if (IS_DIRECTORY "${BASIS_DIR}")
  list (INSERT CMAKE_PREFIX_PATH 0 "${BASIS_DIR}")
endif ()
# if BASIS is not found, set CMAKE_INSTALL_PREFIX to invalid value if not
# set explicitly at the command-line such that BASIS will still initialize
# it to its defaults even though already cached
find_package (BASIS)
if (NOT BASIS_FOUND)
  if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "" CACHE INTERNAL "Installation prefix." FORCE)
  endif ()
  message (FATAL_ERROR "BASIS not found! Please specify installation directory"
                       " of BASIS or the location of the BASISConfig.cmake file"
                       " using the BASIS_DIR variable.")
endif ()

# ----------------------------------------------------------------------------
# configure build system
basis_project_impl ()
