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
  ELSEIF(LINUX)
    SET(ZLIB_LIBRARY "${PROJECT_BINARY_DIR}/ep/lib/zlibstatic.a" CACHE STRING "ZLIB location" FORCE)
  ELSE() # Apple has not been implemented/tested
    MESSAGE(FATAL "macOS has not been implemented or tested. Please use either a Windows/Linux machine")
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
      -DBUILD_SHARED_LIBS:BOOL=OFF
  )

 RETURN()
ENDIF()

# at this point, the dependencies are ready
FILE( GLOB_RECURSE ALL_LIBRARIES "${PROJECT_BINARY_DIR}/ep/lib/*" )
INCLUDE_DIRECTORIES(
  ${PROJECT_BINARY_DIR}/ep/include
  ${PROJECT_BINARY_DIR}/ep/include/nifti
  include
)

ADD_LIBRARY(micolib STATIC src/utilities.cxx src/cross-sectional.cxx src/longitudinal.cxx src/io.cxx)
SET_TARGET_PROPERTIES(micolib PROPERTIES OUTPUT_NAME mico)
TARGET_LINK_LIBRARIES(micolib ${ALL_LIBRARIES})

ADD_EXECUTABLE(mico src/mico.cxx)
TARGET_LINK_LIBRARIES(mico micolib)
