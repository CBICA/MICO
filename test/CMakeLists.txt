##############################################################################
# @file  CMakeLists.txt
# @brief Build configuration of software testing.
#
# Unit tests test only single functions, classes, and modules of the software.
# System tests, on the other side, execute the programs of this package with
# given test input data and compare the output of the programs to the expected
# results. For the definition of system tests which perform image regression
# tests, i.e., compare the images written by the test to given baseline images,
# the basis.testdriver executable target imported from BASIS can be used to
# execute the executable to test and perform the regression tests.
#
# Example:
# @code
# basis_add_test (
#  COMMAND basis.testdriver -v
#    --intensity-tolerance 2
#    --compare "t1_bc.nii.gz"
#              "${PROJECT_TESTING_DIR}/baseline/t1_bc.nii.gz"
#    --
#    $<TARGET_FILE:biascorrect>
#      "${PROJECT_TESTING_DIR}/input/t1.nii.gz"
#      t1_bc.nii.gz
# )
# @endcode
#
# If the tests are implemented in C/C++ as supported by CMake's
# create_test_sourcelist() function, a custom test driver can be generated
# using the basis_add_test_driver() function and this test driver then be used
# to execute the tests and perform the regression tests.
#
# Example:
# @code
# basis_add_test_driver (testdriver test1.cxx test2.cxx test3.cxx)
#
# basis_add_test (
#   COMMAND testdriver -v
#     --intensity-tolerance 2
#     --compare "t1_bc.nii.gz"
#               "${PROJECT_TESTING_DIR}/baseline/t1_bc.nii.gz"
#     --
#     biascorrect # test1.cxx defines int biascorrect(int argc, char* argv[])
#      "${PROJECT_TESTING_DIR}/input/t1.nii.gz"
#      t1_bc.nii.gz
# )
#
# # ...
# @endcode
#
# @sa http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:create_test_sourcelist
#
# The built test executables can be found in @c TESTING_RUNTIME_DIR.
# The test output shall be written to @c TESTING_OUTPUT_DIR. In most cases,
# a subdirectory for each test should be created under this output directory.
# Therefore, the basis_add_test() function creates such subdirectory and sets
# the working directory of the test to this output directory.
#
# Copyright (c) 2011, 2012 University of Pennsylvania. All rights reserved.<br />
# See https://www.med.upenn.edu/sbia/software/license.html or COPYING file.
#
# Contact: SBIA Group <software at cbica.upenn.edu>
##############################################################################

# ============================================================================
# internal tests
# ============================================================================

if (IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/internal")
  add_subdirectory ("internal")
endif ()

# ============================================================================
# public tests
# ============================================================================

set (TIMEOUT 86400) # 24h, i.e., do not use timeout test criterium

# ----------------------------------------------------------------------------
# example - bias correction
basis_add_test (
  biascorrect
  COMMAND
    basis.testdriver
        # enable output of hostname
        --verbose
        # compare bias corrected images
        --intensity-tolerance 5
        --compare "${TESTING_OUTPUT_DIR}/biascorrect/BC01-T1-byte_cbq_bc.nii.gz"
                  "${PROJECT_TESTING_DIR}/baseline/example/3d/BC01-T1-byte_cbq_bc.nii.gz"
        # test command
        -- $<TARGET_FILE:mico>
                --verbose --verbose
                --bias-correct-only
                "${PROJECT_EXAMPLE_DIR}/BC01-T1-byte_cbq.nii.gz"
)

basis_set_tests_properties (biascorrect PROPERTIES TIMEOUT ${TIMEOUT})

# ----------------------------------------------------------------------------
# example - cross-sectional segmentation
basis_add_test (
  segment_cross-sectional
  COMMAND
    basis.testdriver
        # enable output of hostname
        --verbose
        # compare segmentations
        --intensity-tolerance 0
        --max-number-of-differences 5000
        --compare "BC01-T1-byte_cbq_labels.nii.gz"
                  "${PROJECT_TESTING_DIR}/baseline/example/3d/BC01-T1-byte_cbq_labels.nii.gz"
        # test command
        -- $<TARGET_FILE:mico>
                --verbose --verbose
                "${PROJECT_EXAMPLE_DIR}/BC01-T1-byte_cbq.nii.gz"
)

basis_set_tests_properties (segment_cross-sectional PROPERTIES TIMEOUT ${TIMEOUT})

# ----------------------------------------------------------------------------
# example - cross-sectional fuzzy segmentations
basis_add_test (
  segment_cross-sectional_fuzzy
  COMMAND
    basis.testdriver
        # enable output of hostname
        --verbose
        # compare segmentations
        --intensity-tolerance 0
        --max-number-of-differences 5000
        --compare "BC01-T1-byte_cbq_labels.nii.gz"
                  "${PROJECT_TESTING_DIR}/baseline/example/3d/BC01-T1-byte_cbq_labels.nii.gz"
        # compare fuzzy segmentations
        --intensity-tolerance 0.2
        --max-number-of-differences 10
        --compare "BC01-T1-byte_cbq_csf_map.hdr"
                  "${PROJECT_TESTING_DIR}/baseline/example/3d/BC01-T1-byte_cbq_csf.nii.gz"
        --compare "BC01-T1-byte_cbq_gm_map.hdr"
                  "${PROJECT_TESTING_DIR}/baseline/example/3d/BC01-T1-byte_cbq_gm.nii.gz"
        --compare "BC01-T1-byte_cbq_wm_map.hdr"
                  "${PROJECT_TESTING_DIR}/baseline/example/3d/BC01-T1-byte_cbq_wm.nii.gz"
        # test command
        -- $<TARGET_FILE:mico>
                --verbose --verbose
                --fuzzy
                --fuzzy-suffix _<label>_map.hdr
                "${PROJECT_EXAMPLE_DIR}/BC01-T1-byte_cbq.nii.gz"
)

basis_set_tests_properties (segment_cross-sectional_fuzzy PROPERTIES TIMEOUT ${TIMEOUT})

# ----------------------------------------------------------------------------
# example - cross-sectional segmentation, list file input
basis_add_test (
  segment_cross-sectional_list
  COMMAND
    basis.testdriver
        # enable output of hostname
        --verbose
        # compare segmentations
        --intensity-tolerance 0
        --max-number-of-differences 5000
        --compare "BC01-T1-byte_cbq_labels.nii.gz"
                  "${PROJECT_TESTING_DIR}/baseline/example/3d/BC01-T1-byte_cbq_labels.nii.gz"
        --compare "BC02-T1-byte_cbq_flirtedToY16dof_labels.nii.gz"
                  "${PROJECT_TESTING_DIR}/baseline/example/3d/BC02-T1-byte_cbq_flirtedToY16dof_labels.nii.gz"
        --compare "BC03-T1-byte_cbq_flirtedToY16dof_labels.nii.gz"
                  "${PROJECT_TESTING_DIR}/baseline/example/3d/BC03-T1-byte_cbq_flirtedToY16dof_labels.nii.gz"
        
        # test command
        -- $<TARGET_FILE:mico>
                --verbose --verbose
                --inputlist "${CMAKE_CURRENT_LIST_DIR}/exampletest.lst"
)

basis_set_tests_properties (segment_cross-sectional_list PROPERTIES TIMEOUT ${TIMEOUT})

# ----------------------------------------------------------------------------
# example - longitudinal segmentation
basis_add_test (
  segment_longitudinal
  COMMAND
    basis.testdriver
        # enable output of hostname
        --verbose
        # compare segmentations
        --intensity-tolerance 0
        --max-number-of-differences 3500
        --compare "BC01-T1-byte_cbq_labels.nii.gz"
                  "${PROJECT_TESTING_DIR}/baseline/example/4d/BC01-T1-byte_cbq_labels.nii.gz"
        --compare "BC02-T1-byte_cbq_flirtedToY16dof_labels.nii.gz"
                  "${PROJECT_TESTING_DIR}/baseline/example/4d/BC02-T1-byte_cbq_flirtedToY16dof_labels.nii.gz"
        --compare "BC03-T1-byte_cbq_flirtedToY16dof_labels.nii.gz"
                  "${PROJECT_TESTING_DIR}/baseline/example/4d/BC03-T1-byte_cbq_flirtedToY16dof_labels.nii.gz"
        # test command
        -- $<TARGET_FILE:mico>
                --verbose --verbose
                --4d
                "${PROJECT_EXAMPLE_DIR}/BC01-T1-byte_cbq.nii.gz"
                "${PROJECT_EXAMPLE_DIR}/BC02-T1-byte_cbq_flirtedToY16dof.nii.gz"
                "${PROJECT_EXAMPLE_DIR}/BC03-T1-byte_cbq_flirtedToY16dof.nii.gz"
)

basis_set_tests_properties (segment_longitudinal PROPERTIES TIMEOUT ${TIMEOUT})

# ----------------------------------------------------------------------------
# example - longitudinal fuzzy segmentation
basis_add_test (
  segment_longitudinal_fuzzy
  COMMAND
    basis.testdriver
        # enable output of hostname
        --verbose
        # compare segmentations
        --intensity-tolerance 0
        --max-number-of-differences 3500
        --compare "BC01-T1-byte_cbq_labels.nii.gz"
                  "${PROJECT_TESTING_DIR}/baseline/example/4d/BC01-T1-byte_cbq_labels.nii.gz"
        --compare "BC02-T1-byte_cbq_flirtedToY16dof_labels.nii.gz"
                  "${PROJECT_TESTING_DIR}/baseline/example/4d/BC02-T1-byte_cbq_flirtedToY16dof_labels.nii.gz"
        --compare "BC03-T1-byte_cbq_flirtedToY16dof_labels.nii.gz"
                  "${PROJECT_TESTING_DIR}/baseline/example/4d/BC03-T1-byte_cbq_flirtedToY16dof_labels.nii.gz"
        # compare fuzzy segmentations
        --intensity-tolerance 0.2
        --max-number-of-differences 10
        --compare "BC01-T1-byte_cbq_wm.nii.gz"
                  "${PROJECT_TESTING_DIR}/baseline/example/4d/BC01-T1-byte_cbq_wm.nii.gz"
        --compare "BC02-T1-byte_cbq_flirtedToY16dof_gm.nii.gz"
                  "${PROJECT_TESTING_DIR}/baseline/example/4d/BC02-T1-byte_cbq_flirtedToY16dof_gm.nii.gz"
        --compare "BC03-T1-byte_cbq_flirtedToY16dof_csf.nii.gz"
                  "${PROJECT_TESTING_DIR}/baseline/example/4d/BC03-T1-byte_cbq_flirtedToY16dof_csf.nii.gz"
        # test command
        -- $<TARGET_FILE:mico>
                --verbose --verbose
                --4d --fuzzy
                "${PROJECT_EXAMPLE_DIR}/BC01-T1-byte_cbq.nii.gz"
                "${PROJECT_EXAMPLE_DIR}/BC02-T1-byte_cbq_flirtedToY16dof.nii.gz"
                "${PROJECT_EXAMPLE_DIR}/BC03-T1-byte_cbq_flirtedToY16dof.nii.gz"
)

basis_set_tests_properties (segment_longitudinal_fuzzy PROPERTIES TIMEOUT ${TIMEOUT})

# ----------------------------------------------------------------------------
# example - longitudinal segmentation, list file input
basis_add_test (
  segment_longitudinal_list
  COMMAND
    basis.testdriver
        # enable output of hostname
        --verbose
        # compare segmentations
        --intensity-tolerance 0
        --max-number-of-differences 3500
        --compare "BC01-T1-byte_cbq_labels.nii.gz"
                  "${PROJECT_TESTING_DIR}/baseline/example/4d/BC01-T1-byte_cbq_labels.nii.gz"
        --compare "BC02-T1-byte_cbq_flirtedToY16dof_labels.nii.gz"
                  "${PROJECT_TESTING_DIR}/baseline/example/4d/BC02-T1-byte_cbq_flirtedToY16dof_labels.nii.gz"
        --compare "BC03-T1-byte_cbq_flirtedToY16dof_labels.nii.gz"
                  "${PROJECT_TESTING_DIR}/baseline/example/4d/BC03-T1-byte_cbq_flirtedToY16dof_labels.nii.gz"
        # test command
        -- $<TARGET_FILE:mico>
                --verbose --verbose
                --4d
                --inputlist "${CMAKE_CURRENT_LIST_DIR}/exampletest.lst"
)

basis_set_tests_properties (segment_longitudinal_list PROPERTIES TIMEOUT ${TIMEOUT})

# ----------------------------------------------------------------------------
# empty list file
basis_add_test (
  empty_inputlist
  COMMAND
    mico --inputlist "${CMAKE_CURRENT_LIST_DIR}/empty.lst"
)

basis_set_tests_properties (empty_inputlist PROPERTIES WILL_FAIL TRUE)

# ----------------------------------------------------------------------------
# list file with non-existent image
basis_add_test (
  wrong_inputlist
  COMMAND
    mico --inputlist "${CMAKE_CURRENT_LIST_DIR}/wrong.lst"
)

basis_set_tests_properties (wrong_inputlist PROPERTIES WILL_FAIL TRUE)
