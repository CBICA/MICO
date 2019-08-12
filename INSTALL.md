## Intrduction

This document contains the build and installation instructions.

For general build and installation instructions which apply to any software
developed on top of the SBIA Build system And Software Implementation
Standard (BASIS) [1], please refer to the respective how-to guide of the
BASIS documentation [2].



## Obtaining a Copy of the Software

Please see the corresponding section of the BASIS how-to guide [2].


## Installing a Binary Distribution Package

Please see the corresponding section of the BASIS how-to guide [2].


## Building from sources

### Dependencies

The following software has to be installed.

- ZLIB: On Ubuntu, the command is "sudo apt install zlib1g-dev"
- CMake 2.8.12 or higher
- C++ compiler: GCC 4.8+

NOTE: Windows machines are currently not supported.

### Steps

The common steps to build, test, and install software based on CMake,
including this software, are as follows:

1. Extract source files.
2. Create build directory and change to it.
3. Run CMake to configure the build tree - the first round configures the build for NIFTICLIB
3a. Build NIFTICLIB
4a. Run CMake in the same build location once again - this will configure MICO
4b. Build the software using selected build tool.
5. Test the built software.
6. Install the built files.

On Unix-like systems with GNU Make as build tool, these build steps can be
summarized by the following sequence of commands executed in a shell,
where $package and $version are shell variables which represent the name
of this package and the obtained version of the software.

```bash
tar xzf $package-$version-source.tar.gz
mkdir $package-$version-build
cd $package-$version-build
cmake ../$package-$version-source # the first round is for nifticlib
make
cmake ../$package-$version-source # the second round is for mico
make
```

- Press 'c' to configure the build system and 'e' to ignore warnings.
- Set CMAKE_INSTALL_PREFIX and other CMake variables and options.
- Continue pressing 'c' until the option 'g' is available.
- Then press 'g' to generate the configuration files for GNU Make.

$ make
$ make test(optional)
$ make install (optional)

An exhaustive list of minimum build dependencies, including the build tools
along detailed step-by-step build, test, and installation instructions can
be found in the corresponding "Building the Software from Sources" section
of the BASIS how-to guide on software installation [2].

Please refer to this guide first if you are uncertain about above steps or
have problems to build, test, or install the software on your system.
If this guide does not help you resolve the issue, please contact us at
software@cbica.upenn.edu. In case of failing tests, please attach
the output of the following command to your email:
 
$ ctest -V >& test.log

In the following, only package-specific CMake settings available to
configure the build and installation of this software are documented.
