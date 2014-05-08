MMPBSA
======

About
-----

This program organizes the simulation of Molecule Mechanics, solves the Poisson Bolzmann equation, and calculates surface area of molecules.

MMBPSA was developed as part of work done in the Bashford Group in the Structural Biology Department at St Jude Children's Research Hospital. Electrostatic calculations are performed using the MEAD library, which may be found at http://stjuderesearch.org/site/lab/bashford

Molecular dynamics may be performed by a modified version of Sander which must be stored in the path of MMPBSA as "moldyn". There is patch file "sander_mmpbsa.patch" which should be used to modify sander to work with MMPBSA. Simply run "patch < sander_mmpbsa.patch" in the Sander source directory and rebuild sander.

Poisson-Boltzmann and surface area solving may be done using sander or gromacs trajectory and topology files.

MMPBSA may be built to run on a BOINC grid (without requiring further modification to Sander). Simply run configure with the option "--with-boinc=<dir>", where <dir> is the path to the BOINC API install. A separate graphics application, mmpbsa_graphics, has been made to produce graphics during calculation in BOINC. This can be compiled optionally with the flag "--with-graphics=<dir>", where <dir> is an optional path to GL header files.

License
-------

MMPBSA is released under the terms of the GNU General Public License version 3. For details, see COPYING or http://www.gnu.org/licenses/gpl.html.

Install
-------

MMPBSA has been setup to use GNU build tools. Configure and install-sh files are provided. However, if these do not work, one can run "./setup" to recreate those files based on the specific platform, if autotools is installed. Instructions to use configure and make are provided below. For instructions on using GNU build tools, see http://www.gnu.org/software/automake/ .

Building MMPBSA will create a library with all of the objects and functions needed for mmpbsa functions. "make install" will install this library in the directory tree specified by PREFIX provided to configure. Additionally, a program is built to perform the calculations.

Multithreading was added in version 0.10, using pthreads. With version 0.10 also came the introduction of gromacs file format support. To use gromacs file formats, MMPBSA must be build with the gromacs libraries and headers. Building with gromacs libraries is indicated to configure with the flag "--with-gromacs". Currently, the config.h in the gromacs source direction must be placed in the gromacs include directory used by MMPBSA. This will be changed in the future, as it is an undesirable, temporary solution.

If you get a linker error that complains about missing pthread, run configure wiht the flag --enable-static


Build Status
------------

[![Build Status](https://travis-ci.org/kd0kfo/mmpbsa.svg?branch=master)](https://travis-ci.org/kd0kfo/mmpbsa)
