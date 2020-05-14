############################################################################
#
#  Program:         ScaLAPACK
#
#  Module:          SLmake.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   March 20, 1995
#
#  Modified:
#
#  Send bug reports, comments or suggestions to scalapack@cs.utk.edu
#
############################################################################
#
SHELL         = /bin/sh
#
#  The complete path to the top level of ScaLAPACK directory, usually
#  $(HOME)/SCALAPACK
#
home          = /usr/local/SCALAPACK
#
#  The platform identifier to suffix to the end of library names
#
PLAT          = LINUX
#
LAPACKdir     = $(HOME)/LAPACK
#
#  BLACS setup.  All version need the debug level (0 or 1),
#  and the directory where the BLACS libraries are
#
BLACSDBGLVL   = 0
BLACSdir      = $(HOME)/BLACS/LIB
#
#  MPI setup; uncomment and tailor to your system if using MPIBLACS
#  Will need to comment out the default native BLACS setup below below
#
USEMPI        = -DUsingMpiBlacs
SMPLIB        = /src/icl/MPI/mpich1.1/lib/LINUX/ch_p4/libmpi.a
BLACSFINIT    = $(BLACSdir)/blacsF77init_MPI-$(PLAT)-$(BLACSDBGLVL).a
BLACSCINIT    = $(BLACSdir)/blacsCinit_MPI-$(PLAT)-$(BLACSDBGLVL).a
BLACSLIB      = $(BLACSdir)/blacs_MPI-$(PLAT)-$(BLACSDBGLVL).a
TESTINGdir    = $(home)/TESTING

#
#  PVMBLACS setup, comment out next 6 lines if using MPI
#
#USEMPI        =
#SMPLIB        = $(PVM_ROOT)/lib/$(PLAT)/libpvm3.a
#BLACSFINIT    =
#BLACSCINIT    =
#BLACSLIB      = $(BLACSdir)/blacs_PVM-$(PLAT)-$(BLACSDBGLVL).a
#TESTINGdir    = $(HOME)/pvm3/bin/$(PLAT)

CBLACSLIB     = $(BLACSCINIT) $(BLACSLIB) $(BLACSCINIT)
FBLACSLIB     = $(BLACSFINIT) $(BLACSLIB) $(BLACSFINIT)

BLASLIB       = $(LAPACKdir)/blas_$(PLAT).a
#
#  The directories to find the various pieces of ScaLapack
#
PBLASdir      = $(home)/PBLAS
SRCdir        = $(home)/SRC
TESTdir       = $(home)/TESTING
PBLASTSTdir   = $(TESTINGdir)
TOOLSdir      = $(home)/TOOLS
REDISTdir     = $(home)/REDIST
REDISTTSTdir  = $(TESTINGdir)
#
#  The fortran and C compilers, loaders, and their flags
#
F77           = g77
CC            = gcc
NOOPT         =
F77FLAGS      = -O3 $(NOOPT)
F77LOADER     = $(F77)
F77LOADFLAGS  = $(F77FLAGS)
#
#  C preprocessor defs for compilation 
#  (-DNoChange, -DAdd_, -DUpCase, or -Df77IsF2C)
#
CDEFS         = -Df77IsF2C -DNO_IEEE $(USEMPI)
#
#  The archiver and the flag(s) to use when building archive (library)
#  Also the ranlib routine.  If your system has no ranlib, set RANLIB = echo
#
ARCH          = ar
ARCHFLAGS     = cr
RANLIB        = ranlib
#
#  The name of the libraries to be created/linked to
#
PBLASLIB      = $(home)/pblas_$(PLAT).a
SCALAPACKLIB  = $(home)/scalapack_$(PLAT).a
TOOLSLIB      = $(home)/tools_$(PLAT).a
REDISTLIB     = $(home)/redist_$(PLAT).a
#
PBLIBS        = $(PBLASLIB) $(TOOLSLIB) $(FBLACSLIB) $(BLASLIB) $(SMPLIB)
PRLIBS        = $(REDISTLIB) $(TOOLSLIB) $(CBLACSLIB) $(SMPLIB)
RLIBS         = $(SCALAPACKLIB) $(PBLASLIB) $(REDISTLIB) $(TOOLSLIB) \
                $(FBLACSLIB) $(CBLACSLIB) $(BLASLIB) $(SMPLIB)
LIBS          = $(SCALAPACKLIB) $(PBLIBS)
