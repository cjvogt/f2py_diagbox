# -------------------------------------------------------------------------
# FILE    : makefile
# AUTHOR  : Nicolas Gruber
# LAST REV: 05-AUG-96  adapted for program diag_boxmodel
#           17-SEP-99  adapted for SUN and interannual variability version
#                      special version for chem library part
#           14-NOV-01  adapted for lf95
#
# -------------------------------------------------------------------------
#
# ----------------------------------------------------------------------
# -- Absolute path to the root directory of source code tree
# ----------------------------------------------------------------------
#

ROOT = /home/prafter/f2py-diagbox/src_ia

#
# ----------------------------------------------------------------------
# -- Colon separated list of directories which includes include files
#    and modules.
# ----------------------------------------------------------------------
#

DIRS = .:$(ROOT)/incl

#
# ----------------------------------------------------------------------
# -- Colon separated list of directories containing source code.
#    Directories will be searched in the order given
# ----------------------------------------------------------------------
#

SRC_DIRS = .

#
# ----------------------------------------------------------------------
# -- File containing list of CPP define directives
#    (default: $(ROOT)/ifdefs)
# ----------------------------------------------------------------------
#

include $(ROOT)/ifdefs

# -------------------------------------------------------------------------
# -- COMPILATION MACROS 
# -------------------------------------------------------------------------

F_C      = gfortran
#F_C      = lf95
#
#F_FLAGS  =  
F_FLAGS  = -C -g  -fbacktrace -O0 -finit-real=snan -finit-integer=-2147483648 -fno-automatic -std=legacy -I/home/prafter/f2py-diagbox/src_ia/incl -fcheck=all -finit-real=nan -ffpe-trap=invalid,zero,overflow -Wall -Wextra -Wno-compare-reals -Wno-unused-label -Wno-unused-variable -Wno-unused-dummy-argument -Wno-unused-parameter -Wno-line-truncation -Wno-cpp

LFLAGS = -s

LIBNAME = libchem.a

L = $(LIBNAME)

# -------------------------------------------------------------------------
# -- LOCATION OF INCLUDE FILES: (for directives like 'include filename')  
# -- example: ( I_PATH = /users/user_name/proj/include )  
# -------------------------------------------------------------------------

I_PATH = \
	 $(ROOT)/incl

O_FILES = \
	vp_purewater.o\
	vp_seawater.o\
	calcco2.o\
	chem.o\
	estimate.o\
	schmidt_number.o
#
H_FILES = \
 
#
# -------------------------------------------------------------------------
# -- LIBRARY OBJECTS
# -------------------------------------------------------------------------
#
LIBOBJ = \
	$(L)(vp_purewater.o)\
	$(L)(vp_seawater.o)\
	$(L)(calcco2.o)\
	$(L)(chem.o)\
	$(L)(estimate.o)\
	$(L)(schmidt_number.o)

#
# -------------------------------------------------------------------------
# -- SET Default
# -------------------------------------------------------------------------
#
default: all
#
# -------------------------------------------------------------------------
# -- MAIN DEPENDENCY for library
# -------------------------------------------------------------------------
#
$(LIBNAME) : $(LIBOBJ)
	@echo ""
	@echo " *****************************"
	@echo " *** running ranlib        ***"
	@echo " *****************************"
	@echo ""
	ranlib $(LIBNAME)
	@echo ""
	@echo " *****************************"
	@echo " *** lib is now up to date ***"
	@echo " *****************************"
	@echo ""
library: all

all: FIRST_ACTION \
     $(LIBNAME) \
     LAST_ACTION

clean:
	@echo ""
	@echo " ***************************************"
	@echo " *** removing object and *.for files ***"
	@echo " ***************************************"
	@echo ""
	rm -f *.o *% *.for core
#
small:
	@echo ""
	@echo " ***************************************"
	@echo " *** compressing .for and .r files     ***"
	@echo " ***************************************"
	@echo ""
	compress *.r *.for
#
big:
	@echo ""
	@echo " ***************************************"
	@echo " *** uncompressing .Z files"
	@echo " ***************************************"
	@echo ""
	uncompress *.Z
#
FIRST_ACTION:
	@echo ""
	@echo " ***************************************"
	@echo " *** UPDATING $(L) ***"
	@echo " ***************************************"
	@echo ""
 
LAST_ACTION:
	@echo ""
	@echo " *************"
	@echo " *** DONE. ***"
	@echo " *************"

	@echo ""
# -------------------------------------------------------------------------
# -- Default Rules and Suffixes
# -------------------------------------------------------------------------

.SUFFIXES:
.SUFFIXES: .a .o .for .F

.F.for: 
	@echo ""
	@echo " *** PREPROCESSING $<... ***"
	$(F_C) $(D_FLAGS) $(F_FLAGS) -I$(I_PATH) -E $< > $*.for 
#	$(F_C) $(D_FLAGS) $(F_FLAGS) -I$(I_PATH) -P $< -o $*.for

.for.o:   
	@echo ""
	@echo " *** COMPILING $<... ***"
	$(F_C) -c $(F_FLAGS) $*.for -o $*.o

.o.a:
	@echo ""
	@echo " *** ARCHIVING $<... ***"
	ar rv $(L) $*.o
#	rm -f $*.o


# -------------------------------------------------------------------------
# -- Dependencies (explicit) 
#          if you have include files in a subdirectory they must be 
#          indicated together with the path
# -------------------------------------------------------------------------
#
vp_purewater.for : vp_purewater.F
vp_purewater.o : vp_purewater.for
#
vp_seawater.for : vp_seawater.F
vp_seawater.o : vp_seawater.for
#
schmidt_number.for : schmidt_number.F
schmidt_number.o : schmidt_number.for
#
chem.for : chem.F
chem.o : chem.for
#
calcco2.for : calcco2.F
calcco2.o : calcco2.for
#
estimate.for : estimate.F
estimate.o : estimate.for
#
$(L) : $(O_FILES)
