# -------------------------------------------------------------------------
# FILE    : makefile
# AUTHOR  : Nicolas Gruber
# LAST REV: 05-AUG-96  adapted for program diag_boxmodel
#           17-SEP-99  adapted for SUN and interannual variability version
#                      special version for plot library part
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
F_FLAGS  = -C -g  -fbacktrace -O0 -finit-real=snan -finit-integer=-2147483648 -fno-automatic -std=legacy -I/home/data/Gruber-Box-Model-2/src_ia/incl -fcheck=all -finit-real=nan -ffpe-trap=invalid,zero,overflow -Wall -Wextra -Wno-compare-reals -Wno-unused-label -Wno-unused-variable -Wno-unused-dummy-argument -Wno-unused-parameter -Wno-line-truncation -Wno-cpp

LFLAGS = -s

LIBNAME = libplot.a
L = $(LIBNAME)

# -------------------------------------------------------------------------
# -- LOCATION OF INCLUDE FILES: (for directives like 'include filename')  
# -- example: ( I_PATH = /users/user_name/proj/include )  
# -------------------------------------------------------------------------

I_PATH = \
	 $(ROOT)/incl

O_FILES = \
	generate_ncar_plots.o\
	plot_results.o\
	generate_mc_ncar_plots.o\
	calc_grid.o\
	bkgroundcolor.o\
	simplecolor.o
#
H_FILES = \
	${I_PATH}/diagboxmod.h     \
	${I_PATH}/simulation.h     \
	${I_PATH}/input.h          \
	${I_PATH}/slabs.h          \
	${I_PATH}/periods.h        \
	${I_PATH}/processes.h      \
	${I_PATH}/averages.h       \
	${I_PATH}/intval.h         \
	${I_PATH}/solution.h       \
	${I_PATH}/ncar.h           \
	${I_PATH}/plot_sets.h      \
	${I_PATH}/obs_values.h     \
	${I_PATH}/statistics.h     \
	${I_PATH}/mc_plot_slabs.h  
 
#
# -------------------------------------------------------------------------
# -- LIBRARY OBJECTS
# -------------------------------------------------------------------------
#
LIBOBJ = \
	$(L)(generate_ncar_plots.o)\
	$(L)(plot_results.o)\
	$(L)(generate_mc_ncar_plots.o)\
	$(L)(calc_grid.o)\
	$(L)(bkgroundcolor.o)\
	$(L)(simplecolor.o)

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
	@echo " *** removing object files ***"
	@echo " ***************************************"
	@echo ""
	rm -f *.o *% core *.for *.a
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
#
generate_ncar_plots.for : generate_ncar_plots.F $(H_FILES)
generate_ncar_plots.o : generate_ncar_plots.for
#
generate_mc_ncar_plots.for : generate_mc_ncar_plots.F $(H_FILES)
generate_mc_ncar_plots.o : generate_mc_ncar_plots.for
#
plot_results.for : plot_results.F $(H_FILES)
plot_results.o : plot_results.for
#
calc_grid.for : calc_grid.F
calc_grid.o : calc_grid.for
#
bkgroundcolor.for     : bkgroundcolor.F
bkgroundcolor.o     : bkgroundcolor.for
#
simplecolor.for : simplecolor.F
simplecolor.o : simplecolor.for
#
$(L) : $(O_FILES)
