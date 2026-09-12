# -------------------------------------------------------------------------
# FILE    : makefile
# AUTHOR  : Nicolas Gruber
# LAST REV: 05-AUG-96  adapted for program diag_boxmodel
#           17-SEP-99  adapted for SUN and interannual variability version
#                      special version for numrecip library part
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
!F_C      = lf95
#
F_FLAGS  = -C -g  -fbacktrace -O0 -finit-real=snan -finit-integer=-2147483648 -fno-automatic -std=legacy -I/home/data/Gruber-Box-Model-2/src_ia/incl -fcheck=all -finit-real=nan -ffpe-trap=invalid,zero,overflow -Wall -Wextra -Wno-compare-reals -Wno-unused-label -Wno-unused-variable -Wno-unused-dummy-argument -Wno-unused-parameter -Wno-line-truncation -Wno-cpp

LFLAGS = -s

LIBNAME = libnumrecip.a
L = $(LIBNAME)

# -------------------------------------------------------------------------
# -- LOCATION OF INCLUDE FILES: (for directives like 'include filename')  
# -- example: ( I_PATH = /users/user_name/proj/include )  
# -------------------------------------------------------------------------

I_PATH = \
	 $(ROOT)/incl

O_FILES = \
	ran1.o\
	func.o\
	dfunc.o\
	f1dim.o\
	df1dim.o\
	frprmn.o\
	linmin.o\
	dbrent.o\
	mnbrak.o\
	cost_function.o\
	calc_hessian.o\
	invert_hessian.o\
	svdcmp.o\
	pythag.o\
	gasdev.o
#
H_FILES = \
	${I_PATH}/diagboxmod.h     \
	${I_PATH}/simulation.h     \
	${I_PATH}/gasex_params.h   \
	${I_PATH}/diffent_params.h \
	${I_PATH}/adv_params.h     \
	${I_PATH}/ncp_params.h     \
	${I_PATH}/c13_params.h     \
	${I_PATH}/harm_coeff.h     \
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
	${I_PATH}/conj_grad.h      \
	${I_PATH}/conj_grad_common.h \
	${I_PATH}/sg_params.h      \
	${I_PATH}/op_statistics.h  \
	${I_PATH}/mc_statistics.h  \
	${I_PATH}/monte_carlo.h    \
	${I_PATH}/mc_plot_slabs.h  \
	${I_PATH}/sensitivity.h  
#
# -------------------------------------------------------------------------
# -- LIBRARY OBJECTS
# -------------------------------------------------------------------------
#
LIBOBJ = \
	$(L)(ran1.o)\
	$(L)(func.o)\
	$(L)(dfunc.o)\
	$(L)(f1dim.o)\
	$(L)(df1dim.o)\
	$(L)(frprmn.o)\
	$(L)(linmin.o)\
	$(L)(dbrent.o)\
	$(L)(mnbrak.o)\
	$(L)(cost_function.o)\
	$(L)(calc_hessian.o)\
	$(L)(invert_hessian.o)\
	$(L)(svdcmp.o)\
	$(L)(pythag.o)\
	$(L)(gasdev.o)
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
	rm -f *.o *% core *.for
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
func.for : func.F  $(H_FILES)
func.o : func.for
#
dfunc.for : dfunc.F  $(H_FILES)
dfunc.o : dfunc.for
#
f1dim.for : f1dim.F
f1dim.o : f1dim.for
#
df1dim.for : df1dim.F
df1dim.o : df1dim.for
#
dbrent.for : dbrent.F
dbrent.o : dbrent.for
#
mnbrak.for : mnbrak.F
mnbrak.o : mnbrak.for
#
linmin.for : linmin.F
linmin.o : linmin.for
#
frprmn.for : frprmn.F
frprmn.o : frprmn.for
#
ran1.for : ran1.F
ran1.o : ran1.for
#
gasdev.for : gasdev.F
gasdev.o : gasdev.for
#
calc_hessian.for : calc_hessian.F $(H_FILES)
calc_hessian.o : calc_hessian.for
#
svdcmp.for : svdcmp.F
svdcmp.o : svdcmp.for
#
pythag.for : pythag.F
pythag.o : pythag.for
#
invert_hessian.for : invert_hessian.F 
invert_hessian.o : invert_hessian.for
#
cost_function.for : cost_function.F $(H_FILES)
cost_function.o : cost_function.for
#
$(L) : $(O_FILES)
