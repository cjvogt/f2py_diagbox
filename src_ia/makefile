# -------------------------------------------------------------------------
# FILE    : makefile
# AUTHOR  : Nicolas Gruber
# LAST REV: 05-AUG-96  adapted for program diag_boxmodel
#           17-SEP-99  adapted for SUN and interannual variability version
#           28-AUG-01  adapted for LINUX
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
# -- File containing list of CPP define directives
#    (default: $(ROOT)/ifdefs)
# ----------------------------------------------------------------------
#

include $(ROOT)/ifdefs

# -------------------------------------------------------------------------
# -- COMPILATION MACROS 
# -------------------------------------------------------------------------

F_C      = gfortran
NCARF_C = /home/prafter/.conda/envs/ncl_stable/bin/ncargf90
!F_C = lf95
!NCARF_C = ncarglf95
#
#F_FLAGS  = 
F_FLAGS  = -C -g  -fbacktrace -O0 -finit-real=snan -finit-integer=-2147483648 -fno-automatic -std=legacy -I$(ROOT)/incl -fcheck=all -finit-real=nan -ffpe-trap=invalid,zero,overflow -Wall -Wextra -Wno-compare-reals -Wno-unused-label -Wno-unused-variable -Wno-unused-dummy-argument -Wno-unused-parameter -Wno-line-truncation -Wno-cpp
# -------------------------------------------------------------------------
# -- EXTERNAL OBJECTS FILES:
# -- example: ( EXT_O_FILES = /users/user_name/proj/lib/auxmat.o )  
# -------------------------------------------------------------------------

EXT_O_FILES = 

# -------------------------------------------------------------------------
# -- EXTERNAL LIBRARIES:
# -- example: ( EXT_LIB = /users/user_name/proj/lib/libmath.a )  
# -------------------------------------------------------------------------

EXT_LIB = $(ROOT)/diagboxmod/libdiagboxmod.a  $(ROOT)/plot/libplot.a \
          $(ROOT)/numrecip/libnumrecip.a $(ROOT)/chem/libchem.a

# -------------------------------------------------------------------------
# -- INTERNAL LIBRARIES:
# -- example: ( INT_LIB = LIBRARY2 ) environment variable
# -------------------------------------------------------------------------

INT_LIB = 

# -------------------------------------------------------------------------
# -- LOCATION OF INCLUDE FILES: (for directives like 'include filename')  
# -- example: ( I_PATH = /users/user_name/proj/include )  
# -------------------------------------------------------------------------

I_PATH = \
	 $(ROOT)/incl

O_FILES = 
#
I_FILES = 
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
I_FILES_DEBUG = \
	diagboxmod/run_simulation.for    diagboxmod/gasexchange.for        \
	diagboxmod/vert_diff.for         diagboxmod/intvertgrad_ent.for    \
	diagboxmod/masscons_ent.for      diagboxmod/episevent_ent.for      \
	diagboxmod/advection.for         diagboxmod/net_comm_prod.for      \
	diagboxmod/input.for             diagboxmod/advection1.for         \
	diagboxmod/advection2.for        diagboxmod/matear_ent.for
# 
X_FILES = diag_boxmodel_ia.x
#
# -------------------------------------------------------------------------
# -- DEPENDENCIES RULES 
# -------------------------------------------------------------------------

all:	
	@echo ""
	@echo " *****************************"
	@echo " *** UPDATING $(X_FILES)"
	@echo " *****************************"
	@echo "" 
	@echo " *** making all in chem ***" 
	cd chem ; make all 
	@echo ""
	@echo " *** making lib in plot ***"
	cd plot ; make all 
	@echo ""
	@echo " *** making lib in numrecip ***"
	cd numrecip ; make all
	@echo ""
	@echo " *** making lib in diagboxmod ***"
	cd diagboxmod ; make all
	@echo ""
	@echo " *** making $(X_FILES)***"
	make $(X_FILES) 
	mv $(X_FILES) ../bin_ia
	@echo ""
	@echo " *************"
	@echo " *** DONE. ***"
	@echo " *************"
	@echo ""

$(X_FILES): $(O_FILES)


clean:
	@echo ""
	@echo " *** CLEANING UP... ***"            
	-rm -f $(O_FILES) $(X_FILES) $(I_FILES) *.so
	@echo " *** DONE ***"

clean_all:
	@echo ""
	@echo " *** CLEANING UP ALL SUBDIRS ***"            
	@echo " *** clean all files in chem ***" 
	cd chem ; make clean 
	@echo ""
	@echo " *** clean all files in plot ***"
	cd plot ; make clean 
	@echo ""
	@echo " *** clean all files in numrecip ***"
	cd numrecip ; make clean
	@echo ""
	@echo " *** clean all files in diagboxmod ***"
	cd diagboxmod ; make clean

	-rm -f $(O_FILES) $(X_FILES) $(I_FILES)
	@echo " *** DONE ***"
            
clean_debug:
	@echo ""
	@echo " *** CLEANING UP... ***"            
	-rm -f $(I_FILES_DEBUG)
	@echo " *** DONE ***"            
 

# -------------------------------------------------------------------------
# -- Default Rules and Suffixes
# -------------------------------------------------------------------------

.SUFFIXES:
.SUFFIXES: .x .o .for .F

.for.x:
	@echo ""
	@echo " *** COMPILING AND LINKING $<... ***"
	$(NCARF_C) $(F_FLAGS) $< -o $@ \
               $(O_FILES) $(EXT_O_FILES) \
               $(EXT_LIB) $(INT_LIB)
	-rm -f $*.o $*.for

.F.for: 
	@echo ""
	@echo " *** PREPROCESSING $<... ***"
	$(F_C) $(D_FLAGS) $(F_FLAGS) -I$(I_PATH) -E $< > $*.for 
#	$(F_C) $(D_FLAGS) $(F_FLAGS) -I$(I_PATH) -P $< -o $*.for 

.for.o:   
	@echo ""
	@echo " *** COMPILING $<... ***"
#	$(F_C) $(D_FLAGS) $(F_FLAGS) -I$(I_PATH) -E $< > $*.for 
#	$(F_C) $(D_FLAGS) $(F_FLAGS) -I$(I_PATH) -P $< -o $*.for 
#	$(F_C) -c $(F_FLAGS) $*.for -o $*.o

# -------------------------------------------------------------------------
# -- Dependencies (explicit) 
#          if you have include files in a subdirectory they must be 
#          indicated together with the path
# -------------------------------------------------------------------------
#
#
diag_boxmodel_ia.for : diag_boxmodel_ia.F $(H_FILES)  ${I_PATH}/init_params.h
diag_boxmodel_ia.o : diag_boxmodel_ia.for

diagbox: run_simulation.f
	python3 -m numpy.f2py -m diagbox -c run_simulation.f