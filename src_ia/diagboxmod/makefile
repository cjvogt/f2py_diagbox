# -------------------------------------------------------------------------
# FILE    : makefile
# AUTHOR  : Nicolas Gruber
# LAST REV: 05-AUG-96  adapted for program diag_boxmodel
#           17-SEP-99  adapted for SUN and interannual variability version
#                      special version for diagboxmod library part
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

LIBNAME = libdiagboxmod.a
L = $(LIBNAME)

# -------------------------------------------------------------------------
# -- LOCATION OF INCLUDE FILES: (for directives like 'include filename')  
# -- example: ( I_PATH = /users/user_name/proj/include )  
# -------------------------------------------------------------------------

I_PATH = \
	 $(ROOT)/incl

O_FILES = \
	intro.o\
	init_params.o\
	userinput.o\
	run_simulation.o\
	set_intval_to_zero.o\
	assign_params.o\
	set_sums_to_zero.o\
	init_val.o\
	input.o\
	gasexchange.o\
	vert_diff.o\
	intvertgrad_ent.o\
	masscons_ent.o\
	episevent_ent.o\
	matear_ent.o\
	advection1.o\
	advection2.o\
	net_comm_prod.o\
	calc_int_val.o\
	calc_run_stat.o\
	write_std_output.o\
	write_op_output.o\
	write_rates_tofile.o\
	write_fluxes_tofile.o\
	write_solution_tofile.o\
	write_mc_rates_tofile.o\
	write_mc_fluxes_tofile.o\
	write_obs_tofile.o\
	write_rateanomalies_tofile.o\
	write_fluxanomalies_tofile.o\
	write_mc_fluxanomalies_tofile.o\
	compute_anomalies.o\
	compute_runavg.o\
	compute_mc_runavg.o\
	harmonic.o\
	delta_to_ratio.o\
	c13_to_delta.o\
	delta_to_c13.o\
	rates_to_fluxes.o\
	fluxes_to_rates.o\
	getdate.o\
	keep_mc_statistics.o\
	calc_mc_statistics.o\
	write_mc_output.o\
	write_sn_output.o\
	spline.o\
	iplin.o
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
	$(L)(intro.o)\
	$(L)(init_params.o)\
	$(L)(userinput.o)\
	$(L)(run_simulation.o)\
	$(L)(set_intval_to_zero.o)\
	$(L)(assign_params.o)\
	$(L)(set_sums_to_zero.o)\
	$(L)(init_val.o)\
	$(L)(input.o)\
	$(L)(gasexchange.o)\
	$(L)(vert_diff.o)\
	$(L)(intvertgrad_ent.o)\
	$(L)(masscons_ent.o)\
	$(L)(episevent_ent.o)\
	$(L)(matear_ent.o)\
	$(L)(advection1.o)\
	$(L)(advection2.o)\
	$(L)(net_comm_prod.o)\
	$(L)(calc_int_val.o)\
	$(L)(calc_run_stat.o)\
	$(L)(write_std_output.o)\
	$(L)(write_op_output.o)\
	$(L)(write_rates_tofile.o)\
	$(L)(write_fluxes_tofile.o)\
	$(L)(write_solution_tofile.o)\
	$(L)(write_mc_rates_tofile.o)\
	$(L)(write_mc_fluxes_tofile.o)\
	$(L)(write_obs_tofile.o)\
	$(L)(write_rateanomalies_tofile.o)\
	$(L)(write_fluxanomalies_tofile.o)\
	$(L)(write_mc_fluxanomalies_tofile.o)\
	$(L)(compute_anomalies.o)\
	$(L)(compute_runavg.o)\
	$(L)(compute_mc_runavg.o)\
	$(L)(harmonic.o)\
	$(L)(delta_to_ratio.o)\
	$(L)(c13_to_delta.o)\
	$(L)(delta_to_c13.o)\
	$(L)(rates_to_fluxes.o)\
	$(L)(fluxes_to_rates.o)\
	$(L)(getdate.o)\
	$(L)(keep_mc_statistics.o)\
	$(L)(calc_mc_statistics.o)\
	$(L)(write_mc_output.o)\
	$(L)(write_sn_output.o)\
	$(L)(spline.o)\
	$(L)(iplin.o)

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
userinput.for : userinput.F $(H_FILES)
userinput.o : userinput.for
#
intro.for : intro.F
intro.o : intro.for
#
init_params.for : init_params.F
init_params.o : init_params.for
#
run_simulation.for : run_simulation.F  $(H_FILES)
run_simulation.o : run_simulation.for
#
assign_params.for : assign_params.F $(H_FILES)
assign_params.o : assign_params.for 
#
set_intval_to_zero.for : set_intval_to_zero.F $(H_FILES)
set_intval_to_zero.o : set_intval_to_zero.for
#
set_sums_to_zero.for : set_sums_to_zero.F $(H_FILES)
set_sums_to_zero.o : set_sums_to_zero.for
#
init_val.for : init_val.F $(H_FILES)
init_val.o : init_val.for
#
input.for : input.F $(H_FILES)
input.o : input.for
#
gasexchange.for : gasexchange.F $(H_FILES)
gasexchange.o : gasexchange.for
#
vert_diff.for : vert_diff.F $(H_FILES)
vert_diff.o : vert_diff.for
#
intvertgrad_ent.for : intvertgrad_ent.F $(H_FILES)
intvertgrad_ent.o : intvertgrad_ent.for
#
masscons_ent.for : masscons_ent.F $(H_FILES)
masscons_ent.o : masscons_ent.for
#
episevent_ent.for : episevent_ent.F $(H_FILES)
episevent_ent.o : episevent_ent.for
#
matear_ent.for : matear_ent.F $(H_FILES)
matear_ent.o : matear_ent.for 
#
advection1.for : advection1.F $(H_FILES)
advection1.o : advection1.for
#
advection2.for : advection2.F $(H_FILES)
advection2.o : advection2.for
#
net_comm_prod.for : net_comm_prod.F $(H_FILES)
net_comm_prod.o : net_comm_prod.for
#
calc_int_val.for : calc_int_val.F $(H_FILES)
calc_int_val.o : calc_int_val.for
#
calc_run_stat.for : calc_run_stat.F $(H_FILES)
calc_run_stat.o : calc_run_stat.for
#
calc_mc_statistics.for : calc_mc_statistics.F $(H_FILES)
calc_mc_statistics.o : calc_mc_statistics.for
#
keep_mc_statistics.for : keep_mc_statistics.F $(H_FILES)
keep_mc_statistics.o : keep_mc_statistics.for
#
write_std_output.for : write_std_output.F $(H_FILES)
write_std_output.o : write_std_output.for
#
write_op_output.for : write_op_output.F $(H_FILES)
write_op_output.o : write_op_output.for
#
write_mc_output.for : write_mc_output.F $(H_FILES)
write_mc_output.o : write_mc_output.for
#
write_sn_output.for : write_sn_output.F $(H_FILES)
write_sn_output.o : write_sn_output.for
#
write_rates_tofile.for : write_rates_tofile.F $(H_FILES)
write_rates_tofile.o : write_rates_tofile.for
#
write_solution_tofile.for : write_solution_tofile.F $(H_FILES)
write_solution_tofile.o : write_solution_tofile.for
#
write_fluxes_tofile.for : write_fluxes_tofile.F $(H_FILES)
write_fluxes_tofile.o : write_fluxes_tofile.for
#
write_mc_rates_tofile.for : write_mc_rates_tofile.F $(H_FILES)
write_mc_rates_tofile.o : write_mc_rates_tofile.for
#
write_mc_fluxes_tofile.for : write_mc_fluxes_tofile.F $(H_FILES)
write_mc_fluxes_tofile.o : write_mc_fluxes_tofile.for
#
write_obs_tofile.for : write_obs_tofile.F $(H_FILES)
write_obs_tofile.o : write_obs_tofile.for
#
write_rateanomalies_tofile.for : write_rateanomalies_tofile.F $(H_FILES)
write_rateanomalies_tofile.o : write_rateanomalies_tofile.for
#
write_fluxanomalies_tofile.for :  write_fluxanomalies_tofile.F $(H_FILES)
write_fluxanomalies_tofile.o :  write_fluxanomalies_tofile.for
#
write_mc_fluxanomalies_tofile.for :  write_mc_fluxanomalies_tofile.F $(H_FILES)
write_mc_fluxanomalies_tofile.o :  write_mc_fluxanomalies_tofile.for
#
compute_anomalies.for : compute_anomalies.F $(H_FILES)
compute_anomalies.o : compute_anomalies.for
#
compute_runavg.for : compute_runavg.F $(H_FILES)
compute_runavg.o : compute_runavg.for
#
compute_mc_runavg.for : compute_mc_runavg.F $(H_FILES)
compute_mc_runavg.o : compute_mc_runavg.for
#
harmonic.for : harmonic.F
harmonic.o : harmonic.for
#
rates_to_fluxes.for : rates_to_fluxes.F
rates_to_fluxes.o : rates_to_fluxes.for
#
fluxes_to_rates.for : fluxes_to_rates.F
fluxes_to_rates.o : fluxes_to_rates.for
#
delta_to_ratio.for : delta_to_ratio.F
delta_to_ratio.o : delta_to_ratio.for
#
c13_to_delta.for : c13_to_delta.F
c13_to_delta.o : c13_to_delta.for
#
delta_to_c13.for : delta_to_c13.F
delta_to_c13.o : delta_to_c13.for
#
spline.for : spline.F
spline.o : spline.for
#
iplin.for : iplin.F
iplin.o : iplin.for
#
$(L) : $(O_FILES)

