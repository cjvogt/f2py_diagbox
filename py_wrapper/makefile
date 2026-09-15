# Heres what to do, vaguely:
# Included all the archived libraries from src compilation to aid in the module construction

# Still working in src, just writing to this new folder
ROOT = /home/prafter/f2py-diagbox/src_ia

I_PATH = \
	 $(ROOT)/incl

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

EXT_LIB = $(ROOT)/diagboxmod/libdiagboxmod.a  $(ROOT)/plot/libplot.a \
          $(ROOT)/numrecip/libnumrecip.a $(ROOT)/chem/libchem.a

# Now do something with these...

diagbox: run_simulation.for
	python3 -m numpy.f2py -m diagbox -c run_simulation.for