# =====================================================================
# Makefile for the diagbox f2py wrapper.
# Lives alongside src_ia/ and bin_ia/. Links against whichever prebuilt
# static libraries run_simulation's actual call tree requires - it does
# not recompile legacy source itself.
#
# IMPORTANT: LIBS below currently lists diagboxmod and chem, which you
# have *confirmed* are required (chem via tb_uppst, called from
# input.F). Do NOT add numrecip or plot without first grep-ing the full
# call tree the way you did for chem - the earlier trace suggested
# run_simulation's forward-simulation path doesn't reach either of
# them (numrecip is Monte Carlo/optimization-only; plot is NCAR
# graphics output). Adding libraries you don't need still links (may
# even still work) but silently widens what your wrapper depends on
# for no reason - confirm need with grep before adding a directory
# here, same as you did for chem.
# =====================================================================
 
# Pin the interpreter explicitly rather than relying on whatever
# python3 happens to resolve to in the current shell - this is the
# ABI-mismatch lesson from earlier, made permanent instead of
# depending on remembering to `conda activate` the right env first.
PYTHON   = /usr/bin/python3
 
MODNAME  = diagbox
ROOT     = $(abspath ..)
SRC      = $(ROOT)/src_ia/diagboxmod/run_simulation.for
INCLDIR  = $(abspath $(ROOT)/src_ia/incl)
 
# --- Libraries this wrapper actually links against -------------------
# Each entry here is a (dir, name) pair confirmed via `make -n all`
# in that subdirectory and a grep-confirmed symbol dependency - not
# an assumption. Add entries the same way: verify, then add.
LIBDIRS  = $(abspath $(ROOT)/src_ia/diagboxmod) $(abspath $(ROOT)/src_ia/chem) $(abspath $(ROOT)/src_ia/numrecip)
LIBNAMES = diagboxmod chem numrecip
 
LFLAGS   = $(foreach d,$(LIBDIRS),-L$(d)) $(foreach n,$(LIBNAMES),-l$(n))
 
# Each library's archive spelled out explicitly as its own
# prerequisite, rather than derived positionally from LIBDIRS/LIBNAMES
# above - positional pairing across two parallel lists gets fragile
# and hard to read fast; explicit is worth the extra line here.
LIBFILES = $(ROOT)/src_ia/diagboxmod/libdiagboxmod.a \
           $(ROOT)/src_ia/chem/libchem.a \
		   $(ROOT)/src_ia/numrecip/libnumrecip.a 

 
# Flags must match what each library was actually compiled with -
# confirm against src_ia/makefile and each subdirectory's own
# makefile rather than trusting this copy. Mismatched flags between
# this compile and the prebuilt libraries is a silent-bug risk, not
# just a compile-failure risk.
FFLAGS   = -cpp -fno-automatic -std=legacy -I$(INCLDIR)
 
.PHONY: diagbox clean
 
diagbox: $(LIBFILES)
	$(PYTHON) -m numpy.f2py -m $(MODNAME) -c $(SRC) \
		-I$(INCLDIR) \
		--f77flags="$(FFLAGS)" \
		$(LFLAGS)
 
# One explicit rule per library, rather than a clever generic pattern
# rule - each subdirectory's makefile is a separate, independently
# authored build unit (per Nicolas Gruber's own makefile comments),
# so treating them uniformly would hide real differences between them
# (e.g. if flags or DIRS/SRC_DIRS differ between chem's and
# diagboxmod's makefiles).
$(ROOT)/src_ia/diagboxmod/libdiagboxmod.a:
	$(MAKE) -C $(ROOT)/src_ia/diagboxmod
 
$(ROOT)/src_ia/chem/libchem.a:
	$(MAKE) -C $(ROOT)/src_ia/chem
 
clean:
	rm -f $(MODNAME)*.so
	rm -rf build/
 
