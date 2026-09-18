
# =====================================================================
# Makefile for the diagbox f2py wrapper.
# Lives at the same level as src_ia/ and bin_ia/, and consumes the
# static library that src_ia's own makefile produces - it does not
# recompile any legacy model source itself.
#
# Before this will work you need to fill in the three unknowns marked
# below - don't guess at them, go confirm each with `make -n all` in
# src_ia/ and by reading src_ia/makefile directly.
# =====================================================================
 
ROOT = /home/prafter/f2py-diagbox/src_ia

MODNAME  = diagbox
SRC      = ../src_ia/diagboxmod/run_simulation.for
INCLDIR  = ../src_ia/incl
 
# --- Fill in from src_ia/makefile -----------------------------------
# 1) LIBDIR/LIBNAME: confirm via `make -n all` in src_ia/ - watch for
#    the `ar` command and see exactly what path/filename it produces.
LIBDIR   = ../src_ia/diagboxmod
LIBNAME  = diagboxmod
 
# 2) FFLAGS: must match what src_ia/makefile actually passes when
#    compiling diagboxmod's objects (check for -fno-automatic,
#    -std=legacy, and whatever -D defines come from src_ia/ifdefs).
#    Mismatched flags between this file and the prebuilt library is
#    exactly the "compiles fine, behaves differently" trap discussed
#    earlier - don't leave this as a guess.
FFLAGS   = -C -g  -fbacktrace -O0 -finit-real=snan -finit-integer=-2147483648 -fno-automatic -std=legacy -I$(ROOT)/incl -fcheck=all -finit-real=nan -ffpe-trap=invalid,zero,overflow -Wall -Wextra -Wno-compare-reals -Wno-unused-label -Wno-unused-variable -Wno-unused-dummy-argument -Wno-unused-parameter -Wno-line-truncation -Wno-cpp
 
# 3) Confirm the meson-backend equivalent of -L/-l for f2py -c before
#    trusting this line - the distutils-era syntax below may not be
#    honored the same way; see the buildtools/meson.html doc.
LINKFLAGS = -L$(abspath $(LIBDIR)) -l$(LIBNAME)
# ----------------------------------------------------------------------
 
.PHONY: diagbox clean
 
# Explicit dependency edge: force src_ia's library to be current
# before attempting to link against it, rather than relying on
# remembering to build it by hand first.
diagbox: $(abspath $(LIBDIR))/lib$(LIBNAME).a
	python3 -m numpy.f2py -m $(MODNAME) -c $(abspath $(SRC)) \
		-I$(abspath $(INCLDIR)) \
		--f77flags="$(FFLAGS)" \
		$(LINKFLAGS)
 
$(abspath $(LIBDIR))/lib$(LIBNAME).a:
	$(MAKE) -C $$(abspath $(LIBDIR))
 
clean:
	rm -f $(MODNAME)*.so
	rm -rf build/
 
