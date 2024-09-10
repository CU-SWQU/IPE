noinst_LIBRARIES = libspharm.a

libspharm_a_SOURCES = hwm14.f90 spher_harm.for

libspharm_a_FCFLAGS = $(FC_MODINC) .  $(FC_MODINC) ../shared

EXTRA_DIST = hwm93.f90


.PHONY: clean-modules

clean-modules:
	-test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)

clean-local: clean-modules

# dependencies
