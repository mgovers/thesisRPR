#
# ======================================================================
#
#       File:   strangeweb/Makefile.am
#               ~~~~~~~~~~~~~~~~~~~~~~~
#
#       Description:    A subdirectory Makefile.am
#
#       Author: Pieter Vancraeyveld,
#               Pieter.Vancraeyveld@UGent.be   (C) 2009
#
#       Copying: See COPYING file that comes with the distribution
#
#
# ======================================================================

INCLUDES		= $(all_includes) -I@ROOTINCDIR@ \
			-I$(top_srcdir)/src

MAINTAINERCLEANFILES =	Makefile.in

EXTRA_DIST =		strangeweb.cpp.in

# Compile and link strangeweb
##########################
bin_PROGRAMS =		strangeweb
strangeweb_SOURCES	= 	strangeweb.cpp
strangeweb_LDFLAGS	= 	$(all_libraries) \
				-L@ROOTLIBDIR@ \
				-L$(top_builddir)/src/.libs
strangeweb_LDADD =		@ROOTLIBS@ @ROOTAUXLIBS@ \
				-lMathMore \
				-lStrangecalc-@PACKAGE_VERSION@

# Clean
#######
CLEANFILES		= *~ core *.flc


#	End
# ======================================================================
