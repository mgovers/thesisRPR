#
# ======================================================================
#
#   File:   doc/Makefile.am
#           ~~~~~~~~~~~~~~~~~~~~~~~
#
#   Description: Makefile to generate doxygen documentation.
#                The contents was copied from the config/ax_prog_doxygen.m4
#                documentation.
#
#   Author:    Pieter Vancraeyveld <pieter.vancraeyveld@ugent.be>
#
#   Copyright (c) 2012 Pieter Vancraeyveld <pieter.vancraeyveld@ugent.be>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.
#
# ======================================================================

if DX_COND_doc

doxyfile.stamp: $(top_builddir)/src/.libs/libStrangecalc.a
	$(DX_DOXYGEN) Doxyfile
	echo Timestamp > doxyfile.stamp

html-local: doxyfile.stamp

clean-local:
	rm -rf doxyfile.stamp $(top_builddir)/htmldoc/*

endif DX_COND_doc

#	End
# ======================================================================
