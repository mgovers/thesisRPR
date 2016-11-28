# Copyright (C) 2006 Dave Ireland <d.ireland@physics.gla.ac.uk>
#               2010 Pieter Vancraeyveld <Pieter.Vancraeyveld@UGent.be>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

MAINTAINERCLEANFILES 		= Makefile.in

EXTRA_DIST	=	makerootfile.cpp \
			test_strangecalc-wrapper.root \
			data/photo/iso.1/diff.clas.exfwang.05 \
			data/photo/iso.2/diff.clas.exfwang.05 \
			data/photo/iso.3/diff.fwangles.saphir.05 \
			data/photo/iso.1/diff.slac.t.69 \
			data/photo/iso.1/diff.spring8.04 \
			data/photo/iso.2/diff.spring8.04 \
			data/photo/iso.2/pho.fwang.graal.05 \
			data/photo/iso.1/pho.slac.t.79 \
			data/photo/iso.1/pho.spring8.03 \
			data/photo/iso.2/pho.spring8.03 \
			data/photo/iso.1/rec.clas-cmu.exfwang.03 \
			data/photo/iso.2/rec.clas-cmu.exfwang.03 \
			data/photo/iso.1/rec.desy.t.72 \
			data/photo/iso.2/rec.fwang.graal.05 \
			data/photo/iso.3/rec.fwangles.saphir.05 \
			numinput/mass.in \
			numinput/coupl.iso.1 \
			numinput/coupl.iso.2 \
			numinput/varinfo.1 \
			numinput/varinfo.2 \
			fit_specification \
			fit_specification2 \
			ok_strangecalc_br_test \
			ok_strangecalc_chisquared_test \
			ok_testTStrangecalc_test \
			ok_testTStrangecalc2_test \
			TStrangeModel_test_ok \
			Lagrangian_test_ok

INCLUDES		= $(all_includes) \
			-I@ROOTINCDIR@ \
			-I$(top_srcdir)/src


LDFLAGS			= $(all_libraries) \
			-L@ROOTLIBDIR@ \
			-L$(top_builddir)/src/.libs

LDADD			= $(top_builddir)/src/.libs/libStrangecalc.a \
			@ROOTLIBS@ @ROOTAUXLIBS@ @GSL_LIBS@

check_PROGRAMS = 	Lagrangian_test \
			TStrangeModel_test \
			CGLN_test \
			GaugeInvariance_test \
			Speed_test \
			CreateTStrangeModel \
			CCBoundaries \
			pionHighE_test \
			Observable_test \
			testTStrangeCalc
TStrangeModel_test_SOURCES	= TStrangeModel_test.cpp
Lagrangian_test_SOURCES 	= Lagrangian_test.cpp
CGLN_test_SOURCES 		= CGLN_test.cpp
GaugeInvariance_test_SOURCES 	= GaugeInvariance_test.cpp
CreateTStrangeModel_SOURCES 	= CreateTStrangeModel.cpp
Speed_test_SOURCES 		= Speed_test.cpp
CCBoundaries_SOURCES		= CCBoundaries.cpp
pionHighE_test_SOURCES 		= pionHighE_test.cpp
Observable_test_SOURCES 	= Observable_test.cpp
testTStrangeCalc_SOURCES =	testTStrangeCalc.cc

# Clean
#######
CLEANFILES              = *~ core *.flc

#       End
# ======================================================================

