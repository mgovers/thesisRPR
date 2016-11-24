# Copyright (C) 2012 Pieter Vancraeyveld <pieter.vancraeyveld@ugent.be>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

if USE_INTEL_COMPILER
AM_CXXFLAGS 	= 	-wd584 -wd161
endif

EXTRA_DIST	=	version.h.in \
			strangecalc_path.h.in

INCLUDES 	= 	$(all_includes) \
			-I@ROOTINCDIR@ \
			@GSL_CFLAGS@

pkginclude_HEADERS 	=	Structures.h \
				FormFactor.h \
				FourVector.h \
				FourVectorHandle.h \
				GammaStructure.h \
				Lagrangian.h \
				Matrix.h \
				TMatrixElement.h \
				TCachemap.h \
				TCalculateCGLNMatrix.h \
				TCalculateCGLNCoeff.h \
				TCalculateConsistentCoeff.h \
				TCalculateDiagram.h \
				TCurrent.h \
				Tensor.h \
				TensorRank2.h \
				TensorRank3.h \
				TensorRank4.h \
				TMatrixElementOld.h \
				TMatrixElementCGLN.h \
				TCGLNDiagnostic.h \
				BetaIncomplete.h \
				numtoa.h \
				fitting.h \
                                DataHandler.h \
				io_specific.h \
				strange_func.h \
				simanneal.h \
				version.h \
				calcMatrixElement.h \
				FormFactorParametrization.h \
				FormFactorSpecification.h \
				TStrangeCalc.h \
				TKinematics.h \
				TStrangeModel.h \
				TCalcInfo.h \
				TDataset.h \
				TMultiModel.h \
				TGridded.h \
				TAveragedObservable.h \
				strangecalc_path.h

noinst_HEADERS		=	TKinematicsDict.h \
				TStrangeModelDict.h \
				TCalcInfoDict.h \
				TDatasetDict.h \
				TMultiModelDict.h \
				TGriddedDict.h \
				TAveragedObservableDict.h

lib_LTLIBRARIES 	= 	libStrangecalc.la
libStrangecalc_la_SOURCES = 	FormFactor.cpp \
				FourVectorHandle.cpp \
				GammaStructure.cpp \
				Lagrangian.cpp \
				Matrix.cpp \
				TMatrixElement.cpp \
				TMatrixElementOld.cpp \
				TMatrixElementCGLN.cpp \
				TCalculateCGLNCoeff.cpp \
				TCalculateCGLNMatrix.cpp \
				TCalculateConsistentCoeff.cpp \
				TCurrent.cpp \
				TCalculateDiagram.cpp \
				TensorRank2.cpp \
				TensorRank3.cpp \
				TensorRank4.cpp \
				TCGLNDiagnostic.cpp \
				BetaIncomplete.cpp \
				numtoa.cpp \
				fitting.cpp \
                                DataHandler.cpp \
				io_specific.cpp \
				strange_func.cpp \
				simanneal.cpp \
				version.h \
				calcMatrixElement.cpp \
				FormFactorParametrization.cpp \
				FormFactorSpecification.cpp \
				TStrangeCalc.cc \
				TKinematics.cpp \
				TKinematicsDict.cpp \
				TStrangeModel.cpp \
				TStrangeModelDict.cpp \
				TCalcInfo.cpp \
				TCalcInfoDict.cpp \
				TDataset.cpp \
				TDatasetDict.cpp \
				TMultiModel.cpp \
				TMultiModelDict.cpp \
				TGridded.cpp \
				TGriddedDict.cpp \
				TAveragedObservable.cpp \
				TAveragedObservableDict.cpp

libStrangecalc_la_LDFLAGS =	$(all_libraries) \
				-L@ROOTLIBDIR@ \
				@GSL_LIBS@ \
				-release @PACKAGE_VERSION@ 
libStrangecalc_la_LIBADD =  	@ROOTLIBS@ @ROOTAUXLIBS@ \
				-lMathMore

bin_PROGRAMS 		 = 	strangecalc

strangecalc_SOURCES = 	strangecalc.cpp
strangecalc_LDADD = 	libStrangecalc.la \
			@GSL_LIBS@ \
			-lm


# Generate Cint Dictionary
%Dict.cpp %Dict.h:	%.h %LinkDef.h
			@ROOTCINT@ -f $*Dict.cpp -c -p $(INCLUDES) $(AM_CPPFLAGS) $^


# Get rid of Emacs backups as well...
clean-local:
		rm -f *~ core *.flc *Dict.cpp *Dict.h $(BUILT_SOURCES)
