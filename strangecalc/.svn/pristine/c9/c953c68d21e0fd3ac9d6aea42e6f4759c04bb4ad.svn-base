# Copyright (C) 2006 Dave Ireland <d.ireland@physics.gla.ac.uk>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

MAINTAINERCLEANFILES    = Makefile.in

testsubdir              = testSubDir

TESTS_ENVIRONMENT       = top_builddir=$(top_builddir)

TESTS                   = 	strangecalc_chisquared.test \
				test_iso1.test \
				test_iso2+3.test \
				test_iso6.test \
				test_iso11.test \
				Lagrangian.test\
				Lagrangian_gic.test\
				Lagrangian_cgln.test \
				GaugeInvariance_consistent.test \
				GaugeInvariance_lorentz.test \
				GaugeInvariance_lorentz_gaussian.test \
				GaugeInvariance_n.test \
				GaugeInvariance_varcutoff1.test \
				GaugeInvariance_varcutoff2.test \
				GaugeInvariance_y.test \
				Observables.test \
				pionHighE.test
#				TStrangeModel.test

EXTRA_DIST              = defs $(TESTS)

#	Get rid of Emacs backups as well...
clean-local:
		rm -f *~ *.flc
		rm -rf $(testsubdir)
