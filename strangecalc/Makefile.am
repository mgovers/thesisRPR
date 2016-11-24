#
# ======================================================================
#
#	File: 	Makefile.am
#		~~~~~~~~~~~
#
#	Description: 	Top-level Makefile.am. Process this file with
#			automake to produce Makefile.in  
#
#	Author:	Dave Ireland, d.ireland@physics.gla.ac.uk (C) 2006
#
#	Copying: See COPYING file that comes with the distribution
#
#
# ======================================================================

ACLOCAL_AMFLAGS         = -I config

# not a GNU package. You can remove this line, if
# have all needed files, that a GNU package needs
AUTOMAKE_OPTIONS 	= foreign 

SUBDIRS 		= src \
			  tests/setup_files tests doc strangeweb \
			  share


auxdir                  = @ac_aux_dir@
ACLOCAL                 = aclocal -I $(auxdir)


AUX_DIST                = $(ac_aux_dir)/config.guess \
                        $(ac_aux_dir)/config.sub \
                        $(ac_aux_dir)/install-sh \
                        $(ac_aux_dir)/ltmain.sh \
                        $(ac_aux_dir)/missing \
                        $(ac_aux_dir)/mkinstalldirs

EXTRA_DIST              = bootstrap stamp-h

MAINTAINERCLEANFILES    = Makefile.in aclocal.m4 configure in $(AUX_DIST)

pkginclude_HEADERS      = config.h

#       Get rid of Emacs backups as well...
clean-local:
		rm -f *~ *.flc 



#	End
# ======================================================================
