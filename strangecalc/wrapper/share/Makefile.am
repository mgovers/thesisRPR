#
# ======================================================================
#
#       File:   data/Makefile.am
#               ~~~~~~~~~~~~~~~~~~~~~~~
#
#       Description:    A subdirectory Makefile.am
#
#       Author: Pieter Vancraeyveld,
#               Pieter.Vancraeyveld@UGent.be   (C) 2010
#
#       Copying: See COPYING file that comes with the distribution
#
#
# ======================================================================

# The folder where the models will be installed
models_installation_dir	=	@datadir@/@PACKAGE@/models

# The folder where the data will be installed
datafiles_installation_dir	=	@datadir@/@PACKAGE@/data

# List of folders with the RPR models
model_dirs	=	models/rpr-2007 \
			models/rpr-2007-background \
			models/regge-2007 \
			models/rpr-2011 \
			models/rpr-2011-background \
			models/regge-2011 \
			models/VR-2014

# List of folders with the bare data files
data_dirs	=	data/photo \
			data/electro \
			data/radcap \
			data/deuteron

# Specify the ROOT TDatasets
datadir		=	$(datafiles_installation_dir)
data_DATA	=	data/data_iso1_HE.root \
			data/data_iso2_HE.root \
			data/data_iso3.root \
			data/data_iso5.root \
			data/data_iso8.root \
			data/electrodata_iso1.root \
			data/data_iso1.root \
			data/data_iso2.root \
			data/data_iso4.root \
			data/data_iso6.root \
			data/data_radcap.root \
			data/electrodata_iso2.root

# At first we install a dummy file in the models directory. This way
# we benefite from automake's platform-independent stuff.
# In the post-install hook, the dummy file is deleted and the actual
# models are copied.
modelsdir	=	$(models_installation_dir)
models_DATA	=	models/dummy

EXTRA_DIST	= 	$(data_DATA) $(data_dirs) $(model_dirs)

# ======================================================================

# We create a dummy file in the models directory. See supra why we do this.
models/dummy:;
	touch models/dummy

# Because the bare data and the models are held in folders, the default install
# rule would not work. We solve this with a post-install hook
install-data-hook:
	rm -f "$(DESTDIR)$(modelsdir)/../$(models_DATA)" || exit $$?; \
	list='$(model_dirs)'; test -n "$(modelsdir)" || list=; \
        for p in $$list; do \
	  rsync -art --delete --exclude ".svn" --exclude "fit_specification.in" "$$p" "$(DESTDIR)$(modelsdir)/" || exit $$?; \
	done; \
	list='$(data_dirs)'; test -n "$(datadir)" || list=; \
        for p in $$list; do \
	  rsync -art --delete --exclude ".svn" "$$p" "$(DESTDIR)$(datadir)/" || exit $$?; \
	done;

# Because the bare data and models are held in folders, the got installed through
# a hook. Hence we uninstall them in the same way
uninstall-hook:
	list='$(model_dirs)'; test -n "$(modelsdir)" || list=; \
        files=`for p in $$list; do echo $$p; done | sed -e 's|^.*/||'`; \
        test -n "$$files" || exit 0; \
        echo " ( cd '$(DESTDIR)$(modelsdir)' && rm -rf" $$files ")"; \
        cd "$(DESTDIR)$(modelsdir)" && rm -rf $$files \
	list='$(data_dirs)'; test -n "$(datadir)" || list=; \
        files=`for p in $$list; do echo $$p; done | sed -e 's|^.*/||'`; \
        test -n "$$files" || exit 0; \
        echo " ( cd '$(DESTDIR)$(datadir)' && rm -rf" $$files ")"; \
        cd "$(DESTDIR)$(datadir)" && rm -rf $$files

#	End
# ======================================================================
