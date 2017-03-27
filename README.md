# thesisRPR

27/03/2017:
   Structures.h:
	-In mintmanager:
	 	added a field int use_different_tchannel_emffs: should an alternative EMFF for the kaon exchanges? 1 or 0
  	 	added a field FormFactor *tChannelFormFactor[5]: an array containing the alternative EMFFs
	 		for the kaon exchanges (max. 3 trajectories) 
	 	both tested that reduces to Old if:
	 	 	mintmanager.particletwo.FormFactorE = new FormFactor(monopole,1300.)
	 	 	mintmanager.tChannelFormFactor[0] = new FormFactor(monopole,1300.)
	 	 	mintmanager.tChannelFormFactor[1] = new FormFactor(monopole,1300.)
	-In TCurrent.cpp:
	 	added a statement that allows changing the t-channel EMFFs (cfr. adj. in Structures.h)

20/03/2017:
   Structures.h:
	-Added global field external_emff_skip_amound which contains the size of mintmanager.tempskiparray
	-In mintmanager:
	 	added field int tempskiparray[external_emff_skip] to allow plotting switching measured EMFF
	 		to Bonn (as in old model)
   TCurrent.cpp:
	-Implemented tempskiparray

20/03/2017:
   TStrangeModel.h/.cpp:
	-Removed UpdateEMFF() (last update did not work for some reason)
   TCurrent.cpp:
	-Implemented a (very ugly) way to force using the measured EMFF extracted from helicity amplitudes
	TESTED VERSION

17/03/2017:
   ecoupl_reson.f:
	-Added file to directory
   Makefile.am:
	-Added a way to build Fortran code of ecoupl_reson.f (also in wrapper/strangeweb/Makefile.am)
   FormFactorParametrization.h/.cpp:
	-Now includes all form factors obtained from helicity amplitudes and necessary code to include it (cfr. 
	 	ECoupl_ResonWrapper.h
	TESTED VERSION
   FormFactorSpecification.h/.cpp:
	-Now includes all specifications to include form factors obtained from helicity amplitudes (cfr. 
	 	ECoupl_ResonWrapper.h)
   	TESTED VERSION
   TStrangeModel.h/.cpp:
	-Now contains a (public) method UpdateEMFF() to update the EMFF after observ has been changed (which
	 	can only be done after the first initialization in which all FF are specified)

17/03/2017:
   wrapper/share/models/ :
	Added folder NewModelMartijn:
	 	-Contains all code concerning EMFF from measured helicity amplitudes
	-ecoupl_reson.f / .o:
	 	-File containing the Fortran code of https://userweb.jlab.org/~isupov/couplings/
	-ECoupl_ResonWrapper.h:
	 	-File containing the inclusion code for running the Fortran code of ecoupl_reson.f in the
	 	 	C++ code in safe mode with correct units and also for calculating the FF from the 
		 	measured helicity amplitudes. Used as real cpp file instead of header file due to
		 	inclusion problems during running of makefile (library assumes that .o file is included
		 	but it is not included because unknown how to implement this while not removing
		 	ecoupl_reson.o. I tried to fix this, but the straightforward way did not work.)
	 	-external::ecoupl_(int &i, int &j, float &Q2):
		 	Including fortran code (safe environment due to namespace): NOT EXPLICITLY CALLED
		 	EXCEPT IN THE WRAPPER 
	 	-ecoupling(int i,int j, double Q2):
		 	Wrapper for external::ecoupl_(...)
	 	-measured_*INSERT RESONANCE*_*INSERT FORMFACTOR TYPE*(double Q2):
		 	Calculator of EMFF by using the ecoupling wrapper
	 	-specify_external_ff(int resonanceValue, int formfactorType):
		 	Initializer of measured_*INSERT RESONANCE*_*INSERT FORMFACTOR TYPE*(...): returns
		 	currentFF to right EMFF (determined by resonanceValue and formfactorType)
	 	-determineResonanceValue(char* nickname):
		 	Method determining if the measured helicity amplitudes exist for certain resonance name
		 	and which parameter is passed in the Fortran code. The returned value is also used
		 	when initializing via specify_external_ff(...)

08/03/2017:
   MERGE: RPR2011MintImplementation and RPR2011WorkDir1 merged together
	@END OF PHOTOPRODUCTION

08/03/2017:
   TCalculateCGLNCoeff.cpp:
	-Just a commit (forgotten earlier)

02/03/2017:
   Structures.h:
	-In mintmanager:
	 	-Added a field that allows for the alternative Reggeization scheme (multiply all Born terms
	 		by (t-M_K^2)*ReggePropagator(t) )
   TCalculateConsistentCoeff.cpp:
	 	-Fixed bug concerning ReggePropagator2 when it went through 0
	 	-Implemented the alternative Reggeization scheme (tested)

20/02/2017:
   Structures.h:
	-In mintmanager:
	 	-Added fields short FsEqReggeFt, short FsEqFt and short reggeasinTchannel to make various
	 		degrees of difference from origional RPR-2011 model possible
	 	-Added field short onlyBackgroundContributions to calculate only the background contributions
   TCalculateConsistentCoeff.cpp:
	-In CalcA1(...), CalcA2(...), CalcA5(...):
	 	-Implemented FsEqReggeFt, FsEqFt and reggeasinTchannel (see Structures.h (above) )
	 	-Implemented calculating only background contributions (ONLY classindex 1 and classindex 4!!!)

15/02/2017:
   Structures.h:
	-In mintmanager:
	 	-Added a double field xi, implementing the xi freedom
   TCalculateConsistentCoeff.cpp:
	-In CalcA1(...), CalcA2(...), CalcA5(...):
	 	-Added a way to calculate the case that F_s = \mathcal{F}_t
	 	-Implemented the xi freedom
	NOTE:	Tested if the old RPR-2011 model is obtained from the new oneby putting F_s = \mathcal{F}_t and
	 	putting xi to xi=1/2, OK!

13/02/2017:
   TCalculateConsistentCoeff.cpp:
	-In CalcA1(...), CalcA2(...), ... , CalcA6(...):
	 	-Added a statement implementing onlyBornContributions because it didn't work properly
27/12/2016:
   Structures.h:
	-In struct MintManager:
	 	-Added a field short onlyBornContributions to allow showing only those contributions
   TCalculateCGLNCoeff.cpp:
	-In CalcA(...):
	 	-Added an if statement to allow calculating only the contributions by the Born terms

26/12/2016:
   Structures.h:
	-In struct Particle:
	 	-Added a DEBUG comment about particle.H and particle.G => interchanged in implementation in old
	 		model (TCalculateConsistentCoeff.cpp -> CalcA1(...),...,CalcA6(...)) => copied in new
		 	implementation.
   fitting.cpp:
	-In printkinematics(...):
	 	-Added a comment when calculating diffcs/dt: diffcs/dt is given in units mub/GeV^2 instead of
		 	mub/MeV^2 !!!

24/12/2016:
   Added forgotten files from last commit

   TStrangeModel.cpp:
	-In ChiSquared(...):
	 	-Deleted some debug statements

24/12/2016:
   Added all changes that allow printing all inputted experiments to a single file with default format. Doing so
   is done by changing short printkinematics to 1 in struct Observable in Structures.h. NOTE: It is not checked
   whether the file already contains all data. It therefore is important to put observable.printkinematics to 0 
   whenever one is not sure that the data be written out to the file.

   Structures.h:
	-In struct Observable:
	 	-added short printkinematics
   fitting.h:
	-In struct Photo and Elec:
	 	-added char[] experimentname
   fitting.cpp:
	-In photo_chi(...) and electro_chi(...):
	 	-added some preparations and a call to printkinematics(...) in the if (observ->printkinematics)
	 	 statements
	-added printkinematics(...):
	 	-hardcoded!!! :::DEBUG:::FIX OR REMOVE BEFORE PUBLISHING?:::
   DataHandler.cpp:
	-In readData(...):
	 	-added a reading of the experiment name and saving it to the char[] experimentname in Photo and
	 	 Elec structs in fitting.h
   TStrangeModel.cpp:
	-In ChiSquared(...):
	 	-added statements that allow calling printkinematics(...)

12/12/2016:
   strangecalc_path.h:
	-Changed default path to actual path (strangecalcWorkDir instead of strangecalc)

08/12/2016:
   Structures.h:
	-In struct Mintmanager:
	 	-added NoTchannelBornContribution, NoSchannelBornContribution, NoInteractionBornContribution
		 	to allow turning off those contributions

08/12/2016:
   strangecalc_path.h:
	 	-changed path from 'old' installation file to new one (old one created bugg in chisquared(..) 
	 	 	in TStrangeModel.h)

05/12/2016:
   Structures.h:
	-in MintManager:
	 	-added .a , .d , .f
   TCalculateConsistentCoeff.cpp:
	-in CalcA1(...): added a-term, d-term, f-term in New model
   TStrangeModel.h:
	 	-added *Observable GetObserv() to make it possible to easily switch between models
		 	=> :::DEBUG:::DELETED BEFORE PUBLISHING?:::

   	NOTE: -First real results suggest an enery dependend hadronic cut-off parameter for s-channel Born:
	      		-For low energies: low cut-off parameter
		 	-For high energies: large cut-off parameter
	       Can change by means of adding d and f => TEST!!!!
	      -First results suggest f=0. So far, no real predictions can be made about d, but probably negative
	       Both can still change by means of changing the cut-off parameter => TEST!!!!

01/12/2016:
   TCalculateConsistentCoeff.cpp:
	-Fixed a few minor errors
	-Fixed the new implementation so that OK for photoproduction if Mint put to 0 and Ms multiplied by 
	 	ReggePropagator of t-channel (equivalent to the old model RPR-2011) (Not yet tested for electro)
	NOTE: VERY MODEL DEPENDENT, even the differential cross section diffcs for photoproduction!!:
	 	Factor 2 too large, even for costhk = 1.0 (9.849379e-01 instead of 4.679263e-01),
	 	Factor 30 too large 	 for costhk = 0.0 (6.991458e-01 instead of 4.894534e-02),
	 	Up to 100 for backwards (but |t| small not valid anymore... but still a bad result)
	 	=> What to do with this? parameter a cannot fix this (photoprod). Maybe d,f,A,g,h can?

30/11/2016:
   TCalculateConsistentCoeff.cpp: Fixed some very important sign errors.
	NOTE: no implementation of observ.mintmanager (part of arg of CalcAi(...) ) was made yet
	 	=> must be fixed before testing...
	-implemented default (observ.mintmanager not implemented => particletwo==NULL), not sure if work yet
	 	=> TEST!!!!

30/11/2016:
   TCalculateConsistentCoeff.cpp:
	-in CalcA1(...) and CalcA2(...): made some major adjustments to Mint implem. concerning EM charge and FF
	-in CalcA5(...): implemented Mint and included the adjustments like in the above line.
	NOTE: no implementation of observ.mintmanager (part of arg of CalcAi(...) ) was made yet
	 	=> must be fixed before testing...
	-implemented default (observ.mintmanager not implemented => particletwo==NULL), not sure if work yet
	 	=> TEST!!!!

29/11/2016:
   TCalculateConsistentCoeff.cpp:
       -in CalcA1(...) and CalcA2(...): -implemented version with Mint. Debug: inconsistency with old
		 			model in CalcA1(...): in old model, no minus sign appeared for
	 				s-channel Born exchange, but in own calculations one did appear.
	NOTE: CalcA5(...) NOT YET IMPLEMENTED!!!!! => not yet tested
	ALSO: no implementation of observ.mintmanager (part of args of CalcAi(...) ) was made yet
	 	=> must be fixed before testing...

29/11/2016:
   Structures.h:
       -added struct MintManager: contains additional information about the interaction vertex that are not
	passed to CalcA1(...) ... CalcA6(...) in TCalculateConsistentCoeff.cpp
	eg.: information about the s-channel hadronic FormFac, particletwo (instance of struct Properties,
	contains the information about the additional channel exchange), and information about possible
	Reggeization procedures of this particle.
       -in struct Observable: added field mintmanager instance of MintManager

29/11/2016: 
   Structures.h:
       -in struct Observable: added short mint, default initialized as 0 => old model

29/11/2016:
   TCalculateConsistentCoeff.cpp:
       -in CalcA1(...), CalcA2(...), CalcA5(...): -created choice between with/without Mint (with Mint not yet
	                                           implemented!!)

#END OF FILE
