# thesisRPR

24/12/2016:
   Added all changes that allow printing all inputted experiments to a single file with default format. Doing so
   is done by changing short printkinematics to 1 in struct Observable in Structures.h. NOTE: It is not checked
   whether the file already contains all data. It therefore is important to put observable.printkinematics to 0 
   whenever one is not sure that the data be written out to the file.

   Structures.h:
	-In struct Observable:
	 	-added short printkinematics
   Fitting.h:
	-In struct Photo and Elec:
	 	-added char[] experimentname
   Fitting.cpp:
	-In photo_chi(...) and electro_chi(...):
	 	-added some preparations and a call to printkinematics(...) in the if (observ->printkinematics)
	 	 statements
	-added printkinematics(...):
	 	-hardcoded!!! :::DEBUG:::FIX OR REMOVE BEFORE PUBLISHING?:::
   DataHandler.cpp:
	-In readData(...):
	 	-added a reading of the experiment name and saving it to the char[] experimentname in Photo and
	 	 Elec structs in Fitting.h
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
