# thesisRPR

29/11/2016: 
   Structures.h:
       -in struct Observable: added short mint, default initialized as 0 => old model

29/11/2016:
   TCalculateConsistentCoeff.cpp:
       -in CalcA1(...), CalcA2(...), CalcA5(...): -created choice between with/without Mint (with Mint not yet
	                                           implemented!!)

29/11/2016:
   Structures.h:
       -added struct MintManager: contains additional information about the interaction vertex that are not
	passed to CalcA1(...) ... CalcA6(...) in TCalculateConsistentCoeff.cpp
	eg.: information about the s-channel hadronic FormFac, particletwo (instance of struct Properties,
	contains the information about the additional channel exchange), and information about possible
	Reggeization procedures of this particle.
       -in struct Observable: added field mintmanager instance of MintManager

29/11/2016:
   TCalculateConsistentCoeff.cpp:
       -in CalcA1(...) and CalcA2(...): -implemented version with Mint. Debug: inconsistency with old
		 			model in CalcA1(...): in old model, no minus sign appeared for
	 				s-channel Born exchange, but in own calculations one did appear.
	NOTE: CalcA5(...) NOT YET IMPLEMENTED!!!!! => not yet tested
	ALSO: no implementation of observ.mintmanager (part of args of CalcAi(...) ) was made yet
	 	=> must be fixed before testing...

30/11/2016:
   TCalculateConsistentCoeff.cpp:
	-in CalcA1(...) and CalcA2(...): made some major adjustments to Mint implem. concerning EM charge and FF
	-in CalcA5(...): implemented Mint and included the adjustments like in the above line.
	NOTE: no implementation of observ.mintmanager (part of arg of CalcAi(...) ) was made yet
	 	=> must be fixed before testing...
	-implemented default (observ.mintmanager not implemented => particletwo==NULL), not sure if work yet
	 	=> TEST!!!!

30/11/2016:
   TCalculateConsistentCoeff.cpp: Fixed some very important sign errors.
	NOTE: no implementation of observ.mintmanager (part of arg of CalcAi(...) ) was made yet
	 	=> must be fixed before testing...
	-implemented default (observ.mintmanager not implemented => particletwo==NULL), not sure if work yet
	 	=> TEST!!!!

