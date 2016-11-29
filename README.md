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
