// ***************************************************************************************
// strangecalc-test.cpp: Optimizer for Cut-off parameter
// ***************************************************************************************
#define OUTPUTWIDTH 102

#include "TStrangeModel.h"
#include "FormFactorParametrization.h"
#include "FormFactor.h"
#include "/home/mgovers/Software/strangecalcWorkDir/wrapper/share/models/NewModelMartijn/ECoupl_ResonWrapper.h"

using std::ostream;
using std::cout; 
using std::cerr; 
using std::endl;
using std::setw;
using std::setprecision;
using std::string;
using std::ofstream;
using std::atoi;
using std::atof;
using std::to_string;
using std::strcmp;
using std::strncmp;
using std::strcpy;
using std::replace;
using std::fixed;
using std::scientific;


int main(){
  string filename = "temp";
  ofstream output(filename);
  char nickname[] = "N1";
  int resVal = determineResonanceValue(nickname);
  FormFactor* ff = new FormFactor(monopole, 0., 0., 0.,
				    0., 0., 0., 0., 0., 0.,
				    0., 0., 0., 0., 0., 0.);
  FormFactor* tempff = specify_external_ff(resVal,1);
  if (tempff != NULL) ff = tempff;
  
  double current = 1.;

  double Q2 = 0.;
  while (Q2 < 5000000){
    current = ff->value(1.,Q2);
    output << Q2 << '\t' << current << '\t' << endl;
    Q2+=1000;
  }
  output.close();
  return 1; 
}

//g++ -I ~/Software/strangecalcWorkDir/wrapper/src -o ecouplFFWriter ecoupl_reson.o CppFileECouplResonWrapperTester.cpp $(root-config --cflags --libs) -L ~/Software/strangecalcWorkDir/wrapper/src/.libs -lStrangecalc-10.3 -I ~/Software/strangecalcWorkDir/wrapper/src

