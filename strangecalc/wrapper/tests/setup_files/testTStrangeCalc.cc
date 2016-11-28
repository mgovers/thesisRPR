#include "TStrangeCalc.h"
using std::cerr;
using std::cout;
using std::endl;
int main(int argc, char *argv[])
{
  /* This program requires one argument:
   * The pad of a "fit_specification" file */
  if(argc != 2)
    {
      cerr << "Error: give location 'fit_specification' as argument!" << endl;
      exit(1);
    }
    
    //	Open logfile...
    FILE *log;
    if ( (log = fopen("tmplog", "w")) == NULL ) {
	cerr << "Error: cannot open " << "log" << "\n";
	exit (1);
    }
    
    double vertex[MAXNDIM] = {0};
    double *array_ptr;
    int isospin = 1; // will have to fix this

    TStrangeCalc c1(log,argv[1]);
    TStrangeCalc c2 = c1;   // test assignment operator
    TStrangeCalc c(c2);     // test copy constructor
    
    cout << "\nIntial GetChiSquared():\t" << c.GetChiSquared() << endl;
    c.PrintParticles();
    

    array_ptr = c.GetVertex();
    for ( int i=0 ; i<MAXNDIM ; i++ ) {
	vertex[i] = *(array_ptr+i);
    }

    for(int i=0; i<3; i++) {
	switch (i) {
	case 0:	
	    if(isospin == 1)
		vertex[0] = -3.11355;	
	    else	
		vertex[0] = 0.9;
	    break;
	case 1:
	    if(isospin == 1)
		vertex[0] = -3.0;
	    else
		vertex[0] = 1.0;
	    break;
	case 2:
	    if(isospin == 1)
		vertex[0] = -4.0;
	    else
		vertex[0] = 1.3;
	    break;
	default:;
	}
      
	c.SetVertex(vertex);
	c.SetChiSquared();
	cout << "\nGetChiSquared():\t" << c.GetChiSquared() << endl;
	c.PrintParticles();

    }
 
    cout << "\nGetDim():\t" << c.GetDim() << endl;
    cout.flush();
    fclose(log);
    

    return(0);
}
