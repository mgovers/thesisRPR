/*!
 * \file CGLN_test.cpp
 *
 * Very low-level test for the diagrams calculated in Lagrangian.cpp.
 *
 * \author Lesley De Cruz <lesley.decruz@ugent.be>
 
 * \copyright
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 */

#include <Lagrangian.h>
#include <Structures.h>
#include <TCurrent.h>
#include <TCalculateCGLNCoeff.h>
#include <iostream>
#include <strange_func.h>
#include <cstring>
#include <vector>
using namespace std;

int main(int argc, char *argv[])
{
    if (argc==1)
    {
        cout << "Please provide a model: old or new\n"<<endl;
        exit(1);
    }
    else if ((strcmp(argv[1],"new") && strcmp(argv[1],"old") ) )
    {
        cout << " Please enter \"old\" or \"new\"\n" << endl;
        exit(1);
    }

    Properties fakeparticle = {0};
    strcpy(fakeparticle.nickname,"N99");
    fakeparticle.mass=1700.;
    fakeparticle.width=100.;
    fakeparticle.E=2.0; /* electric charge divided by |e_electron|  */
    fakeparticle.G=3.0;  /* g_KYN, mult. by factor */
    fakeparticle.H=4.0; /* kappa, mult. by factor */
    fakeparticle.I=5.0;
    fakeparticle.J=6.0;
    fakeparticle.X=7.0; /* off-shell parameter */
    fakeparticle.Y=8.0; /* " */
    fakeparticle.Z=9.0; /* " */
    fakeparticle.r_kappa_n_p=2.0; /* ratio of trans. mag. moms k_{nN*}/k_{pN*} (1/2 N*)*/
    fakeparticle.r_kappa_1_n_p=3.0; /* k^(1)_{nN*}/k^(1)_{pN*} (3/2 N*) */
    fakeparticle.r_kappa_2_n_p=4.0; /* k^(2)_{nN*}/k^(2)_{pN*} (3/2 N*) */

    double mN=1000.;
    double mK=500.;
    double mY=1500.;

    double k = 1300.;
    double w = 1300.;
    double pk = 1180.12;
    double costheta_k=0.9;


    /*
    TKinematics tk("tk","tk",1,"wcm:costhkcm:qsquared",w,costheta_k,0);
    tk.Print();
    */
    FourVector<double> k4vect;
    FourVector<double> p4vect;
    FourVector<double> pK4vect;
    FourVector<double> pY4vect;

    k4vect = FourVector<double>(w,0,0,k);                // photon 4vector
    p4vect = FourVector<double>(sqrt(k*k+mN*mN),0,0,-k); // nucleon 4vector
    pK4vect = FourVector<double>(sqrt(pk*pk+mK*mK),
                                 pk*sqrt(1-costheta_k*costheta_k),
                                 0,pk*costheta_k);       // kaon 4vector
    pY4vect = FourVector<double>(sqrt(pk*pk+mY*mY),
                                 -1*pk*sqrt(1-costheta_k*costheta_k),
                                 0,-1*pk*costheta_k);    // hyperon 4vector

    // Calculate the mandelstam variables
    double s = ( p4vect+k4vect  ) * ( p4vect+k4vect );
    double t = ( k4vect-pK4vect ) * ( k4vect-pK4vect );
    double u = ( k4vect-pY4vect ) * ( k4vect-pY4vect );


    if (!strcmp(argv[1],"new"))
    {
        cout << "New model" << "\n" ;
        Observable obs = {0};
	obs.photoprod=1;
	Class emptyclass = {0};
	Class* particles = new Class[CLASSMAX];
	for (int i=0; i<CLASSMAX; i++)
	  particles[i] = emptyclass;
	
	particles[0].partic[0].mass = mN;
	particles[1].partic[0].mass = mK;
	particles[2].partic[0].mass = mY;
	particles[0].partic[0].E = 1.;// nucleon_charge
	particles[2].partic[0].E = 0.; // hyperon_charge
	
	TCalculateCGLNCoeff* calculator = TCalculateCGLNCoeff::GetCGLNCoeffCalculator();	
	TCurrent current(calculator,particles,&obs,s,t,u); // 4 components: epsilon (1,0,0,0),(0,1,0,0) etc
	
        std::vector < FourVector < complex< double > > > unitEpsilon;
        unitEpsilon.push_back(FourVector< complex <double> >(1.0,0.0,0.0,0.0));
        unitEpsilon.push_back(FourVector< complex <double> >(0.0,-1.0,0.0,0.0));
        unitEpsilon.push_back(FourVector< complex <double> >(0.0,0.0,-1.0,0.0));
        unitEpsilon.push_back(FourVector< complex <double> >(0.0,0.0,0.0,-1.0));

        for (int j=0;j<1000;j++)
	{
            for (int i = 0; i < CLASSMAX; i++)
            {
                // Test the new code
                if ( i==3 || (i >= 6 && i <=9 ))
                    fakeparticle.E=0.0;
                for (int mu = 0; mu < 4; mu++)
                {
                    current.Clear();
                    current.AddDiagram(i,fakeparticle);
                    GammaStructure currentGS = current.GetCurrentGS(unitEpsilon[mu],k4vect,p4vect,pY4vect);
                    if (j==100 || j==0)
                        cout << "Current["<<mu<<"] ="  << currentGS << "\n";
                }

                fakeparticle.E = 2.0;
            }
	}

	delete [] particles;
    }
    else
    {
        cout << "Old model" << "\n" ;
        for (int j=0;j<1000;j++)
            for (int i = 0; i < CLASSMAX; i++)
            {
                // Test the old code
                if ( i==3 || (i >= 6 && i <=9 ))
                    fakeparticle.E=0.0;
//                cout << "class " << i << "\n";
                FourVector<GammaStructure> diagramContribution = calculateDiagram[i](fakeparticle,
                        k4vect, p4vect, pK4vect, pY4vect);
                if (j==100 || j==0)
                    for (int mu = 0; mu < 4; mu++)
                        cout << "Current["<<mu<<"] ="  << diagramContribution[mu] << "\n";
            }

        fakeparticle.E = 2.0;
    }
    return 0;
}
