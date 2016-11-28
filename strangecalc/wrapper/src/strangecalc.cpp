/*!
 * \file strangecalc.cpp
 * 
 * The 'strangecalc' program allows to do chi-squared calculations
 * and fitting processes for
 * photo- or electroproduction process, or a kaon capture reaction.
 * The input for chi-squared calculation 
 * and fitting processes is taken from the directories 
 * strange/calc/input/numinput and strange/calc/data.
 * The output is written in the directories strange/calc/output/##/iso.# 
 * 
 * A theoretical description of the processes and the model used in
 * this calculations can be found in the papers and thesisses found
 * at http://inwpent5.ugent.be/Publication.
 *
 * \remark The use of this program is limited and most of its functionality 
 * has been replaced by other parts of strangecalc.
 *
 * \author Tamara Corthals
 * \author Lesley De Cruz <lesley.decruz@ugent.be>
 * \author Stijn Janssen
 * \author Pieter Vancraeyveld <pieter.vancraeyveld@ugent.be>
 * \author Tom Vrancx <tom.vrancx@ugent.be>
 
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

#include <cstdio>
#include "strange_func.h"

int main(void) 
{  

  Observable observ = {0}; 

  FILE *ifp = stdin;
  FILE *ofp1 = stdout;
  FILE *ofp2 = stdout;
 
  observablespecification(&observ, ifp, ofp1, ofp2);

  if(observ.fitting || observ.chisquare) 
    fittingprocess(&observ); 
  else  
    error_exit("This functionality of strangecalc is depricated. Use the ROOT classes instead.\n");
 
  return 0;
 
}




