/*!
 * \file numtoa.cpp
 * \ingroup wrapper
 *
 * NUMbers TO Anything (and vice versa)
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

#include "numtoa.h"

/*!
 * Generates a unique key from an array of doubles of length 'elements'. 
 * Dimension of key should be at least 2*sizeof(double)*elements + 1 (for terminating char \0)
 */
void doublestochar(int elements, void * dArray, char* key)
{
  int keylen=sizeof(double)*elements;
  char* p = (char*)dArray;
  int nulchars[keylen];
  int nuls = 0;
  int i;
  for(i=0; i < keylen; ++i)
  {
    if (*p == '\0')
    {
      // replace '\0' (=end-of-string) by something more sensible.
      key[i]='_';
      nulchars[nuls]=i+1; // keep track of the '\0' location.
      ++nuls;
    }
    else{
      key[i] = *p;
    }
    ++p;
  }
  p=0;
  for(i=0; i < nuls; ++i)
  {
    // suffix the key with the locations of the replacements, converted to char.
    key[keylen + i] = ((char) nulchars[i]=='\0')?'A':(char) nulchars[i];
  }
  // Finish with an end-of-string character.
  key[keylen + nuls] = '\0';
}

/*!
 * Generates a unique key from an array of floats of length 'elements'. 
 * Dimension of key should be at least 2*sizeof(double)*elements + 1 (for terminating char \0)
 */
void floatstochar(int elements, void * fArray, char* key) 
// makes more sense if floating point errors cause small fluctuations in cc values
{
  int keylen = sizeof(float)*elements;
  char* p = (char*)fArray;
  int nulchars[keylen];
  int nuls = 0;
  int i;
  for(i=0; i < keylen; ++i)
  {
    if (*p == '\0')
    {
      // replace '\0' (=end-of-string) by something more sensible.
      key[i]='_';
      nulchars[nuls]=i+1; // keep track of the '\0' location.
      ++nuls;
    }
    else{
      key[i] = *p;
    }
    ++p;
  }
  p=0;
  for(i=0; i < nuls; ++i)
  {
    // suffix the key with the locations of the replacements, converted to char.
    key[keylen + i] = ((char) nulchars[i]=='\0')?'A':(char) nulchars[i];
  }
  // Finish with an end-of-string character.
  key[keylen + nuls] = '\0';
}
