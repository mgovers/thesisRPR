/*!
 * \file TCachemap.h
 * \ingroup wrapper
 *
 * \author Pieter Vancraeyveld <pieter.vancraeyveld@ugent.be>
 
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

#ifndef __TCACHEMAP__
#define __TCACHEMAP__

#include <map>
#include <string>
#include <GammaStructure.h>
#include <FourVector.h>

//________________________________________________________________________
/*! \class TCachemap
 * \brief Implements a cached map for objects of type FourVector.
 *
 * The container size will never exceed TCachemap::kMaxCacheSize.
 * TCachemap is implemented as a wrapper of map<string, FourVector<T> >.
 *
 * We implement the container methods find and insert.
 *
 */

template< typename T >
class TCachemap
{
 public:
  //________________________________________________________________________
  /*! \brief Constructor
   * 
   * Initializes the container
   */
 TCachemap() : fCacheSize(0), fCachemap(std::map<std::string,T >()) {}

  //________________________________________________________________________
  /*! \brief Copy constructor
   * 
   * Copies the container
   */ 
 TCachemap(const TCachemap& rhs) : fCacheSize(rhs.fCacheSize), fCachemap(rhs.fCachemap) {}
  
  //________________________________________________________________________
  /*! \brief Assignment
   * 
   * Copies the container
   */
  TCachemap& operator=(const TCachemap& rhs) {
    if( this!=&rhs ) {
      fCacheSize = rhs.fCacheSize;
      fCachemap = rhs.fCachemap;
    }
    return *this;
  }

  //________________________________________________________________________
  /*! \brief Destructor
   */
  ~TCachemap() {}
  
  //________________________________________________________________________
  /*! \brief Searches the container for a key
   * 
   * \param key the key one wishes to find.
   * \param value will hold the associated value if the key is found.
   * \return true if the key is found.
   */
  bool find(const std::string& key, T& value) const {
    bool keyFound = false;
        
    typename std::map < std::string, T >::const_iterator it = fCachemap.find(key);
    if( it != fCachemap.end() ) {
      keyFound = true;
      value = it->second;
    }
  
    return keyFound;
  }

  //________________________________________________________________________
  /*! \brief Inserts a key-value pair in the container
   * 
   * This function makes sure the cache never contains more than
   * kMaxCacheSize elements.
   *
   * \param key the key one wishes to add.
   * \param value the value associated with the key.
   */
  void insert(const std::string& key, const T& value) {

    // Make sure the cacheMap doesn't grow out of control.
    if( fCacheSize==kMaxCacheSize ) {
      fCacheSize=0;
      fCachemap.clear();
    } else {
      ++fCacheSize;
    }
  
    // Insert the key-value pair in the cacheMap.
    fCachemap.insert(std::pair<std::string,T >(key,value));
  }

 private:
  //________________________________________________________________________
  // TCachemap data members
  static const int kMaxCacheSize = 9500; //!< clear container when cache size exceeds
  int fCacheSize; //!< Number of elements in cache
  std::map<std::string,T > fCachemap; //!< the wrapper container

}; // class TCachemap


#endif // __TCACHEMAP__
