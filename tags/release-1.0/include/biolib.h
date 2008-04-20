/*********************************************************************************************************

 Library name     :   Biolib++ (Rutines for biological applications)
 File name        :   biolib.h
 Defined Objects  :   matrix  

 Programmer       :   Victor E. Bazterra 
 Directed by      :   Julio C. Facelli
                      
 Initialized      :   Wed Nov 12 10:59:47 MST 2003  
 Last updated     :   Fri Nov 11 02:41:57 MST 2005
   
 DESCRIPTION:
    Header file for global definition in biolib.
 
**********************************************************************************************************/

#ifndef BIOLIB_BIOLIB_INCLUDED
#define BIOLIB_BIOLIB_INCLUDED

#include <time.h>
#include <vector>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <ostream>
#include <iostream>
#include <assert.h>

//! Main namespace in biolib.
namespace biolib { 

/*! \brief Basic matrix object. 

 It is a typedef to vectos<vector<float>> type.
 */
typedef std::vector<std::vector<float> > matrix;

}

#endif 
