/*********************************************************************************************************

 Library name     :   Biolib++ (Rutines for biological applications)
 File name        :   alphabet.h
 Defined Objects  :   alphabet  

 Programmer       :   Victor E. Bazterra 
 Directed by      :   Julio C. Facelli
                      
 Initialized      :   Wed Nov 12 10:59:47 MST 2003  
 Last updated     :   Wed Nov 12 10:59:47 MST 2003
   
 DESCRIPTION:
    Class to define an alphabet.
 
**********************************************************************************************************/

#ifndef BIOLIB_ALPHABET_INCLUDED
#define BIOLIB_ALPHABET_INCLUDED

#include <biolib.h>
#include <container.h>

namespace biolib {

/*! \brief Class to define an alphabet. 

 It is a specialitation of a container class by typedef statement.
 */
typedef container<char, std::string> alphabet;

//! Overload operator<< for printing an alphabet.  
std::ostream& operator<< (std::ostream& o, alphabet const & a);

}

#endif
