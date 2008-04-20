/*********************************************************************************************************

 Library name     :   Biolib++ (Rutines for biological applications)
 File name        :   amino_acids.h
 Defined Objects  :   amino_acid_alphabet  

 Programmer       :   Victor E. Bazterra 
 Directed by      :   Julio C. Facelli
                      
 Initialized      :   Wed Nov 12 10:59:47 MST 2003  
 Last updated     :   Wed Nov 12 10:59:47 MST 2003
   
 DESCRIPTION:
    It defines amino acid alphabet and its sustitution matrices.
 
**********************************************************************************************************/

#ifndef BIOLIB_AMINO_ACIDS_ALPHABET_INCLUDED 
#define BIOLIB_AMINO_ACIDS_ALPHABET_INCLUDED 

#include <biolib.h>
#include <alphabet.h>

namespace biolib {

/*! \brief Amino acid alphabet. 

 This class defines the amino acid alphabet and its sustitution matrices. The complete list 
 of matrices is:
 
   - blosum30mt
   - blosum35mt 
   - blosum40mt  
   - blosum45mt  
   - blosum50mt  
   - blosum55mt  
   - blosum62mt  
   - blosum62mt2  
   - blosum65mt  
   - blosum70mt  
   - blosum75mt  
   - blosum80mt  
   - blosum85mt  
   - blosum90mt  

   -  pam20mt  
   -  pam60mt  
   -  pam120mt  
   -  pam160mt  
   -  pam250mt  
   -  pam350mt  

   -  md_40mt  
   -  md_120mt  
   -  md_250mt  
   -  md_350mt  

   -  idmat  

   -  gon40mt  
   -  gon80mt  
   -  gon120mt  
   -  gon160mt  
   -  gon250mt  
   -  gon300mt  
   -  gon350mt  
  
 
 \note See Biological sequence analysis: Probabilistic models of proteins and nucleic acids. 
 R. Durbin, S. Eddy, A. Krogh and G. Mitchison. Cambridge University Press 1998.
 */
class amino_acid_alphabet : public alphabet {
public:
  // Void constructor.
  amino_acid_alphabet();
};

/*! \example alphabet.C. 
 Usage of alphabet and amino_acid_alphabet objects.   
 */

extern short blosum30mt[];
extern short blosum35mt[];
extern short blosum40mt[];
extern short blosum45mt[];
extern short blosum50mt[];
extern short blosum55mt[];
extern short blosum62mt[];
extern short blosum62mt2[];
extern short blosum65mt[];
extern short blosum70mt[];
extern short blosum75mt[];
extern short blosum80mt[];
extern short blosum85mt[];
extern short blosum90mt[];

extern short pam20mt[];
extern short pam60mt[];
extern short pam120mt[];
extern short pam160mt[];
extern short pam250mt[];
extern short pam350mt[];

extern short md_40mt[];
extern short md_120mt[];
extern short md_250mt[];
extern short md_350mt[];

extern short idmat[];

extern short gon40mt[];
extern short gon80mt[];
extern short gon120mt[];
extern short gon160mt[];
extern short gon250mt[];
extern short gon300mt[];
extern short gon350mt[];


}

#endif

