/*********************************************************************************************************

 Library name     :   Biolib++ (Rutines for biological applications)
 File name        :   trees7.C
 Defined Objects  :   --  

 Programmer       :   Victor E. Bazterra 
 Directed by      :   Julio C. Facelli
                      
 Initialized      :   Wed Nov 12 10:59:47 MST 2003  
 Last updated     :   Wed Nov 12 10:59:47 MST 2003
   
 DESCRIPTION:
 
**********************************************************************************************************/

#include <iostream>
#include <amino_acids.h>
#include <scores/alignment_score.h>

int main() 
{
  // Define a alphabet of aminoacids.
  biolib::amino_acid_alphabet alphabet;

  // Create a aligment object that can hold the alphabet.
  biolib::alignment alignment(alphabet);
  
  // Add a given aligment of two sequences.
  alignment.add("VSP1_AGKCO","VIGGDECNINEHRFLGMHNLKVLNKDALRRFPKEKYFCLN--------------------");
  alignment.add("VSP1_AGKRH","-IGGDECNINEHRFLVAVYEGTCARRRMNLVFGMHRKSEKFDDEQERYPKKRYFIRCNKT"); 

  // Print the aligment.
  std::cout << "Aligment to be evaluated:" << std::endl;
  std::cout << alignment;

  // Define the simalirity definition using BLOSUM 50 matrix.
  biolib::similarity similarity(alphabet, biolib::blosum50mt, "BLOSUM 50");

  // Create an object to evaluate the aligment score.
  biolib::scores::alignment_score score(similarity, 12, 2);

  // Evaluate the score.
  score.evaluate(alignment);

  return 0;  
}

