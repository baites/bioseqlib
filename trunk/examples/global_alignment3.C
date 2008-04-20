/*********************************************************************************************************

 Library name     :   Biolib++ (Rutines for biological applications)
 File name        :   alignment3.C
 Defined Objects  :   --  

 Programmer       :   Victor E. Bazterra 
 Directed by      :   Julio C. Facelli
                      
 Initialized      :   Thu Nov 13 10:17:20 MST 2003  
 Last updated     :   Thu Nov 13 10:17:20 MST 2003
   
 DESCRIPTION:
    Example of global alignment in linear space.
 
**********************************************************************************************************/

#include <iostream>
#include <io.h>
#include <amino_acids.h>
#include <global_alignment/two_sequences_simple.h>
#include <global_alignment/two_sequences_linear.h>

main() {

  biolib::amino_acid_alphabet alphabet;

  biolib::alignment alignment(alphabet);

  biolib::load("test2.msf", alignment);  
 
  biolib::global_alignment::two_sequences_simple_gap alignerA(alignment);

  biolib::global_alignment::two_sequences_simple_gap_linear lalignerA(alignment);

  biolib::similarity sim(alphabet, biolib::blosum50mt, "BLOSUM 50");

  alignerA.align_sequences(0,3,sim,8);

  std::cout << alignerA.get_alignment();

  lalignerA.align_sequences(0,3,sim,8);
  
  std::cout << lalignerA.get_alignment();

  biolib::global_alignment::two_sequences_affine_gap alignerB(alignment);

  biolib::global_alignment::two_sequences_affine_gap_linear lalignerB(alignment);

  alignerB.align_sequences(0,3,sim,8,8);

  std::cout << alignerA.get_alignment();

  lalignerB.align_sequences(0,3,sim,8,8);

  std::cout << lalignerB.get_alignment();

}
