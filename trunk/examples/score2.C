/*********************************************************************************************************

 Library name     :   Biolib++ (Rutines for biological applications)
 File name        :   alignment5.C
 Defined Objects  :   --  

 Programmer       :   Victor E. Bazterra 
 Directed by      :   Julio C. Facelli
                      
 Initialized      :   Thu Nov 13 10:17:20 MST 2003  
 Last updated     :   Thu Nov 13 10:17:20 MST 2003
   
 DESCRIPTION:
    Example of distances estimation between sequences and their related tree.
 
**********************************************************************************************************/

#include <iostream>
#include <amino_acids.h>
#include <global_alignment/two_sequences_simple.h>
#include <scores/alignment_score.h>

main() {
  
  biolib::amino_acid_alphabet alphabet;

  biolib::alignment alignment(alphabet);
  
  alignment.add("SeqA","HEAGAWGHEE");
  alignment.add("SeqB","PAWHEAE"   );

  std::cout << alignment;

  biolib::similarity similarity(alphabet, biolib::blosum50mt, "BLOSUM 50");
      
  biolib::global_alignment::two_sequences_affine_gap aligner(alignment);
 
  aligner.align(similarity, 8, 8);

  std::cout << aligner.score(similarity, 8, 8) << std::endl;

  std::cout <<  aligner.get_alignment();

  biolib::alignment nalignment(aligner.get_alignment()); 

  std::cout << nalignment << std::endl;

  biolib::scores::alignment_score score(similarity, 8, 8);

  score.evaluate(nalignment);
}

