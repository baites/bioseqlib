/*********************************************************************************************************

 Library name     :   Biolib++ (Rutines for biological applications)
 File name        :   alignment3.C
 Defined Objects  :   --  

 Programmer       :   Victor E. Bazterra 
 Directed by      :   Julio C. Facelli
                      
 Initialized      :   Thu Nov 13 10:17:20 MST 2003  
 Last updated     :   Thu Nov 13 10:17:20 MST 2003
   
 DESCRIPTION:
    Example of global alignment between two sequences.
 
**********************************************************************************************************/

#include <iostream>
#include <amino_acids.h>
#include <global_alignment/two_sequences_simple.h>

main() 
{  
  biolib::amino_acid_alphabet alphabet;
  
  biolib::alignment alignment(alphabet);
  
  alignment.add("SeqA","HEAGAWGHEE");
  alignment.add("SeqB","PAWHEAE"   );
  alignment.add("SeqC","VLSSDPADK" );
  alignment.add("SeqD","HLASAESK"  );
  
  std::cout << alignment;

  biolib::global_alignment::two_sequences_simple_gap alignerA(alignment);

  biolib::similarity sim(alphabet, biolib::blosum50mt, "BLOSUM 50");

  alignerA.align(sim,8);

  std::cout << alignerA.get_alignment();

  std::cout << "Score : " << alignerA.score(sim,8) << std::endl;

  alignerA.align_sequences(2,3,sim,8);

  std::cout << alignerA.get_alignment();

  std::cout << "Score : " << alignerA.score_sequences(2,3,sim,8) << std::endl;

  biolib::global_alignment::two_sequences_affine_gap alignerB(alignment);
  
  alignerB.align(sim,8,8);

  std::cout << alignerB.get_alignment();

  std::cout << "Score : " << alignerB.score(sim,8,8) << std::endl;

  alignerB.align_sequences(2,3,sim,8,8);

  std::cout << alignerB.get_alignment();

  std::cout << "Score : " << alignerB.score_sequences(2,3,sim,8,8) << std::endl;
  
  return 0;
}
