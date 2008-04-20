/*********************************************************************************************************

 Library name     :   Biolib++ (Rutines for biological applications)
 File name        :   alignment3.C
 Defined Objects  :   --  

 Programmer       :   Victor E. Bazterra 
 Directed by      :   Julio C. Facelli
                      
 Initialized      :   Thu Nov 13 10:17:20 MST 2003  
 Last updated     :   Thu Nov 13 10:17:20 MST 2003
   
 DESCRIPTION:
    Example of global alignmentment in linear space.
 
**********************************************************************************************************/

#include <iostream>
#include <amino_acids.h>
#include <global_alignment/two_sequences_simple.h>
#include <global_alignment/two_sequences_linear.h>

main() {
  
  biolib::amino_acid_alphabet alphabet;

  biolib::alignment alignment(alphabet);
  
  alignment.add("SeqA","HEAGAWGHEE");
  alignment.add("SeqB","PAWHEAE"   );
  alignment.add("SeqC","VLSSDPADK" );
  alignment.add("SeqD","HLASAESK"  );

  std::cout << alignment;
    
  std::cout << std::endl;
  std::cout << "*** Using simple algorithm ***" << std::endl;
  
  biolib::similarity similarity(alphabet, biolib::blosum50mt, "BLOSUM 50");

  biolib::global_alignment::two_sequences_simple_gap alignerA(alignment);

  alignerA.align(similarity,8);

  std::cout << "Score :" << alignerA.score(similarity,8) << std::endl;
  std::cout << alignerA.get_alignment();

  alignerA.align_sequences(2,3,similarity,8);

  std::cout << "Score :" << alignerA.score_sequences(2,3,similarity,8) << std::endl;
  std::cout << alignerA.get_alignment();  

  biolib::global_alignment::two_sequences_affine_gap alignerB(alignment);

  alignerB.align(similarity,12,2);
  std::cout << "Score :" << alignerB.score(similarity,12,2) << std::endl;

  std::cout << alignerB.get_alignment();

  std::cout << "Score :" << alignerB.score_sequences(2,3,similarity,12,2) << std::endl;
  alignerB.align_sequences(2,3,similarity,12,2);

  std::cout << alignerB.get_alignment(); 

 
  std::cout << std::endl;
  std::cout << "*** Using linear-space algorithm ***" << std::endl;


  biolib::global_alignment::two_sequences_simple_gap_linear alignerC(alignment);

  alignerC.align(similarity,8);

  std::cout << "Score :" << alignerC.score(similarity,8) << std::endl;
  std::cout << alignerC.get_alignment();

  alignerC.align_sequences(2,3,similarity,8);

  std::cout << "Score :" << alignerC.score_sequences(2,3,similarity,8) << std::endl;
  std::cout << alignerC.get_alignment();  

  biolib::global_alignment::two_sequences_affine_gap_linear alignerD(alignment);

  alignerD.align(similarity,12,2);
  std::cout << "Score :" << alignerD.score(similarity,12,2) << std::endl;

  std::cout << alignerD.get_alignment();

  std::cout << "Score :" << alignerD.score_sequences(2,3,similarity,12,2) << std::endl;
  alignerD.align_sequences(2,3,similarity,12,2);

  std::cout << alignerD.get_alignment(); 

}
