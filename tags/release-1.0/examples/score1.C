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
#include <scores/two_sequences_clusterw.h>

main() {
  
  biolib::amino_acid_alphabet alphabet;

  biolib::alignment alignment(alphabet);
  
  alignment.add("SeqA","HEAGAWGHEE");
  alignment.add("SeqB","PAWHEAE"   );
  alignment.add("SeqC","VLSSDPADK" );
  alignment.add("SeqD","HLASAESK"  );

  std::cout << alignment;

  biolib::similarity similarity(alphabet, biolib::blosum50mt, "BLOSUM 50");
      
  biolib::scores::two_sequences_clusterw dist(similarity, 8, 8);  

  std::cout << "Distance between seq. 0 and 1 is " << dist.evaluate(alignment) << std::endl;

  biolib::matrix mat = dist.evaluate_matrix(alignment);

  std::cout.setf(std::ios::fixed);

  std::cout << "Distances between sequences : " << std::endl;
  for(std::size_t i=0; i<mat.size(); i++) {
    for(std::size_t j=0; j<mat[i].size(); j++)
      std::cout << std::setw(6) << std::setprecision(2) << mat[i][j];
    std::cout << std::endl;
  }    

}

