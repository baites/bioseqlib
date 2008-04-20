/*********************************************************************************************************

 Library name     :   Biolib++ (Rutines for biological applications)
 File name        :   trees7.C
 Defined Objects  :   --  

 Programmer       :   Victor E. Bazterra 
 Directed by      :   Julio C. Facelli
                      
 Initialized      :   Wed Nov 12 10:59:47 MST 2003  
 Last updated     :   Wed Nov 12 10:59:47 MST 2003
   
 DESCRIPTION:
    A examples of weighting a tree using branch_proportional_weighting object.
 
**********************************************************************************************************/

#include <iostream>
#include <vector>
#include <iomanip>
#include <alignment.h>
#include <amino_acids.h>
#include <scores/alignment_score.h>
#include <scores/two_sequences_clusterw.h>
#include <trees/neighbour_joining_tree.h>
#include <trees/branch_proportional_weighting.h>

int main() 
{
  biolib::amino_acid_alphabet alphabet;

  biolib::alignment alignment(alphabet);
  biolib::alignment nalignment(alphabet);


  alignment.add("SeqA","HEAGAWGHEE");
  alignment.add("SeqB","PAWHEAE"   );
  alignment.add("SeqC","PAWGAWE"   );

  nalignment.add("SeqA","HEAGAWGHEE");
  nalignment.add("SeqB","PAW-H-E-AE");
  nalignment.add("SeqC","PAWG---AWE");

  std::cout << alignment;
  std::cout << nalignment;

  biolib::similarity similarity(alphabet, biolib::blosum50mt, "BLOSUM 50");

  biolib::scores::two_sequences_clusterw scoring(similarity, 8, 8);  
  biolib::matrix distances(scoring.evaluate_matrix(alignment));
  
  std::cout << "Distances: " << std::endl;
  for(int i=0; i<distances.size(); i++) {
    for(int j=0; j<distances[i].size(); j++)
      std::cout << std::setw(8) << distances[i][j];
    std::cout << std::endl; 
  }  

  biolib::trees::neighbour_joining_tree tree(distances);

  std::cout << "Helo" << std::endl;

  tree.rooted();
  tree.print();

  biolib::trees::branch_proportional_weighting<int> weights(&tree);
  weights.solve();
  
  std::vector<float> const & w = weights.get_weights();
  for(int i=0; i<w.size(); i++) 
    std::cout << w[i] << std::endl;

  biolib::scores::alignment_score score(similarity, 12, 2);

  score.set_weights(w);
  score.evaluate(nalignment);

  return 0;  
}

