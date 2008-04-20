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
#include <trees/neighbour_joining_tree.h>
#include <trees/branch_proportional_weighting.h>

int main() 
{
  std::vector<std::vector<float> > dists(6, std::vector<float>(6,0));
  
  dists[0][1] = dists[1][0] =  5;
  dists[0][2] = dists[2][0] =  6;
  dists[0][3] = dists[3][0] = 11;
  dists[0][4] = dists[4][0] =  9;
  dists[0][5] = dists[5][0] =  9;

  dists[1][2] = dists[2][1] =  5;
  dists[1][3] = dists[3][1] = 10;
  dists[1][4] = dists[4][1] =  8;
  dists[1][5] = dists[5][1] =  8;

  dists[2][3] = dists[3][2] =  7;
  dists[2][4] = dists[4][2] =  5;
  dists[2][5] = dists[5][2] =  5;

  dists[3][4] = dists[4][3] =  8;
  dists[3][5] = dists[5][3] =  8;

  dists[4][5] = dists[5][4] =  2;

  std::cout << "Distances: " << std::endl;
  for(int i=0; i<dists.size(); i++) {
    for(int j=0; j<dists.size(); j++)
      std::cout << std::setw(4) << dists[i][j];
    std::cout << std::endl; 
  }  

  std::cout << std::endl;
  
  biolib::trees::neighbour_joining_tree tree(dists);

  std::cout << "Unrooted tree: " << std::endl;
  tree.print();

  tree.rooted();
  
  std::cout << "Rooted tree: " << std::endl;
  tree.print();

  biolib::trees::branch_proportional_weighting<int> weights(&tree);

  weights.solve();
 
  std::vector<float> const & w = weights.get_weights();
 
  for(int i=0; i<w.size(); i++) 
    std::cout << w[i] << std::endl;
  
  return 0;  
}

