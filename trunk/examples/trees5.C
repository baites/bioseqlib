/*********************************************************************************************************

 Library name     :   Biolib++ (Rutines for biological applications)
 File name        :   trees5.C
 Defined Objects  :   --  

 Programmer       :   Victor E. Bazterra 
 Directed by      :   Julio C. Facelli
                      
 Initialized      :   Wed Nov 12 10:59:47 MST 2003  
 Last updated     :   Wed Nov 12 10:59:47 MST 2003
   
 DESCRIPTION:
    First example of rooted tree using neighbour_joining_tree algorithm.
 
**********************************************************************************************************/

#include <vector>
#include <iomanip>
#include <trees/neighbour_joining_tree.h>

int main() 
{

  std::vector<std::vector<float> > 
    dists(4, std::vector<float>(4, 0));
  
  dists[0][1] = dists[1][0] = 0.3;
  dists[0][2] = dists[2][0] = 0.5;
  dists[0][3] = dists[3][0] = 0.6;
  dists[1][2] = dists[2][1] = 0.6;
  dists[1][3] = dists[3][1] = 0.5;
  dists[2][3] = dists[3][2] = 0.9;

  std::cout << "Distances: " << std::endl;
  for(size_t i=0; i<dists.size(); i++) {
    for(size_t j=0; j<dists.size(); j++)
      std::cout << std::setw(4) << dists[i][j];
    std::cout << std::endl; 
  }  

  std::cout << std::endl;
  
  biolib::trees::neighbour_joining_tree tree(dists);

  tree.print();
  
  tree.rooted();

  tree.print();
  
  return 0;  
}
