/*********************************************************************************************************

 Library name     :   Biolib++ (Rutines for biological applications)
 File name        :   trees1.C
 Defined Objects  :   --  

 Programmer       :   Victor E. Bazterra 
 Directed by      :   Julio C. Facelli
                      
 Initialized      :   Wed Nov 12 10:59:47 MST 2003  
 Last updated     :   Wed Nov 12 10:59:47 MST 2003
   
 DESCRIPTION:
    Example of creation and access to a binary_tree object.
 
**********************************************************************************************************/

#include <iostream>
#include <trees/binary_tree.h>

int main() 
{
   biolib::trees::binary_tree<int,int> treeX(4,4);

   biolib::trees::binary_tree<int,int> treeY(3,3);

   treeY.left(treeX);
 
   treeX.initialize(1,1);
   
   treeX.right(treeY);
    
   treeY.initialize(2,2);

   treeX.left(treeY);

   std::cout << "Direct access : " << std::endl;
   
   std::cout << std::endl;
   std::cout << "Data : " << treeX.data() << std::endl;
   std::cout << "Edge : " << treeX.edge() << std::endl;

   std::cout << std::endl;
   std::cout << "Moving to the left : " << std::endl;
   treeX.left();
   std::cout << std::endl;

   std::cout << "Data : " << treeX.data() << std::endl;
   std::cout << "Edge : " << treeX.edge() << std::endl;

   std::cout << std::endl;
   std::cout << "Moving in preorder : " << std::endl;
   std::cout << std::endl;
 
   for (int i=0; i<3; i++) {
     std::cout << "Data : " << treeX.data() << std::endl;
     std::cout << "Edge : " << treeX.edge() << std::endl;
     std::cout << std::endl;
     treeX.preorder();
   }
   
   std::cout << std::endl;
   std::cout << "Moving in postorder : " << std::endl;
   std::cout << std::endl;

   for (int i=0; i<8; i++) {
     std::cout << "Data : " << treeX.data() << std::endl;
     std::cout << "Edge : " << treeX.edge() << std::endl;
     std::cout << std::endl;
     treeX.postorder();
   } 

   biolib::trees::binary_tree<int,int>::iterator itr(treeX.begin());

   std::cout << "Moving using iterators : " << std::endl;
   std::cout << std::endl;
   
   for (int i=0; i<7; i++) {
     std::cout << "Data : " << (*itr).data() << std::endl;
     std::cout << "Edge : " << (*itr).edge() << std::endl;
     std::cout << std::endl;
     ++itr;
   }
   
   for (int i=0; i<8; i++) {
     std::cout << "Data : " << (*itr).data() << std::endl;
     std::cout << "Edge : " << (*itr).edge() << std::endl;
     std::cout << std::endl;
     --itr;
   }
   
   return 0;  
}
