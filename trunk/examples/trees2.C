/*********************************************************************************************************

 Library name     :   Biolib++ (Rutines for biological applications)
 File name        :   trees2.C
 Defined Objects  :   --  

 Programmer       :   Victor E. Bazterra 
 Directed by      :   Julio C. Facelli
                      
 Initialized      :   Wed Nov 12 10:59:47 MST 2003  
 Last updated     :   Wed Nov 12 10:59:47 MST 2003
   
 DESCRIPTION:
    Example of creation and releasing of a binary_tree object.
 
**********************************************************************************************************/

#include <trees/binary_tree.h>

int main() 
{
   biolib::trees::binary_tree<int,int> treeX(4,4);

   biolib::trees::binary_tree<int,int> treeY(3,3);

   treeY.left(treeX);
 
   treeX.initialize(1);
   
   treeX.right(treeY);
    
   treeY.initialize(2,2);

   treeX.left(treeY);
   
   treeX.print();   
 
   biolib::trees::binary_tree<int,int> treeZ(100);
 
   treeZ.release(treeX);

   treeZ.print();   
 
   if (treeX.empty())
     std::cout << "The tree is empty!" << std::endl;

   return 0;  
}
