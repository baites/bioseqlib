/*********************************************************************************************************

 Library name     :   Biolib++ (Rutines for biological applications)
 File name        :   trees3.C
 Defined Objects  :   --  

 Programmer       :   Victor E. Bazterra 
 Directed by      :   Julio C. Facelli
                      
 Initialized      :   Wed Nov 12 10:59:47 MST 2003  
 Last updated     :   Wed Nov 12 10:59:47 MST 2003
   
 DESCRIPTION:
    Example of copy and root function defined binary_tree object.
 
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
   
   biolib::trees::binary_tree<int,int> clone(treeX);

   std::cout << "Original Tree" << std::endl;
   treeX.print();

   std::cout << "Copy Tree" << std::endl;
   clone.print();

// Insertion of a root

   treeX.postorder();

   treeX.root();

   std::cout << "Root insertion 1" << std::endl;
   treeX.print();

// Insertion of a root

   treeX.left();

   treeX.root(5);

   std::cout << "Root insertion 2" << std::endl;
   treeX.print();

   return 0;  
}
