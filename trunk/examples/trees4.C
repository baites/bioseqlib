/*********************************************************************************************************

 Library name     :   Biolib++ (Rutines for biological applications)
 File name        :   trees4.C
 Defined Objects  :   --  

 Programmer       :   Victor E. Bazterra 
 Directed by      :   Julio C. Facelli
                      
 Initialized      :   Wed Nov 12 10:59:47 MST 2003  
 Last updated     :   Wed Nov 12 10:59:47 MST 2003
   
 DESCRIPTION:
    Example of left_rotation function defined binary_tree object.
 
**********************************************************************************************************/

#include <trees/binary_tree.h>

int main() 
{		
   biolib::trees::binary_tree<int,float> treeX(3,0.3);

   biolib::trees::binary_tree<int,float> treeY(2,0.2);

   treeY.left(treeX);
 
   treeX.initialize(1,0.1);
   
   treeX.right(treeY);
      
   treeX.print();
       
   treeX.preorder();
   
   treeX.left_rotation();
                  
   treeX.print();  
  
   return 0;  
}
