/*********************************************************************************************************

 Library name     :   Biolib++ (Rutines for biological applications)
 File name        :   alignment1.C
 Defined Objects  :   --  

 Programmer       :   Victor E. Bazterra 
 Directed by      :   Julio C. Facelli
                      
 Initialized      :   Thu Nov 13 10:17:20 MST 2003  
 Last updated     :   Thu Nov 13 10:17:20 MST 2003
   
 DESCRIPTION:
    Example of sequence object and its container alignment.
 
**********************************************************************************************************/

#include <iostream>
#include <alignment.h>
#include <amino_acids.h>

main() 
{  
  biolib::amino_acid_alphabet alphabet;
  
  biolib::alignment alignment(alphabet);
  
  alignment.add("SeqA" ,"HAEGAWGHEE");
  alignment.add("SeqB" ,"PAWHEAE"   );

  std::cout << alignment;  
}
