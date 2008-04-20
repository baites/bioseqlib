/*********************************************************************************************************

 Library name     :   Biolib++ (Rutines for biological applications)
 File name        :   alignment2.C
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

main() {
  
  biolib::amino_acid_alphabet alphabet;
  
  biolib::alignment alignment(alphabet);

  alignment.add("SeqA" ,"HAEGAWGHEE");
  alignment.add("SeqB" ,"PAWHEAE"   );

  std::cout << alignment;

  biolib::alignment::id_array_t ids(alignment.get_id_array());
  biolib::alignment::data_array_t seqs(alignment.get_data_array());

  alignment.remove("SeqB");

  std::cout << alignment;  
  
  alignment.set(ids,seqs);

  std::cout << alignment;

  biolib::alignment alignmentB(ids,seqs,alphabet);

  std::cout << alignment;
}
