#include <iostream>
#include <similarity.h>
#include <amino_acids.h>

main()
{
  biolib::amino_acid_alphabet alphabet;
  biolib::similarity similarity(alphabet, biolib::blosum50mt, "BLOSUM 50");
  
  std::cout << alphabet;
  std::cout << similarity;
  
  std::cout << "The score for the pair H,E is: ";
  std::cout << similarity.get_similarity('H','E') << std::endl; 

  return 0;
}
