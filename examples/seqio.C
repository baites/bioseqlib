#include <io.h>
#include <alignment.h>
#include <amino_acids.h>

main() {
  
  biolib::amino_acid_alphabet alphabet;
  
  biolib::alignment alignment(alphabet);
  
  biolib::load("test2.msf", alignment);

  std::cout << alignment;

  return 0;  
}
