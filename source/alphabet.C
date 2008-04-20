#include <alphabet.h>

namespace biolib {

std::ostream& operator<< (std::ostream& o, alphabet const & a)
{ 
  alphabet::id_array_t const & id_array = a.get_id_array();
  alphabet::data_array_t const & data_array = a.get_data_array();
  
  o << "'" << a.get_id() << "' alphabet:\n";
  for (std::size_t i=0; i<id_array.size(); i++) 
    o << id_array[i] << ") " << data_array[i] << "\n"; 

  return o;
}

}
