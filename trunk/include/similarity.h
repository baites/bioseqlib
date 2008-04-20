#ifndef BIOLIB_SIMILARITY_INCLUDED
#define BIOLIB_SIMILARITY_INCLUDED

#include <biolib.h>
#include <amino_acids.h>

namespace biolib {

//! Class similarity contains information about replacement matrices.
class similarity {

  std::string          _name;
  std::size_t          _translate[256];
  std::vector<short>   _score_matrix;
  alphabet::id_array_t _id_array;

public:

  //! Void constructor.
  similarity() { 
  	set_name("unknown"); 
  }

  //! Constructor by alphabet and replacement matrices.
  similarity (alphabet const & a, short* const m) {
	set_name("unknown");
  	set_similarity(a,m); 
  }

  //! Constructor by alphabet, replacement matrices and name.
  similarity (alphabet const & a, short* const m, std::string const & s) {
        set_name(s);
  	set_similarity(a,m); 
  }

  //! Copy constructor.
  similarity (similarity const &);

  //! Sets the name for the similarity object.
  void  set_name(std::string const & name)
  {
    assert(name.length() > 0);
    _name = name;
  }
  
  //! Sets the alphabet and replacement matrices.
  void  set_similarity(alphabet const &, short* const);

  //! Gets the similarity for two given symbols.
  short get_similarity(char, char) const;

  //! Gets the name of the similarity object.
  std::string const & get_id() const
  {
    return _name;
  }

  //! Print all the information of the similarity object.
  friend std::ostream& operator<< (std::ostream& o, similarity const & a);
  
};


}

#endif 
