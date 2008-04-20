#ifndef BIOLIB_ALIGMENT_INCLUDED
#define BIOLIB_ALIGMENT_INCLUDED

#include <biolib.h>
#include <alphabet.h>

namespace biolib {

//! An aligment object is defined as a container of sequences.
class alignment : public container<std::string, std::string> {

  typedef container<std::string, std::string> alignment_base;

  alphabet _alphabet; 

  bool _alphabet_flag;

public:

  //! Void constructor.
  alignment() : alignment_base() 
  {
    _alphabet_flag == false;
  }

  //! Constructor by name.
  alignment(std::string const & name) : alignment_base(name) 
  {
    _alphabet_flag == false;
  }

  //! Constructor by an alphabet.
  alignment(alphabet const & alph) 
    : alignment_base()
  {
    set_alphabet(alph);
  }

  //! Constructor by name and alphabet.
  alignment(std::string const & name, alphabet const & alph) 
    : alignment_base(name)
  {
    set_alphabet(alph);
  }
 
  //! Constructor by id, data array and alphabet.
  alignment(id_array_t const & id_array, 
            data_array_t const & data_array,
	    alphabet const & alph)
    : alignment_base(id_array, data_array)
  {
    set_alphabet(alph);
  }

  //! Constructor by name, id, data array and alphabet.
  alignment(std::string const & name,
            id_array_t const & id_array, 
            data_array_t const & data_array,
	    alphabet const & alph)  
    : alignment_base(name, id_array, data_array)
  {
    set_alphabet(alph);
  }

  //! Sets the alphabet of the sequences.
  void set_alphabet(alphabet alph)
  {
    assert(alph.size() > 0);
    _alphabet = alph;
    _alphabet_flag = true;
  }

  //! Gets the alphabet of the sequences.
  alphabet get_alphabet() const
  {
    return _alphabet;
  }

  //! Verifies if all the symbols of the sequences belong to the alphabet.
  bool check(id_t const &, data_t const &) const; 

  //! Print all the sequences.
  friend std::ostream& operator<< (std::ostream & o, alignment const & a);

}; 

 
}

#endif
