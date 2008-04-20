#ifndef BIOLIB_GLOBAL_ALIGNMENT_BASE_INCLUDED
#define BIOLIB_GLOBAL_ALIGNMENT_BASE_INCLUDED

#include <biolib.h>
#include <alignment.h>
#include <similarity.h>

namespace biolib { 
//! Contains all global aligment algorithms.
namespace global_alignment {

//! Base class for all global aligment algorithms.
class base {

public:

  //! Void contructor.
  base()
  {
    _name = "global";
    _ga_flag = false;    
  }

  //! Copy constructor.
  base(alignment const & align)
  {
    _name = "global";
    _ga_flag = false;    
    set_alignment(align);
  }
  
  //! Sets a group of sequences for aligment.
  void set_alignment(alignment const & align)
  { 
    _alphabet = align.get_alphabet();
    _ids_in = align.get_id_array();
    _align_in  = align.get_data_array();
  }

  //! Gets a group of aligned sequences.
  alignment get_alignment() const
  {
    assert(_ga_flag == true); 
    return alignment(_name, _ids_out, _align_out, _alphabet);
  }

protected:

  bool _ga_flag;

  std::string _name;

  alphabet _alphabet;

  alignment::id_array_t   _ids_in;
  alignment::id_array_t   _ids_out;

  alignment::data_array_t _align_in;
  alignment::data_array_t _align_out; 

};

//! Abstract base class for simple gap aligment algorithms.
class simple_gap_base : public base {

public:

  //! Void constructor.
  simple_gap_base() : base() {}

  //! Copy constructor from a set of sequences.
  simple_gap_base(alignment const & align) 
    : base(align) {} 
 
  //! Aligns a given group of sequences.
  virtual void align (similarity const &, int) = 0;

  //! Score of the best aligment of a group of sequences.
  virtual long score (similarity const &, int) = 0;
  
};

//! Abstract base class for affine gap aligment algorithms.
class affine_gap_base : public base {

public:

  //! Void constructor.
  affine_gap_base() : base() {}

  //! Copy constructor from a set of sequences.
  affine_gap_base(alignment const & align)
    : base(align) {} 
 
  //! Aligns a given group of sequences.
  virtual void align (similarity const &, int, int) = 0;

  //! Score of the best aligment of a group of sequences.
  virtual long score (similarity const &, int, int) = 0;
  
};


}

}

#endif
