#ifndef BIOLIB_GLOBAL_ALIGNMENT_TWO_SEQUENCES_BASE_INCLUDED
#define BIOLIB_GLOBAL_ALIGNMENT_TWO_SEQUENCES_BASE_INCLUDED

#include <global_alignment/base.h>

namespace biolib {  namespace global_alignment {

//! Abstract base class for simple gap aligment algorithms for two sequences. 
class two_sequences_simple_gap_base : public simple_gap_base {

public:

  //! Void constructor.
  two_sequences_simple_gap_base() 
    : simple_gap_base() 
  {
    _ids_out.resize(2);
    _align_out.resize(2); 
  }
  
  //! Constructor from a set of sequences.
  two_sequences_simple_gap_base(alignment const & align) 
    : simple_gap_base(align) 
  {
    _ids_out.resize(2);
    _align_out.resize(2);
  }
  
  //! Aligns two sequences (by defauls the first two of the group).
  void align (similarity const & sim, int gap)
  {
    align_sequences(0, 1, sim, gap); 
  }

  //! Scores the best alignment between two sequences (by default the first two of the group).
  long score (similarity const & sim, int gap)
  {
    return score_sequences(0, 1, sim, gap);
  }

  //! Aligns two arbitrary sequences from a group of sequences.
  virtual void align_sequences (size_t, size_t, similarity const &, int) = 0;

  //! Scores two arbitrary sequences from a group of sequences.
  virtual long score_sequences (size_t s1, size_t s2, similarity const & sim, int gap) = 0;

};

//! Abstract base class for affine gap aligment algorithms for two sequences.
class two_sequences_affine_gap_base : public affine_gap_base {

public:

  //! Void constructor.
  two_sequences_affine_gap_base() 
    : affine_gap_base() 
  {
    _ids_out.resize(2);
    _align_out.resize(2); 
  }
  
  //! Constructor from a set of sequences.
  two_sequences_affine_gap_base(alignment const & align) 
    : affine_gap_base(align) 
  {
    _ids_out.resize(2);
    _align_out.resize(2);
  }

  //! Aligns two sequences (by defauls the first two of the group).
  void align (similarity const & sim, int alpha, int beta)
  {
    align_sequences(0, 1, sim, alpha, beta); 
  }

  //! Scores the best alignment between two sequences (by default the first two of the group).
  long score (similarity const & sim, int alpha, int beta)
  {
    return score_sequences(0, 1, sim, alpha, beta);
  }

  //! Aligns two arbitrary sequences from a group of sequences.
  virtual void align_sequences (size_t, size_t, similarity const &, int, int) = 0;

  //! Scores two arbitrary sequences from a group of sequences.
  virtual long score_sequences (size_t s1, size_t s2, similarity const & sim, int alpha, int beta) = 0;  
};


}

}

#endif
