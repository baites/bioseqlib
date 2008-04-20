#ifndef BIOLIB_GLOBAL_ALIGNMENT_TWO_SEQUENCES_LINEAR_INCLUDED
#define BIOLIB_GLOBAL_ALIGNMENT_TWO_SEQUENCES_LINEAR_INCLUDED

#include <global_alignment/two_sequences_base.h>


namespace biolib {  namespace global_alignment {

//! This class implements a simple gap aligment algorithm for two sequences in O(n).
class two_sequences_simple_gap_linear : public two_sequences_simple_gap_base {

  std::vector<long> Sf, Sb;

  void global_block (size_t, size_t, size_t, 
                     size_t, size_t, size_t,
		     similarity const &, int); 

  //! Global score forward in O(n).
  long global_score_forward  (size_t, size_t, size_t, 
  		              size_t, size_t, size_t,
			      similarity const &, int);

  //! Global score backward in O(n).
  long global_score_backward (size_t, size_t, size_t, 
  			      size_t, size_t, size_t,
		              similarity const &, int); 

public:

  //! Void constructor.
  two_sequences_simple_gap_linear() 
    : two_sequences_simple_gap_base() {}
  
  //! Constructor from a set of sequences.
  two_sequences_simple_gap_linear(alignment const & align) 
    : two_sequences_simple_gap_base(align) {} 

  //! Aligns two arbitrary sequences from a group of sequences.  
  void align_sequences (size_t, size_t, similarity const &, int);

  //! Scores two arbitrary sequences from a group of sequences.
  long score_sequences (size_t s1, size_t s2, similarity const & sim, int gap);

};

//! This class implements a affine gap aligment algorithm for two sequences in O(n).
class two_sequences_affine_gap_linear : public two_sequences_affine_gap_base {

  std::vector<long> Sf, Ff, Sb, Fb;

  void global_block (size_t, size_t, size_t, 
                     size_t, size_t, size_t,
		     similarity const &, int, int, int, int); 

  //! Global score forward in O(n).
  long global_score_forward  (size_t, size_t, size_t, 
  		              size_t, size_t, size_t,
			      similarity const &, int, int, int);

  //! Global score backward in O(n).
  long global_score_backward (size_t, size_t, size_t, 
  			      size_t, size_t, size_t,
		              similarity const &, int, int, int); 

public:

  //! Void constructor.
  two_sequences_affine_gap_linear() 
    : two_sequences_affine_gap_base() {}
  
  //! Constructor from a set of sequences.
  two_sequences_affine_gap_linear(alignment const & align) 
    : two_sequences_affine_gap_base(align) {} 
  
  //! Aligns two arbitrary sequences from a group of sequences.
  void align_sequences (size_t, size_t, similarity const &, int, int);

  //! Scores two arbitrary sequences from a group of sequences.
  long score_sequences (size_t s1, size_t s2, similarity const & sim, int alpha, int beta);

};


}

}

#endif

