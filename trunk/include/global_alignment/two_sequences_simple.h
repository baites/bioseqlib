#ifndef BIOLIB_GLOBAL_ALIGNMENT_TWO_SEQUENCES_SIMPLE_INCLUDED
#define BIOLIB_GLOBAL_ALIGNMENT_TWO_SEQUENCES_SIMPLE_INCLUDED

#include <global_alignment/two_sequences_base.h>


namespace biolib {  namespace global_alignment {

//! This class implements a simple gap aligment algorithm for two sequences in O(n^2).
class two_sequences_simple_gap : public two_sequences_simple_gap_base {

  std::vector<std::vector<long> > D; 

  //! Global score in O(n^2).
  long global_score (size_t, size_t, size_t, 
                     size_t, size_t, size_t,
		     similarity const &, int); 

  //! Traceback aligment reconstruction.
  void traceback_aligment (size_t, size_t, size_t, 
                           size_t, size_t, size_t,
		           similarity const &, int); 
public:

  //! Void constructor.
  two_sequences_simple_gap() 
    : two_sequences_simple_gap_base() {}
  
  //! Constructor from a set of sequences.
  two_sequences_simple_gap(alignment const & align) 
    : two_sequences_simple_gap_base(align) {}
  
  //! Aligns two arbitrary sequences from a group of sequences.
  void align_sequences (size_t, size_t, similarity const &, int);

  //! Scores two arbitrary sequences from a group of sequences. 
  long score_sequences (size_t s1, size_t s2, similarity const & sim, int gap)
  {
    return global_score(s1, s2, 0, 0, _align_in[s1].size(), _align_in[s2].size(), sim, gap);  
  }

};

//! This class implements a affine gap aligment algorithm for two sequences in O(n^2).
class two_sequences_affine_gap : public two_sequences_affine_gap_base {

  std::vector<std::vector<long> > D, E, F; 

  //! Global score in O(n^2).
  long global_score (size_t, size_t, size_t, 
                     size_t, size_t, size_t,
		     similarity const &, int, int); 

  //! Traceback aligment reconstruction.
  void traceback_aligment (size_t, size_t, size_t, 
                           size_t, size_t, size_t,
		           similarity const &, int, int); 
public:

  //! Void constructor.
  two_sequences_affine_gap() 
    : two_sequences_affine_gap_base() {}
  
  //! Constructor from a set of sequences.
  two_sequences_affine_gap(alignment const & align) 
    : two_sequences_affine_gap_base(align) {}
  
  //! Aligns two arbitrary sequences from a group of sequences.
  void align_sequences (size_t, size_t, similarity const &, int, int);

  //! Scores two arbitrary sequences from a group of sequences.
  long score_sequences (size_t s1, size_t s2, similarity const & sim, int alpha, int beta)
  {
    return global_score(s1, s2, 0, 0, _align_in[s1].size(), _align_in[s2].size(), sim, alpha, beta);  
  }

}; 


}

}

#endif
