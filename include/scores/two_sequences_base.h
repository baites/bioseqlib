#ifndef BIOLIB_SCORES_TWO_SEQUENCES_BASE_INCLUDED
#define BIOLIB_SCORES_TWO_SEQUENCES_BASE_INCLUDED

#include <biolib.h>
#include <scores/base.h>

namespace biolib {  namespace scores {

//! Base class for all score algorithms.
class two_sequences_base : public base {

public:

  //! Void constructor.
  two_sequences_base()
    :base() {}

  //! Constructor by similarity and the penalization for open and extend gaps.
  two_sequences_base(similarity const & sim, int alpha, int beta)
    :base(sim, alpha, beta) {}

  //! Copy constructor.
  two_sequences_base(two_sequences_base const & original)  
    :base(original) {}

  //! Evaluates score between two sequences (by default the first two of the group).
  float evaluate(alignment const & align) const
  {
    return evaluate_sequences(0,1,align); 
  }
  
  //! Evaluates score between any two sequences of a group.
  virtual float evaluate_sequences(
    std::size_t, std::size_t ,
    alignment const & align) const = 0;

  //! Evaluates the score matrix between any two sequences of a group.
  matrix evaluate_matrix(alignment const &) const;    

};


}


}

#endif
