#ifndef BIOLIB_SCORES_ALIGNMENT_SCORE_INCLUDED
#define BIOLIB_SCORES_ALIGNMENT_SCORE_INCLUDED

#include <scores/base.h>

namespace biolib {  namespace scores {

//! Score of a given sequence aligment.
class alignment_score : public base {

  std::vector<float> _weights;

public:

  //! Void constructor.
  alignment_score()
    :base() {}

  //! Constructor by similarity and the penalization for open and extend gaps. 
  alignment_score(similarity const & sim, int alpha, int beta)
    :base(sim, alpha, beta) {}
 
  //! Copy constructor.
  alignment_score(alignment_score const & original)  
    :base(original) {}

  //! Sets sequence weights.
  void set_weights(std::vector<float> const & weights)
  {
     _weights = weights;
  }

  //! Evaluates score of an aligment.
  float evaluate(alignment const & align) const;
};



}


}

#endif

