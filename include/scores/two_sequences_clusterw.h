#ifndef BIOLIB_SCORES_TWO_SEQUENCES_CLUSTERW_INCLUDED
#define BIOLIB_SCORES_TWO_SEQUENCES_CLUSTERW_INCLUDED

#include <scores/two_sequences_base.h>
#include <global_alignment/two_sequences_linear.h>

namespace biolib {  namespace scores {

//! Evaluates distances between two sequences.
class two_sequences_clusterw : public two_sequences_base {

public:

  //! Void constructor.
  two_sequences_clusterw():two_sequences_base() {}

  //! Constructor by similarity and the penalization for open and extend gaps. 
  two_sequences_clusterw(similarity const & sim, int alpha, int beta)
    :two_sequences_base(sim, alpha, beta) {}

  //! Copy constructor.
  two_sequences_clusterw(two_sequences_clusterw const & original)  
    :two_sequences_base(original) {}
  
  //! Evaluates the two sequences distances.
  float evaluate_sequences(
    std::size_t, std::size_t ,
    alignment const & align) const;
  };


}


}

#endif
