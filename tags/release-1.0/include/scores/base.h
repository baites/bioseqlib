#ifndef BIOLIB_SCORES_BASE_INCLUDED
#define BIOLIB_SCORES_BASE_INCLUDED

#include <biolib.h>
#include <alignment.h>
#include <similarity.h>

namespace biolib {  
//! Contains all score algorithms.
namespace scores {

//! Base class for all score algorithms.
class base {

public:

  //! Void constructor.
  base() 
  {
    _parameters_flag = false;  
  }

  //! Constructor by similarity and the penalization for open and extend gaps.
  base(similarity const & sim, int alpha, int beta)
  {
    set_parameters(sim, alpha, beta);
  }

  //! Copy constructor.
  base(base const & original)
  { 
    _alpha = original._alpha;
    _beta = original._beta;
    _similarity = original._similarity;
  } 

  //! Sets the similarity and the penalization for open and extend gaps.   
  void set_parameters(similarity const & sim, int alpha, int beta) 
  {
    _similarity = sim;
    _alpha = alpha;
    _beta = beta;
    _parameters_flag = true;
  }
  
  //! Evaluates the score of an aligment.  
  virtual float evaluate(alignment const &) const = 0;

protected:

  bool _parameters_flag;

  int _alpha, _beta;
  similarity _similarity;

};


}


}

#endif
