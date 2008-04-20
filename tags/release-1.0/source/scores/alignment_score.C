#include <scores/alignment_score.h>


namespace biolib {  namespace scores {


float alignment_score::evaluate(alignment const & align) const
{
  //Initialization
  std::size_t totalgaps = 0;

  float score = 0;
  float totalscore = 0;
  std::size_t gaps = 0;
    
  // Score evaluation
  char ibefore, kibefore;

  std::size_t alignsize = align.size();

  alignment::data_array_t const & seqs = align.get_data_array();

  std::size_t seqssize  = seqs[0].size();

  for (std::size_t i = 1; i < alignsize; i++)
    assert(seqs[i].size() == seqs[i-1].size());

  for (std::size_t i = 0; i < alignsize; i++) {
    for (std::size_t ki = i + 1; ki < alignsize; ki++) {
      for (int j = 0; j < seqssize; j++) {
           if (seqs[i][j] != '-' && seqs[ki][j] != '-' )
		score  += _similarity.get_similarity(seqs[i][j], seqs[ki][j]);
           // Nancy score function (not penalize initial and final gaps)
	   // else if  (seqs[i][j] == '-' && seqs[ki][j] != '-')
           //    score  -= _beta;
           // else if  (seqs[i][j] != '-' && seqs[ki][j] == '-')
           //    score  -= _beta;

           if (j == 0) {
           	ibefore   = seqs[i][j];
           	kibefore  = seqs[ki][j];
           	continue;
           } else if ( seqs[i][j] == '-' && seqs[ki][j] == '-'
            			         && j != seqssize - 1)
	       continue;
           else if ( ibefore  != '-' && seqs[i][j]  == '-' &&
                     kibefore != '-' && seqs[ki][j] != '-'
                                     && j != seqssize - 1)
               gaps ++;
           else if ( ibefore  != '-' && seqs[i][j]  != '-' &&
                     kibefore != '-' && seqs[ki][j] == '-'
                                     && j != seqssize - 1)
               gaps ++;
           else if ( ibefore  != '-' && seqs[i][j]  == '-' &&
                     kibefore == '-' && seqs[ki][j] != '-'
                                     && j != seqssize - 1)
               gaps ++;
	   else if ( ibefore  == '-' && seqs[i][j]  != '-' &&
                     kibefore != '-' && seqs[ki][j] == '-'
                                     && j != seqssize - 1)
               gaps ++;
           else if ((ibefore  == '-' && seqs[i][j]  == '-' ||
                     kibefore == '-' && seqs[ki][j] == '-' )
                                     && j == seqssize - 1)
               gaps --;
           ibefore   = seqs[i][j];
           kibefore  = seqs[ki][j];
        }
        score -= gaps * (_alpha - _beta);
	if(_weights.empty()) 
	  totalscore += score;
	else
	  totalscore += _weights[i]*_weights[ki]*score;	  
	totalgaps += gaps;
        score = 0.0;
        gaps = 0;
      }
  }

#ifdef BIOLIB_DEBUG
  std::cout.setf(std::ios::fixed);
  std::cout << " *** ALIGN EVALUATION *** " << std::endl  << std::endl;
  std::cout << "SCORE                   : " << totalscore << std::endl;
  std::cout << "NUMBER OF GAPS          : " << totalgaps  << std::endl;
  std::cout << "GAP OPEN PENALTY        : " << _alpha << std::endl;
  std::cout << "GAP EXTENTION PENALTY   : " << _beta << std::endl;
#endif

  return totalscore;
}


}


}

