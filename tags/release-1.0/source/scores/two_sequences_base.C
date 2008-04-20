#include <scores/two_sequences_base.h>


namespace biolib {  namespace scores {


matrix two_sequences_base::evaluate_matrix(alignment const & align) const
{
  matrix distances;   
 
  std::size_t nseqs = align.size();
 
  distances.resize(nseqs, std::vector<float>(nseqs,0));

  for (std::size_t i=0; i<nseqs; i++) 
    for (std::size_t j=i+1; j<nseqs; j++)  
      distances[i][j] = distances[j][i] = evaluate_sequences(i, j, align); 

  return distances; 
}


}


}
