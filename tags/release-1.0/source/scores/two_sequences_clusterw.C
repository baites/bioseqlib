#include <scores/two_sequences_clusterw.h>


namespace biolib {  namespace scores {


float two_sequences_clusterw::evaluate_sequences(
    std::size_t indxA, std::size_t indxB,
    alignment const & align) const
{
  assert(_parameters_flag && align.size() > 1); 

  global_alignment::two_sequences_affine_gap_linear aligner(align);
  
  aligner.align_sequences(indxA, indxB, _similarity, _alpha, _beta); 

  alignment naling = aligner.get_alignment();
    
  alignment::data_array_t const & seq = naling.get_data_array();

  assert(seq[0].size() == seq[1].size()); 

  std::size_t size = seq[0].size();

  float count = 0;

  for (std::size_t i=0; i<size; i++) 
  {
    if (seq[0][i] == '-' || seq[1][i] == '-') continue;
      
    if (seq[0][i] == seq[1][i]) count++;
  }

  return 1.0 - count/size; 
}


}


}
