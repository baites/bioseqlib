#include <global_alignment/two_sequences_simple.h>


namespace biolib {  namespace global_alignment {


void two_sequences_simple_gap::align_sequences(size_t s1, size_t s2, similarity const & sim, int gap)
{
  assert(s1 < _align_in.size() && s2 < _align_in.size()); 

  _ga_flag = true;
  
  _ids_out[0] = _ids_in[s1];
  _ids_out[1] = _ids_in[s2];    
  
  _align_out[0].clear();
  _align_out[1].clear();    
 
  std::size_t maxsize = std::max(_align_in[s1].size(), _align_in[s2].size());
  
  _align_out[0].reserve(maxsize);
  _align_out[1].reserve(maxsize);    
  
  global_score(s1, s2, 0, 0, _align_in[s1].size(), _align_in[s2].size(), sim, gap);
  traceback_aligment(s1, s2, 0, 0, _align_in[s1].size(), _align_in[s2].size(), sim, gap);

  std::reverse(_align_out[0].begin(), _align_out[0].end());
  std::reverse(_align_out[1].begin(), _align_out[1].end());
}


long two_sequences_simple_gap::global_score(size_t s1, size_t s2,
				            size_t i1, size_t j1, 
				            size_t i2, size_t j2,
				            similarity const & sim, int gap)
{  
  std::size_t i = i2 - i1;   
  std::size_t j = j2 - j1;   

  std::vector<long> T(j+1,0);
  
  D.clear();
  D.resize(i+1, T);

  for (std::size_t l=1; l<=i; l++)
    D[l][0] = - gap * l;

  for (std::size_t k=1; k<=j; k++) 
    D[0][k] = - gap * k;

  for (std::size_t l=1; l<=i; l++) 
    for (std::size_t k=1; k<=j; k++) {
      D[l][k] = std::max(D[l-1][k-1] + sim.get_similarity(_align_in[s1][l+i1-1], _align_in[s2][k+j1-1]), 
		std::max(D[l-1][k] - gap, D[l][k-1] - gap));      
    }

  return D[i][j];
}


void two_sequences_simple_gap::traceback_aligment(size_t s1, size_t s2,
				                  size_t i1, size_t j1, 
				                  size_t i2, size_t j2,
				                  similarity const & sim, int gap)
{
  std::size_t i = i2 - i1;   
  std::size_t j = j2 - j1;   
  
  while(1) {

    if (i == 0 && j == 0)
       break;

    if (i == 0) {
       for (std::size_t h=j;; h--) {
         _align_out[0].push_back('-');
         _align_out[1].push_back(_align_in[s2][h+j1-1]);
         if(h == 1) break;
       }
       break;
    }

    if (j == 0) {
       for (std::size_t h=i;; h--) {
         _align_out[0].push_back(_align_in[s1][h+i1-1]);
         _align_out[1].push_back('-');
         if(h == 1) break;
       }
       break;
    }
       
    if (D[i][j] - D[i][j-1] == - gap)
    {
       _align_out[0].push_back('-');
       _align_out[1].push_back(_align_in[s2][j+j1-1]);
       j--;  
    } 
    else if (D[i][j] - D[i-1][j-1] == sim.get_similarity(_align_in[s1][i+i1-1], _align_in[s2][j+j1-1])) 
    {
       _align_out[0].push_back(_align_in[s1][i+i1-1]);
       _align_out[1].push_back(_align_in[s2][j+j1-1]);
       j--; i--; 
    } 
    else  
    {
       _align_out[0].push_back(_align_in[s1][i+i1-1]);
       _align_out[1].push_back('-');
       i--;  
    }

  }
}  


void two_sequences_affine_gap::align_sequences(size_t s1, size_t s2, similarity const & sim, int alpha, int beta)
{
  assert(s1 < _align_in.size() && s2 < _align_in.size()); 

  _ga_flag = true;
  
  _ids_out[0] = _ids_in[s1];
  _ids_out[1] = _ids_in[s2];    
  
  _align_out[0].clear();
  _align_out[1].clear();    
 
  std::size_t maxsize = std::max(_align_in[s1].size(), _align_in[s2].size());
  
  _align_out[0].reserve(maxsize);
  _align_out[1].reserve(maxsize);    
  
  global_score(s1, s2, 0, 0, _align_in[s1].size(), _align_in[s2].size(), sim, alpha, beta);
  traceback_aligment(s1, s2, 0, 0, _align_in[s1].size(), _align_in[s2].size(), sim, alpha, beta);

  std::reverse(_align_out[0].begin(), _align_out[0].end());
  std::reverse(_align_out[1].begin(), _align_out[1].end());
}


long two_sequences_affine_gap::global_score(size_t s1, size_t s2,
				            size_t i1, size_t j1, 
				            size_t i2, size_t j2,
				            similarity const & sim, int alpha, int beta)
{  
  std::size_t i = i2 - i1;   
  std::size_t j = j2 - j1;   

  std::vector<long> T(j+1,0);
  
  D.clear();
  D.resize(i+1, T);
  E.clear();
  E.resize(i+1, T);
  F.clear(); 
  F.resize(i+1, T); 

  for (std::size_t l=1; l<=i; l++)
    E[l][0] = D[l][0] = - alpha - beta*(l-1);

  for (std::size_t k=1; k<=j; k++) 
    F[0][k] = D[0][k] = - alpha - beta*(k-1);

  for (std::size_t l=1; l<=i; l++) 
    for (std::size_t k=1; k<=j; k++) {
      E[l][k] = std::max(D[l][k-1] - alpha, E[l][k-1] - beta);
      F[l][k] = std::max(D[l-1][k] - alpha, F[l-1][k] - beta);       

      D[l][k] = std::max(D[l-1][k-1] + sim.get_similarity(_align_in[s1][l+i1-1], _align_in[s2][k+j1-1]), 
		std::max(E[l][k], F[l][k]));      
    }

  return D[i][j];
}


void two_sequences_affine_gap::traceback_aligment(size_t s1, size_t s2,
				                  size_t i1, size_t j1, 
				                  size_t i2, size_t j2,
				                  similarity const & sim, int alpha, int beta)
{
  std::size_t i = i2 - i1;   
  std::size_t j = j2 - j1;   
  
  while(1) {

    if (i == 0 && j == 0)
       break;

    if (i == 0) {
       for (std::size_t h=j;; h--) {
         _align_out[0].push_back('-');
         _align_out[1].push_back(_align_in[s2][h+j1-1]);
         if(h == 1) break;
       }
       break;
    }

    if (j == 0) {
       for (std::size_t h=i;; h--) {
         _align_out[0].push_back(_align_in[s1][h+i1-1]);
         _align_out[1].push_back('-');
         if(h == 1) break;
       }
       break;
    }
       
    if (D[i][j] == E[i][j])
    {
       _align_out[0].push_back('-');
       _align_out[1].push_back(_align_in[s2][j+j1-1]);
       j--;  
    }
    else if (D[i][j] - D[i-1][j-1] == sim.get_similarity(_align_in[s1][i+i1-1], _align_in[s2][j+j1-1])) 
    {
       _align_out[0].push_back(_align_in[s1][i+i1-1]);
       _align_out[1].push_back(_align_in[s2][j+j1-1]);
       j--; i--; 
    } 
    else
    {
       _align_out[0].push_back(_align_in[s1][i+i1-1]);
       _align_out[1].push_back('-');
       i--;  
    } 
  }
}  


}


}

