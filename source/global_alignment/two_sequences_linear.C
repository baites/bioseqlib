#include <global_alignment/two_sequences_linear.h>


namespace biolib {  namespace global_alignment {
 

long two_sequences_simple_gap_linear::score_sequences (size_t s1, size_t s2, similarity const & sim, int gap)
{
  assert(s1 < _align_in.size() && s2 < _align_in.size()); 

  Sf.clear();
  Sf.reserve(_align_in[s2].size() + 1);
  return global_score_forward(
    s1, s2, 0, 0, _align_in[s1].size(), _align_in[s2].size(), sim, gap
  );
}


void two_sequences_simple_gap_linear::align_sequences(size_t s1, size_t s2, similarity const & sim, int gap)
{
  assert(s1 < _align_in.size() && s2 < _align_in.size()); 

  _ga_flag = true;
  
  _ids_out[0] = _ids_in[s1];
  _ids_out[1] = _ids_in[s2];    
  
  Sf.clear();
  Sb.clear();  
  Sf.reserve(_align_in[s2].size() + 1);
  Sb.reserve(_align_in[s2].size() + 1);

  _align_out[0].clear();
  _align_out[1].clear();    
  std::size_t maxsize = _align_in[s1].size() + 
                        _align_in[s2].size();
  _align_out[0].reserve(maxsize);
  _align_out[1].reserve(maxsize);    

  global_block(s1, s2, 0, 0, _align_in[s1].size(), _align_in[s2].size(), sim, gap);
}   


void two_sequences_simple_gap_linear::global_block(size_t s1, size_t s2,
				                   size_t i1, size_t j1,
			                           size_t i2, size_t j2,
			                           similarity const & sim, int gap)
{
  std::size_t i = i2 - i1;   
  std::size_t j = j2 - j1;   

  long score;
  long gap_score;
  long max_score;
  std::size_t max_index;

  if (!j) {
    for (std::size_t l=1; l<=i; l++) {
      _align_out[0].push_back(_align_in[s1][l+i1-1]);
      _align_out[1].push_back('-');
    }
    return;
  }

  if (!i) {
    for (std::size_t l=1; l<=j; l++) {
      _align_out[0].push_back('-');
      _align_out[1].push_back(_align_in[s2][l+j1-1]);
    }
  }
   
  if (i == 1) {    
    gap_score = - gap * (j+1);
      
    for (std::size_t l=1; l<=j; l++) {      
      score = sim.get_similarity(_align_in[s1][i1], _align_in[s2][l+j1-1]);
      if (l == 1) {
        max_score = score;
        max_index = l;
      } else if (score > max_score) { 
        max_score = score;
        max_index = l;
      }
    }
    
    max_score -= gap * (j-1);
        
    if (gap_score < max_score)
    {
      for (std::size_t l=1; l<max_index; l++) 
      {     
        _align_out[0].push_back('-');
        _align_out[1].push_back(_align_in[s2][l+j1-1]); 
      }
      _align_out[0].push_back(_align_in[s1][i1]);
      _align_out[1].push_back(_align_in[s2][max_index+j1-1]); 
      for (std::size_t l=max_index+1; l<=j; l++) 
      {     
        _align_out[0].push_back('-');
        _align_out[1].push_back(_align_in[s2][l+j1-1]); 
      }
      
    } else {  
    
      _align_out[0].push_back(_align_in[s1][i1]);
      _align_out[1].push_back('-'); 

      for (std::size_t l=2; l<=j; l++) {     
        _align_out[0].push_back('-');
        _align_out[1].push_back(_align_in[s2][l+j1-1]); 
      }

    }  
    return;
  } 
    
  std::size_t mid = (i1+i2)/2;

  global_score_forward (s1, s2, i1 , j1, mid, j2, sim, gap);

  global_score_backward(s1, s2, mid, j1,  i2, j2, sim, gap);

  for (std::size_t h=j1; h<=j2; h++) {
     score = Sf[h] + Sb[h];
     if (h == j1) {
       max_score = score;
       max_index = h;
     } else if (score > max_score) { 
       max_score = score;
       max_index = h;
     }
  }
  
  global_block(s1, s2, i1, j1, mid, max_index, sim, gap);

  global_block(s1, s2, mid, max_index, i2, j2, sim, gap);

} 


long two_sequences_simple_gap_linear::global_score_forward(size_t s1, size_t s2,
					 	           size_t i1, size_t j1,
						           size_t i2, size_t j2,
						           similarity const & sim, 
						           int gap)
{
  long s, c;
  
  Sf[j1] = 0;     

  for (std::size_t j=j1+1; j<=j2; j++) 
    Sf[j] = - gap * (j-j1);
  
  for (std::size_t i=i1+1; i<=i2; i++) { 
    s      = Sf[j1];
    Sf[j1] = c = Sf[j1] - gap;
    for (std::size_t j=j1+1; j<=j2; j++) { 
    	c     = std::max(s + sim.get_similarity(_align_in[s1][i-1], _align_in[s2][j-1]), 
	        std::max(Sf[j] - gap , c - gap));
	s = Sf[j];
	Sf[j] = c;
    }
  } 
  return Sf[j2];
}



long two_sequences_simple_gap_linear::global_score_backward(size_t s1, size_t s2,
				  	 	            size_t i1, size_t j1,
						            size_t i2, size_t j2,
 						            similarity const & sim, 
						            int gap)
{
  long s, c;

  Sb[j2] = 0;     

  for (std::size_t j=j2-1;; j--) {
    Sb[j] = - gap * (j2-j);
    if (j == j1) break;
  }
    
  for (std::size_t i=i2-1;; i--) {
    s      = Sb[j2];
    Sb[j2] = c = Sb[j2] - gap;
    for (std::size_t j=j2-1;; j--) { 
    	c     = std::max(s + sim.get_similarity(_align_in[s1][i], _align_in[s2][j]), 
	        std::max(Sb[j] - gap, c - gap));
	s = Sb[j];
	Sb[j] = c;
        if (j == j1) break;
    }
    if (i == i1) break;
  }
  return Sb[j1];
}

/*************************************************************************************************/

long two_sequences_affine_gap_linear::score_sequences (size_t s1, size_t s2, similarity const & sim, int alpha, int beta)
{
  assert(s1 < _align_in.size() && s2 < _align_in.size()); 

  Sf.clear();
  Ff.clear();
  Sf.reserve(_align_in[s2].size() + 1);
  Ff.reserve(_align_in[s2].size() + 1);

  return global_score_forward(
    s1, s2, 0, 0, _align_in[s1].size(), _align_in[s2].size(), sim, alpha, beta, alpha
  );
}


void two_sequences_affine_gap_linear::align_sequences(size_t s1, size_t s2, similarity const & sim, int alpha, int beta)
{
  assert(s1 < _align_in.size() && s2 < _align_in.size()); 

  _ga_flag = true;
  
  _ids_out[0] = _ids_in[s1];
  _ids_out[1] = _ids_in[s2];    
  
  Sf.clear();
  Ff.clear();
  Sb.clear();  
  Fb.clear();  

  Sf.reserve(_align_in[s2].size() + 1);
  Ff.reserve(_align_in[s2].size() + 1);
  Sb.reserve(_align_in[s2].size() + 1);
  Fb.reserve(_align_in[s2].size() + 1);

  _align_out[0].clear();
  _align_out[1].clear();    
  std::size_t maxsize = _align_in[s1].size() + _align_in[s2].size();  
  _align_out[0].reserve(maxsize);
  _align_out[1].reserve(maxsize);    

  global_block(s1, s2, 0, 0, _align_in[s1].size(), _align_in[s2].size(), sim, alpha, beta, alpha, alpha);
}   


void two_sequences_affine_gap_linear::global_block(size_t s1, size_t s2,
				                   size_t i1, size_t j1,
			                           size_t i2, size_t j2,
			                           similarity const & sim, 
						   int alpha, int beta,
						   int gammaf,int gammab)
{
  std::size_t i = i2 - i1;   
  std::size_t j = j2 - j1;   

  long score;
  long gap_score;
  long max_score;
  std::size_t max_index;

  if (!i) {
    for (std::size_t l=1; l<=j; l++) {
      _align_out[0].push_back('-');
      _align_out[1].push_back(_align_in[s2][l+j1-1]);
    }
  }

  if (!j) {
    for (std::size_t l=1; l<=i; l++) {
      _align_out[0].push_back(_align_in[s1][l+i1-1]);
      _align_out[1].push_back('-');
    }
    return;
  }
   
  if (i == 1) {    
    gap_score = - alpha - beta * (j+1) - std::max(gammaf,gammab);
      
    for (std::size_t l=1; l<=j; l++) {      
      score = sim.get_similarity(_align_in[s1][i1], _align_in[s2][l+j1-1]);	             
      if (l > 1) score -= alpha; 
      if (l < j) score -= alpha;
      if (l == 1) {
        max_score = score;
        max_index = l;
      } else if (score > max_score) { 
        max_score = score;
        max_index = l;
      }
    }
    
    max_score -= beta * (j-1);
        
    if (gap_score < max_score)
    { 
      for (std::size_t l=1; l<max_index; l++) 
      {     
        _align_out[0].push_back('-');
        _align_out[1].push_back(_align_in[s2][l+j1-1]); 
      }
      _align_out[0].push_back(_align_in[s1][i1]);
      _align_out[1].push_back(_align_in[s2][max_index+j1-1]); 
      for (std::size_t l=max_index+1; l<=j; l++) 
      {    
        _align_out[0].push_back('-');
        _align_out[1].push_back(_align_in[s2][l+j1-1]); 
      }
      
    } 
    else 
    {  
      _align_out[0].push_back(_align_in[s1][i1]);
      _align_out[1].push_back('-');
      for (std::size_t l=2; l<=j; l++) {     
        _align_out[0].push_back('-');
        _align_out[1].push_back(_align_in[s2][l+j1-1]); 
      }
     
    }  
    return;
  } 
    
  std::size_t mid = (i1+i2)/2;

  global_score_forward (s1, s2, i1 , j1, mid, j2, sim, alpha, beta, gammaf);

  global_score_backward(s1, s2, mid, j1,  i2, j2, sim, alpha, beta, gammab);

  long score1, score2; 
  bool flag;

  for (std::size_t h=j1; h<=j2; h++) {
     score1 = Sf[h] + Sb[h];
     score2 = Ff[h] + Fb[h] + alpha - beta;
     score = std::max(score1,score2);    
     if (h == j1) {
       if(score1 < score2) flag = false;
       else flag = true;
       max_score = score;
       max_index = h;
     } else if (score > max_score) { 
       if(score1 < score2) flag = false;
       else flag = true;
       max_score = score;
       max_index = h;
     }
  }

 if (flag) 
  { 
    global_block(s1, s2, i1, j1, mid, max_index, sim, alpha, beta, gammaf, alpha);

    global_block(s1, s2, mid, max_index, i2, j2, sim, alpha, beta, alpha, gammab);
  } 
  else
  {
    global_block(s1, s2, i1, j1, mid - 1, max_index, sim, alpha, beta, gammaf, 0);

    _align_out[0].push_back(_align_in[s1][mid]);
    _align_out[1].push_back('-'); 
    _align_out[0].push_back(_align_in[s1][mid+1]);
    _align_out[1].push_back('-'); 
    
    global_block(s1, s2, mid + 1, max_index, i2, j2, sim, alpha, beta, 0, gammab);  
  }

} 


long two_sequences_affine_gap_linear::global_score_forward(size_t s1, size_t s2,
						           size_t i1, size_t j1,
						           size_t i2, size_t j2,
						           similarity const & sim, 
						           int alpha, int beta, int gamma)
{
  long s, c, e;
  
  Sf[j1] = 0;     
  Ff[j1] = 0;     

  long gap = - alpha;

  for (std::size_t j=j1+1; j<=j2; j++) {
    Sf[j] = gap; 
    Ff[j] = gap;
    gap -= beta;
  }

  gap = - gamma;
   
  for (std::size_t i=i1+1; i<=i2; i++) 
  {
    s      = Sf[j1];
    e      = gap;
    Sf[j1] = gap;
    Ff[j1] = gap;
    gap -= beta;

    for (std::size_t j=j1+1; j<=j2; j++) { 
        e        = std::max(Sf[j-1] - alpha, e - beta);
	Ff[j]    = std::max(Sf[j] - alpha, Ff[j] - beta);
    	c        = std::max(s + sim.get_similarity(_align_in[s1][i-1], _align_in[s2][j-1]), 
	           std::max(e , Ff[j]));
	s = Sf[j];
	Sf[j] = c;
    }

  }  
  return Sf[j2];
}


long two_sequences_affine_gap_linear::global_score_backward(size_t s1, size_t s2,
				  	 	            size_t i1, size_t j1,
						            size_t i2, size_t j2,
 						            similarity const & sim, 
						            int alpha, int beta, int gamma)
{
  long s, c, e;

  Sb[j2] = 0;     
  Fb[j2] = 0;

  long gap = - alpha;

  for (std::size_t j=j2-1;; j--) 
  {
    Sb[j] = gap;
    Fb[j] = gap;
    gap -= beta;
    if (j == j1) break;
  }
    
  gap = - gamma;

  for (std::size_t i=i2-1;; i--) 
  {
    s      = Sb[j2];
    e      = gap;
    Sb[j2] = gap;
    Fb[j2] = gap;
    gap -= beta;

    for (std::size_t j=j2-1;; j--) { 
      e     = std::max(Sb[j+1] - alpha, e - beta);
      Fb[j] = std::max(Sb[j]   - alpha, Fb[j] - beta);
      c     = std::max(s + sim.get_similarity(_align_in[s1][i], _align_in[s2][j]), 
	      std::max(e , Fb[j]));
      s = Sb[j];
      Sb[j] = c;
      if (j == j1) break;
    }
    if (i == i1) break;
  }
  return Sb[j1];
}



}


} 
