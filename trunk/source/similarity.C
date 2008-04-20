#include <similarity.h>


namespace biolib {


similarity::similarity(similarity const & original)
{
  _name = original._name;
  _score_matrix = original._score_matrix;

  for(int i=0; i<256; i++)
    _translate[i] = original._translate[i];
}


void similarity::set_similarity (alphabet const & a, short* const m) 
{
  _id_array = a.get_id_array();
 
  for (std::size_t i=0; i<_id_array.size(); i++) 
    _translate[static_cast<int> (_id_array[i])] = i;

  std::size_t dim = (_id_array.size()/2)*(_id_array.size()+1);
  for (std::size_t i=0; i<dim; i++) 
    _score_matrix.push_back(m[i]);  
}

    
short similarity::get_similarity (char a, char b) const 
{
  std::size_t maxcode = std::max(_translate[static_cast<int>(a)], _translate[static_cast<int>(b)]);
  std::size_t mincode = std::min(_translate[static_cast<int>(a)], _translate[static_cast<int>(b)]);
  std::size_t code = maxcode * (maxcode+1)/2 + mincode;
 
  assert(code < _score_matrix.size());

  return _score_matrix[code];
}


std::ostream& operator<< (std::ostream& o, similarity const & a)
{ 
  std::size_t k = 0;

  o << "'" << a._name << "' score matrix :\n";

  for (std::size_t i=0; i<a._id_array.size(); i++) {
    o << a._id_array[i] << ") "; 
    for (std::size_t j=0; j<i+1; j++)  
      o << std::setw(4) << a._score_matrix[k++];
    o << std::endl;
  }

  return o;
}


}
