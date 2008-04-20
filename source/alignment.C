#include <alignment.h>


namespace biolib {


bool alignment::check(id_t const &, data_t const & data) const
{  
  if(data.length() <= 0 || _alphabet_flag == false)
    return false;

  for (std::size_t i=0; i<data.length(); i++) 
    if(_alphabet.is_id(data[i])==false)
      return false;

  return true;
}


std::ostream& operator<< (std::ostream & o, alignment const & a)
{  
  std::size_t k = 0;
  std::size_t nseqs = 0;
  std::size_t title = 20;
  std::size_t nraw = 80; 
  std::size_t scrlen = nraw - title;

  std::vector<char> buffer(scrlen);

  std::size_t maxsize = 0;

  for(std::size_t i=0; i<a.size(); i++)
    maxsize = std::max(maxsize, a._data_array[i].size());

  o << "Alignment " << "'" << a._name << "'" << " with ";
  o << a.size() << " sequences : " << std::endl;

  while(nseqs < a.size()) 
  {
    o << std::endl; 
    for (std::size_t i=0; i<a.size(); i++) 
    {      
	if(i)
	{
  	  o << std::setw(title) << " ";  
          for (std::size_t j=k*scrlen; j<scrlen + k*scrlen; j++) 
	  {	              
            if (k*scrlen >= a._data_array[i].size() || 
	        j >= a._data_array[i].size()) break;

            if (buffer[j-k*scrlen] == '-' || a._data_array[i][j] == '-')
	      o << ' ';
	    else if (a._data_array[i][j] == buffer[j-k*scrlen]) 
	      o << '!'; 	         
	    else 
	      o << '.';         
	  }
          o << std::endl; 
	}
	else
	{
  	  o << std::setw(title-2) << k*scrlen << "  " ;  
          for (std::size_t j=k*scrlen; j<scrlen + k*scrlen; j++) 
	  {	
	    if (j >= maxsize)
	      break;              
	    if ((j+1) % 10 == 0)
	      o << ':';	
	    else if ((j+1) % 5 == 0)
	      o << '\'';		    
	    else
	      o << ' ';		      
	  }
          o << std::endl; 	
	}
	
	buffer = std::vector<char>(scrlen,'-');	   
	o << std::setw(title-2) << a._id_array[i] << ": ";  
        for (std::size_t j=k*scrlen; j<scrlen + k*scrlen; j++) 
	{	              
           if (k*scrlen >= a._data_array[i].size()) 
       	     break;

	   o << a._data_array[i][j];
           buffer[j-k*scrlen] = a._data_array[i][j];
	   
	   if (j >= a._data_array[i].size() - 1) 
	   { 
	     nseqs++;
       	     break;
	   }
	}
        o << std::endl; 
    }
    k++;
    o << std::endl << std::endl; 
  }
  o << std::endl << std::endl;

  return o;
}

}
