/*********************************************************************************************************

 Library name     :   Biolib++ (Rutines for biological applications)
 File name        :   container.h
 Defined Objects  :   container  

 Programmer       :   Victor E. Bazterra 
 Directed by      :   Julio C. Facelli
                      
 Initialized      :   Wed Nov 12 10:59:47 MST 2003  
 Last updated     :   Wed Nov 12 10:59:47 MST 2003
   
 DESCRIPTION:
    Main container class for biolib.
 
**********************************************************************************************************/

#ifndef BIOLIB_CONTAINER_INCLUDED
#define BIOLIB_CONTAINER_INCLUDED

#include <biolib.h>

namespace biolib {

/*! \brief Generic container in biolib. 

  The container class defines a simple generic container based on random access object std::vectors.
  This implementaion was guided by the following goals:

    - Template container style.
    - Fast random access.
    - Dictionary like access in linear time.
    - Light implementation style.

  The random access is the most important feature because most of the biological sequence 
  algorithms rely on it. The basic dictionary like operation are develop only for basic
  manipulation of sequences and symbols.
  
  It is basically to std::vector structures for _id_array (keywords associated to each 
  entry) and _data_array (data given in each entry).  
 */
template<typename idType, typename dataType> 
class container {

public:

  typedef idType id_t;
  typedef dataType data_t;  
  typedef std::vector<idType> id_array_t;
  typedef std::vector<dataType> data_array_t;
     
  //! Void constructor.
  container()
  {
    _name = "";
  }  

  //! Constructor by container name.
  container(std::string const & name)
  {
    set_name(name);
  }

  //! Constructor by an array of id and data.
  container(id_array_t const & id_array, 
	    data_array_t const & data_array)
  {
    _name = "";
    set(id_array,data_array);
  }

  //! Constructor by an array of id, data and a container name.
  container(std::string const & name,  
            id_array_t const & id_array, 
	    data_array_t const & data_array)
  {
    set_name(name);
    set(id_array,data_array);
  }

  //! Copy constructor.
  container(container const & original)
  {
    _name = original._name;
    _id_array = original._id_array;
    _data_array = original._data_array;
  }

  //! Sets a name to the container.
  void set_name(std::string const & name) 
  {
    assert(name.length() > 0);
    _name=name;
  }

  //! Gets the name to the container.
  std::string const & get_id() const 
  {
    return _name;
  }

  //! Sets container from an array of id and data. 
  void set(id_array_t const &, data_array_t const &);

  //! Adds a new entry to the container. 
  bool add(idType const &, dataType const &);

  //! Allocates memory. 
  bool reserve(std::size_t num)
  {
    _id_array.reserve(num);
    _data_array.reserve(num);
  }

  //! Removes an entry by id. 
  bool remove(idType const &);

  //! Returns the number of entries. 
  std::size_t size() const
  {
    return _id_array.size();
  }

  //! Clean the container. 
  void clear() 
  {
    _id_array.clear();
    _data_array.clear();    
  }

  //! Empty verification. 
  bool empty() const
  {
    return _id_array.empty();
  }

  /*! \brief Verify if the content of container is valid. 
   
   This function requieres to be redefine for each new object that use it as base clase. The
   purpose of this function is the verification of the consistence of the information store
   in the container.
   */
  virtual bool check(idType const &, dataType const &) const
  {
    return true;
  }

  //! Verifies an container entry by its id.  
  bool is_id(idType const & id) const
  {
    if (_id_array.empty()) return false;
    return find(id) != _id_array.size();
  }
 
  //! Finds an container entry by its id.  
  std::size_t find(idType const &) const;
     
  //! Access by constant reference to id array.  
  id_array_t const & get_id_array() const
  {
    return _id_array;
  }
  
  //! Access by constant reference to data array.  
  data_array_t const & get_data_array() const
  {
    return _data_array;
  }

protected:

  //! Name of the container.
  std::string _name;

  //! Id array.
  id_array_t _id_array;

  //! Data array.
  data_array_t _data_array;
};


template<typename idType, typename dataType>
bool container<idType, dataType>::add(idType const & id, dataType const & data)
{  
  assert(check(id,data));
  if (!is_id(id))
  {
    _id_array.push_back(id);
    _data_array.push_back(data);
    return true;
  }   
  return false; 
}


template<typename idType, typename dataType>
void container<idType, dataType>::set(
     container<idType, dataType>::id_array_t const & id_array,
     container<idType, dataType>::data_array_t const & data_array)
{  
  assert(id_array.size() == data_array.size());  
  
  for (std::size_t i=0; i<_id_array.size();i++)
    assert(check(id_array[i], data_array[i]));  

  _id_array = id_array;
  _data_array = data_array;
}


template<typename idType, typename dataType>
bool container<idType, dataType>::remove(idType const & id)
{
  typename id_array_t::iterator i;
  typename data_array_t::iterator j=_data_array.begin();
   
  for (i=_id_array.begin(); i!=_id_array.end(); i++)
  {
    if (*i == id)
    {
      _id_array.erase(i);
      _data_array.erase(j);
      return true;
    }
    j++; 
  }
  return false;
}  


template<typename idType, typename dataType>
std::size_t container<idType, dataType>::find(idType const & id) const
{
  for (std::size_t i=0; i<_id_array.size();i++)
    if (_id_array[i] == id) return i;
  return _id_array.size();
}


}

#endif
