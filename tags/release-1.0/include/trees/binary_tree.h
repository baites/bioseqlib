/*********************************************************************************************************

 Library name     :   Biolib++ (Rutines for biological applications)
 File name        :   binary_tree.h
 Defined Objects  :   binary_tree, binary_tree::vertex, binary_tree::iterator  

 Programmer       :   Victor E. Bazterra 
 Directed by      :   Julio C. Facelli
                      
 Initialized      :   Wed Nov 12 10:59:47 MST 2003  
 Last updated     :   Fri Nov 11 02:41:57 MST 2005
   
 DESCRIPTION:
    This objects implement all the necesary tools to create and manipulate binary trees.
 
**********************************************************************************************************/

#ifndef BIOLIB_TREES_BINARY_TREE_INCLUDED
#define BIOLIB_TREES_BINARY_TREE_INCLUDED

#include <biolib.h>

#define __need_NULL
#include <stddef.h>


namespace biolib { 
//! Contains all the algorithms related to trees.
namespace trees {

/*! \brief Binary tree representation and manipulation. 

  The binary_tree class defines a group of three classes to represent and manipulate binary trees. 
  The implementation was guided by the following goals:

    - General and simple implementation.
    - Avoidance to create cycles.
    - Template container style.
    - Access through iterators.
    - Exception handler.

  A tree is a cycle-free connected graph. As any graph, it has a collection of vertices and edges 
  that establish relationships between pairs of vertices. A binary tree has all vertices of degree 
  one or of degree three. Those vertices of degree one are called the terminal or labeled vertices 
  that for short we will call them leaves. 

  Frequently a tree is represented with a hierarchy up - down where for given vertex it is possible 
  to classify its neighbor vertices as its parent and/or left and right children see figure. In theses 
  cases, the node that does not have a parent is called the root of the tree.
  
  \image html biolib-binary_trees_fig1.jpg "Rooted and unrooted trees"
  
  In the case that a hierarchy is not established the resulting tree is known as an unrooted tree. 
  It is always possible to "insert" a root to convert an unrooted tree into a rooted one. Moreover, 
  it  will be shown that an unrooted tree can be represented as an rooted tree.
  
  It would be convenient to extend the definition of binary tree allowing that its vertices have a 
  degree of 0 to 3. These more general binary trees are called incomplete trees. The case of a 
  vertex of degree zero is a tree with only one vertex.

  The binary_tree class implements the followed functional qualities:

    - Creation and edition of a tree.
    - Access to vertices and their upper edges objects.
    - Memory allocation and deallocation of the vertices.

  This class is basically implemented by two vertex pointers _pointer that points to the root node 
  and _current the points to the actual vertex. A tree is a empty tree when it is not controling 
  (pointing) any vertices. This case is implemented saying that both pointers are null.
  
  All the vertices of a tree are created and destroyed through the binary_tree interface, keeping all 
  the vertices interconnected through the chain of pointers to finish in _pointer. Therefore, the 
  deallocation of all the vertices can be done controlling the deallocation of _pointer. This can be 
  triggered by clear() member function or by an exception that could be thrown by the program 
  deallocating automatically all the vertices. This property makes this implementation free of 
  memory leaks.
  
  In this implementation, the edition and creation of the trees is based on subtrees composition. The
  left or right side of a tree can be extended by pointing to a new tree. This can be done using left()
  and right() member function under the only condition the composition is given between different trees. 
  This last condition forbids the creation of graph with cycles that are not by definition a tree.
 */
template<typename I, typename J>
class binary_tree {

public:

 class iterator;

 /*! \brief Vertex for a binary tree. 

  The nested class vertex implements a vertex for a binary tree. It is basically three pointers _parent,
  _left, _right and two objects _data and _edge. The pointers establish the relationship between vertex 
  objects playing as edges between vertices. It is often to associate to each vertex and edge with a 
  particular information. This information can be put in _data and _edge objects, where _data and _edge 
  are the information associate to the vertex and its upper edge respectively, see figure.
  
  \image html biolib-binary_trees_fig2.jpg "Structure of the vertex class"
  
  The template parameter I defines the object type associated to the vertex and J the object type 
  associated to the edge. The vertex objects are meant to not be used directly. All the client programmer 
  interaction with the vertex objects will be carried by binary_tree objects or iterators in the way that
  all vertex are created and destroyed by trees. This will be allowed that if any exception ocurrs, the 
  binary_tree objects will deallocate automatically all the vertices from memory.
 */ 
 class vertex {

   friend class iterator;
   friend class binary_tree<I,J>;

   //! Object type I associated with the vertex.
   I       _data;

   //! Object type J associated with the upper edge of vertex.
   J	   _edge;

   //! Pointer to the left vertex.
   vertex* _left;

   //! Pointer to the right vertex.
   vertex* _right;

   //! Pointer to the parent of the vertex.
   vertex* _parent;

// Modification member functions

   //! Sets the left pointer of the vertex.
   void left  (vertex*);

   //! Sets the right pointer of the vertex.
   void right (vertex*);

   /*! \brief Transforms parent to child vertex relation.

    It also changes the edge values to keep the relationship between verteces. 
   */
   void transf_parent_child(const J&);

   /*! \brief Constructor by data and edge information.
   
    This is the constructor of the class that initialice the values of the objects related to the vertex 
    and its upper edge. By default, the assigned value is given by the default constructor of the objects
    I and J.

    \note All the function related to create or modify the relationship between vertices are private 
    making imposible for the user declares or modifies the tree topology directly.
   */
   vertex(const I& data = I(), 
          const J& edge = J())
   {
     _data    = data;
     _edge    = edge;
     _left    = NULL;
     _right   = NULL;
     _parent  = NULL;
   }

   //! Copy constructor.
   vertex(const vertex&);

 public:
   
   /*! \brief Vertex destructor.
   
    It destroys the actual vertex together with all the vertices that are interconnected to its left and 
    right pointers.
   */
   ~vertex();
     
// Modification member functions

   //! Removes all the vertices connected to the left.
   void remove_left();

   //! Removes all the vertices connected to the right.
   void remove_right(); 

   //! Sets the object associated to the vertex.
   void data (const I& data) { _data = data; }

   //! Sets the object associated to the upper edge.
   void edge (const J& edge) { _edge = edge; }
 
// Constant member functions

   //! Gets the object associated to the vertex.
   const I& data()    const { return _data;   }

   //! Gets the object associated to the upper edge.
   const J& edge()    const { return _edge;   }

   //! Gets the left pointer of the vertex.
   vertex*  left ()   const { return _left;   }

   //! Gets the right pointer of the vertex.
   vertex*  right ()  const { return _right;  }

   //! Gets the parent of the vertex.
   vertex*  parent () const { return _parent; }
   
   /*! \brief Verified if the vertex is a leaf. 
   
    It returns true when the vertex has the left and right pointers NULL.
    */
   bool leaf() const 
   { return (_left == NULL) && (_right == NULL); }

   /*! \brief Verified if the vertex is a terminal. 

    It returns true when the vertex is degree one.
    */
   bool terminal() const 
   { return (_left == NULL   && _right  == NULL) || 
            (_right == NULL  && _parent == NULL) ||
            (_parent == NULL && _left   == NULL); }
	          
 };

protected:
 
 //! Points to the root of the tree.
 vertex* _pointer;

 //! Points to the actual vertex.
 vertex* _current;
 
public:

 //! Void constructor.
 binary_tree() { 
   _pointer = NULL; _current = NULL;
 }
 
 //! Non-empty tree constructor.
 binary_tree(const I& data, const J& edge = J()) {
   _pointer = _current = new vertex(data, edge);
 }

 //! Copy constructor by a vertex pointer.
 binary_tree(const vertex* pvertex) {
   _pointer = _current = new vertex(*pvertex);
 }
    
 //! Copy constructor.
 binary_tree(const binary_tree& tree); 

 //! Binary_tree destructor.
 ~binary_tree() {
   clear();
 }

 //! Deallocates all the vertices of a tree makign an empty one.
 void clear();

 //! Creates the first vertex on an empty tree.
 void initialize(const I& data = I(), const J& edge = J()) {
   assert(empty()); _pointer = _current = new vertex(data, edge);
 }

// Constant member functions

 /*! \brief Sets the current pointer into the left.
  
  If left vertex pointer of _current is not null then sets the current pointer 
  into its left vertex and returns true. Otherwise _current is not modified and 
  the function returns false.
  */
 bool left ();

 /*! \brief Sets the current pointer into the right.
  
  If right vertex pointer of _current is not null then sets the current pointer 
  into its right vertex and returns true. Otherwise _current is not modified and 
  the function returns false.
  */
 bool right ();

 /*! \brief Sets the current pointer into the parent.
  
  If parent vertex pointer of _current is not null then sets the current pointer 
  into its parent vertex and returns true. Otherwise _current is not modified and 
  the function returns false.
  */
 bool parent ();

 //! Sets current pointer into next vertex in preorder.
 void preorder();
 
 //! Sets current pointer into next vertex in postorder.
 void postorder();
 
 //! Sets current pointer into root.
 void reset() { _current = _pointer; }

 //! Returns true if the tree is empty.
 bool empty() const { return _pointer == NULL; }

 //! Gets data for the current vertex.
 const I& data() const 
 { assert (!empty()); return _current->data(); }

 //! Gets edge for the current vertex.
 const J& edge() const
 { assert (!empty()); return _current->edge(); }

 /*! \brief Verified if the current vertex is a leaf. 
   
  It returns true when the vertex has the left and right pointers NULL.
  */
 bool leaf() const 
 { assert (!empty()); return _current->leaf(); }

 /*! \brief Verified if the current vertex is a terminal. 

  It returns true when the vertex is degree one.
  */
 bool terminal() const
 { assert (!empty()); return _current->terminal(); }
   
// Modification member functions

 /*! \brief Creates a new root vertex. 

  It creates a root nodes cutting the upper edge of the current vertex.
  */
 void root(const I&);

 /*! \brief Makes cycle permutation to the left.

  This is an operation that makes a cycle permutation to the left of the neighbor vertices.

  \image html biolib-binary_trees_fig3.jpg "Left and right rotation of a tree"
 */
 void left_rotation();

 /*! \brief Makes cycle permutation to the right.

  This is an operation that makes a cycle permutation to the right of the neighbor vertices.

  \image html biolib-binary_trees_fig3.jpg "Left and right rotation of a tree"
 */
 void right_rotation();

 //! Set the left pointer of the current vertex.
 void left  (binary_tree<I,J>&);

 //! Set the right pointer of the current vertex.
 void right (binary_tree<I,J>&);

 /*! \brief Releases the vertices of a given tree.
 
  It releases the vertices of a given tree an put them under the control of the calling tree. 
  If the callign tree is not empty it will deallocate automatically all its vertices to take 
  cantrol over the new ones.
  */ 
 void release(binary_tree<I,J> &tree);

 //! Sets data for the current vertex.
 void data  (const I& data) 
 { assert (!empty()); _current->data(data); }

 //! Sets edge for the current vertex.
 void edge  (const J& edge) 
 { assert (!empty()); _current->edge(edge); }
  
 /*! \brief An iterator for a binary tree. 

  This class implements a pointer like object to the vertices of a tree. The goal is to access all the 
  nodes of a tree in preorder. This access has the constrains given by the vertex interface allowing the 
  modification of the objects related to the nodes and its upper edge but forbidden to modify the left 
  and right pointers.

  Another restriction is that an iterator can be created only by another iterator. The initial iterator 
  is given by the begin() member function of the tree that contains the group of vertices to be pointed.
  */ 
 class iterator {

   friend class binary_tree<I,J>;
   
   //! Pointer to the actual vertex.
   vertex* _current;
 
 public:
 
   //! Void constructor.
   iterator() 
   { _current = NULL; } 

   //! Constructor by a vertex pointer.
   iterator(vertex* node) 
   { _current = node; } 

   //! Copy constructor.  
   iterator(const iterator& node) 
   { operator=(node); } 

   //! Increment iterator in prorder.
   iterator& operator++();

   //! Decrement iterator in prorder.
   iterator& operator--(); 

   //! Increment iterator in prorder returning the caller vertex.
   iterator operator++(int )
   { iterator aux(*this); operator++(); return aux;}

   //! Decrement iterator in prorder returning the caller vertex.
   iterator operator--(int )
   { iterator aux(*this); operator--(); return aux;}

   //! Returns a vertex pointer to the vertex that it is been pointed.
   vertex* pointer()    const { return  _current; }
   //! Returns a vertex pointer to the vertex that it is been pointed.
   vertex* operator->() const { return  _current; }
   //! Returns the vertex that it is been pointed.
   vertex& operator*()  const { return *_current; }
 
   //! Define the asignation operator.
   void operator= (const iterator& i)
   { _current = i._current; }

   //! Define the equal comparator operator.
   bool operator== (const iterator& i) const
   { return _current == i._current; }

   //! Define the no-equal comparator operator.
   bool operator!= (const iterator& i) const
   { return _current != i._current; }

 };
 
 //! Returnrs the number of vertices of the tree.    
 int size() const;
 
 //! Basic printing of the tree structure.
 void print() const;
  
 //! Returns an iterator to the right most vertex of the tree.
 iterator end() const;

 //! Returns an iterator pointing to the root of the tree.
 iterator begin() const { return iterator(_pointer); }
 
};


//**********************************************************
//                 binary_tree functions
//**********************************************************


template<typename I, typename J>
binary_tree<I,J>::binary_tree(const binary_tree<I,J>& orign)
{ 
  if (orign._pointer == NULL) {
    _pointer = NULL;
    _current = NULL;
    return;
  }

  _pointer = _current = new vertex(*orign._pointer);
}    

template<typename I, typename J>
void binary_tree<I,J>::clear()
{ 
  if (_pointer != NULL) {
    delete _pointer;
    _pointer = NULL;
  }
}

template<typename I, typename J>
bool binary_tree<I,J>::left() 
{
  if (_current == NULL)
    return false;
  else if (_current->_left == NULL)
    return false;
  
  _current = _current->_left;    

  return true;
}

template<typename I, typename J>
bool binary_tree<I,J>::right() 
{
  if (_current == NULL)
    return false;
  else if (_current->_right == NULL)
    return false;
  
  _current = _current->_right;    

  return true;
}

template<typename I, typename J>
bool binary_tree<I,J>::parent() 
{
  if (_current == NULL)
    return false;
  else if (_current->_parent == NULL)
    return false;
  
  _current = _current->_parent;    

  return true;
}

template<typename I, typename J>
void binary_tree<I,J>::preorder() 
{
  if (_current == NULL) return;

  if (_current->_left != NULL)
    _current = _current->_left;
  else if (_current->_right != NULL)
    _current = _current->_right;
  else {
    vertex* old;
    do {       
      old = _current;
      if (_current->_parent == NULL) return;  
      _current = _current->_parent;
    } while( _current->_right == NULL || _current->_right == old );
    _current = _current->_right;
  }
  return;
}

template<typename I, typename J>
void binary_tree<I,J>::postorder() 
{
  if (_current == NULL) return;


  if (_current->_right != NULL)
    _current = _current->_right;
  else if (_current->_left != NULL)
    _current = _current->_left;
  else {
    vertex* old;
    do {       
      old = _current;
      if (_current->_parent == NULL) return;  
      _current = _current->_parent;
    } while( _current->_left == NULL || _current->_left == old );
    _current = _current->_left;
  }
  return;
}

template<typename I, typename J>
void binary_tree<I,J>::left (binary_tree<I,J>& tree)
{
  assert(_pointer != tree._pointer);

  if(_current->_left != NULL) {
    delete _current->_left;
    _current->_left = NULL;
  }
  
  _current->left(tree._pointer); 
  
  tree._pointer = NULL;
  tree._current = NULL;
}

template<typename I, typename J>
void binary_tree<I,J>::right (binary_tree<I,J>& tree)
{
  assert(_pointer != tree._pointer);

  if(_current->_right != NULL) {
    delete _current->_right;
    _current->_right = NULL;
  }
    
  _current->right(tree._pointer); 

  tree._pointer = NULL;
  tree._current = NULL;
}

template<typename I, typename J>
void binary_tree<I,J>::release(binary_tree<I,J> &tree)
{
  assert(_pointer != tree._pointer);

  if(!empty()) clear();
 
  _pointer = _current = tree._pointer;
  
  tree._pointer = NULL;   
  tree._current = NULL;
}   

template<typename I, typename J>
void binary_tree<I,J>::left_rotation()
{
  assert(!empty());

  vertex* node = _current;

// Go to the most right node
  
  while(node->_right != NULL)
    node = node->_right;

// Rotation between pointers

  node->transf_parent_child();
  node->right(node->_left);
  node->left(node->_parent);    
  node->_parent = NULL;
  _pointer = node;
}  

template<typename I, typename J>
void binary_tree<I,J>::right_rotation()
{
  assert(!empty());

  vertex* node = _current;

// Go to the most left node
  
  while(node->_left != NULL)
    node = node->_left;

// Rotation between pointers
    
  node->transf_parent_child();
  node->left(node->_right);
  node->right(node->_parent);    
  node->_parent = NULL;
  _pointer = node;
}  

template<typename I, typename J>
void binary_tree<I,J>::root(const I& data = I())
{
  assert(!empty());
  
  if (_current->_parent == NULL) return;

  vertex* root = new vertex(data);
  if (_current == _current->_parent->_left) {    
    _current->transf_parent_child();    
    root->right(_current->_parent);    
    root->left(_current);
  } else {
    _current->transf_parent_child();    
    root->left(_current->_parent);    
    root->right(_current);
  }    	
  
  _pointer = root;
}

template<typename I, typename J>
int binary_tree<I,J>::size() const
{
  int mysize = 1;
  
  if (empty()) return 0;
 
  iterator itr(begin());
  iterator ende(end());
 
  for (; itr != ende; ++itr)
    mysize++;
  
  return mysize;
}

template<typename I, typename J>
typename binary_tree<I,J>::iterator binary_tree<I,J>::end() const
{ 
  assert(!empty());

  vertex* node = _pointer;

  while (!node->leaf()) {
    if(node->_right != NULL)	
      node = node->_right;
    else
      node = node->_left;
  }
  
  return iterator(node);
}

template<typename I, typename J>
void binary_tree<I,J>::print() const
{
   iterator itr(begin());

   std::cout << "Tree size is " << size() << std::endl; 
   std::cout << "Tree in preorder :" << std::endl;
   for (itr = begin(); itr != end(); ++itr) {
     std::cout <<   "Node's data : " << itr->data() << std::endl;
     if (itr->parent() != NULL)
       std::cout << "Node's edge : " << itr->edge() << std::endl;
     else
       std::cout << "Node's edge : NULL" << std::endl;
     if (itr->left() != NULL)
       std::cout << "Left   data : " << itr->left()->data() << std::endl;
     else
       std::cout << "Left   data : NULL" << std::endl;
     if (itr->right() != NULL)
       std::cout << "Right  data : " << itr->right()->data() << std::endl;
     else
       std::cout << "Right  data : NULL" << std::endl;
     std::cout << std::endl;
   }
   std::cout <<   "Node's data : " << itr->data() << std::endl;
   if (itr->parent() != NULL)
     std::cout << "Node's edge : " << itr->edge() << std::endl;
   else
     std::cout << "Node's edge : NULL" << std::endl;
   if (itr->left() != NULL)
     std::cout << "Left   data : " << itr->left()->data() << std::endl;
   else
     std::cout << "Left   data : NULL" << std::endl;
   if (itr->right() != NULL)
     std::cout << "Right  data : " << itr->right()->data() << std::endl;
   else
     std::cout << "Right  data : NULL" << std::endl;
     
   std::cout << std::endl;
   std::cout << std::endl;
}


//**********************************************************
//             binary_tree::vertex functions
//**********************************************************


template<typename I, typename J>
binary_tree<I,J>::vertex::vertex(const vertex& orign)
{ 
  _parent = NULL;
  _data   = orign._data;
  _edge   = orign._edge;
    
  if (orign._left != NULL) {
    _left = new vertex(*orign._left);
    _left->_parent = this;
  } else _left = NULL;
  
  if (orign._right != NULL) {
    _right = new vertex(*orign._right);
    _right->_parent = this;
  } else _right = NULL;
}

template<typename I, typename J>
binary_tree<I,J>::vertex::~vertex()
{
  if (this != NULL) {
    delete _left;
    _left = NULL;
    delete _right;
    _right = NULL;
  }
}

template<typename I, typename J>
void binary_tree<I,J>::vertex::remove_left ()
{
  if (_left != NULL) {
    delete _left;
    _left = NULL;
  }
}

template<typename I, typename J>
void binary_tree<I,J>::vertex::remove_right ()
{
  if (_right != NULL) {
    delete _right;
    _right = NULL;
  }
}

template<typename I, typename J>
void binary_tree<I,J>::vertex::left (vertex* node)
{ 
  _left = node; 
  if (node != NULL) 
    node->_parent = this;
}

template<typename I, typename J>
void binary_tree<I,J>::vertex::right (vertex* node)
{ 
  _right = node; 
  if (node != NULL) 
    node->_parent = this;
}
    
template<typename I, typename J>
void binary_tree<I,J>::vertex::transf_parent_child(const J& last= J())
{
// Change the edge values to keep 
// the relation ship between nodes 

  J itr;
  
  itr   = _edge; 
  _edge = last;      

  if (_parent != NULL) {
    _parent->transf_parent_child(itr);
    if (this == _parent->_left) { 
      _parent->left(_parent->_right);
      _parent->right(_parent->_parent);
    } else {
      _parent->right(_parent->_left);
      _parent->left(_parent->_parent);
    }      
  } 
}


//**********************************************************
//             binary_tree::iterator functions
//**********************************************************


template<typename I, typename J>
typename binary_tree<I,J>::iterator& binary_tree<I,J>::iterator::operator++() 
{
  if (_current == NULL)
    return *this;

  if (_current->_left != NULL)
    _current = _current->_left;
  else if (_current->_right != NULL)
    _current = _current->_right;
  else {
    vertex* old;
    do {       
      old = _current;
      if (_current->_parent == NULL)
        return *this;  
      _current = _current->_parent;
 
    } while( _current->_right == NULL || _current->_right == old );
    _current = _current->_right;
  }
  return *this;
}

template<typename I, typename J>
typename binary_tree<I,J>::iterator& binary_tree<I,J>::iterator::operator--() 
{
  if (_current == NULL)
    return *this;

  if (_current->_parent != NULL)  {
    if (_current == _current->_parent->_left || 
      _current->_parent->_left == NULL)
      _current = _current->_parent;
    else {
      _current = _current->_parent->_left;             
      while (!_current->leaf()) {
        if (_current->_right != NULL) 
          _current = _current->_right;
       else  
          _current = _current->_left;
      }    
    }
  } else {
    while (!_current->leaf()) {
      if (_current->_right != NULL) 
        _current = _current->_right;
      else  
        _current = _current->_left;
    }    
  }
  return *this;
}


}


}

#endif
