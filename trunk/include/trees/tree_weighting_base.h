#ifndef BIOLIB_TREES_TREE_WEIGHTING_BASE_INCLUDED
#define BIOLIB_TREES_TREE_WEIGHTING_BASE_INCLUDED

#include <biolib.h>
#include <trees/binary_tree.h>


namespace biolib { namespace trees {

//! Based class for weighting trees.
template<typename I>
class tree_weighting_base {

public:

 //! Void constructor.
 tree_weighting_base() { _tree = NULL; }

 //! Copy constructor by binary tree pointer.
 tree_weighting_base(binary_tree<I,float>* tree)
 { _tree = tree; }
 
 //! Sets a tree by a binary tree pointer. 
 void tree(binary_tree<I,float>* tree)
 { _tree = tree; }

 //! Gets the weights.
 std::vector<float> const & get_weights() const
 {
   return weights;
 } 

 //! Calculates the weights.
 virtual void solve() = 0;

protected:

 std::vector<float> weights;

 typename binary_tree<I,float>::binary_tree* _tree;

};


}


}

#endif
