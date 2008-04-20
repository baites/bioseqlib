#ifndef BIOLIB_TREES_BRANCH_PROPORTIONAL_WEIGHTING_INCLUDED
#define BIOLIB_TREES_BRANCH_PROPORTIONAL_WEIGHTING_INCLUDED

#include <trees/tree_weighting_base.h>

namespace biolib { namespace trees {

//! Implementation of branch proportional weighting algorithms.
template<typename I>
class branch_proportional_weighting : public tree_weighting_base<I> {

public:

 //! Void constructor.
 branch_proportional_weighting() : tree_weighting_base<I>() {} 

 //! Constructor by binary tree pointer.
 branch_proportional_weighting(binary_tree<I,float>* tree) : tree_weighting_base<I>(tree) {} 

 //! Evaluates the weights associated to the tree.
 virtual void solve();
 
};

template<typename I>
void branch_proportional_weighting<I>::solve() 
{
  if(tree_weighting_base<I>::_tree == NULL) {
    tree_weighting_base<I>::weights.clear();
    return;
  }
 
  typename binary_tree<I,float>::iterator itr (tree_weighting_base<I>::_tree->begin());  
  typename binary_tree<I,float>::iterator stop(tree_weighting_base<I>::_tree->end());  
  typename binary_tree<I,float>::vertex* node;  

  std::vector<std::size_t> index;
  std::vector<float> weights;
  
  float weight;
  unsigned counter;
        
  for(; itr != stop; ++itr) {
    counter = 1; weight = 0;
    if (itr->leaf()) {
      node = itr.pointer();
      while(node->parent() != NULL) {
        weight += node->edge()/counter;
        node = node->parent();
      	counter++;
      }
      weights.push_back(weight);
      index.push_back(itr->data());
    }
  }
  counter = 1; weight = 0;
  if (itr->leaf()) {
    node = itr.pointer();
    while(node->parent() != NULL) {
      weight += node->edge()/counter;
      node = node->parent();
      counter++;
    }
    weights.push_back(weight);
    index.push_back(itr->data());
  }

  weight = 0;
  for(int i=0; i < weights.size(); i++)
    weight += weights[i];

  tree_weighting_base<I>::weights.resize(weights.size());
  
  for(int i=0; i < weights.size(); i++)
    tree_weighting_base<I>::weights[index[i]] = weights[i] / weight;
} 	

}


}

#endif

