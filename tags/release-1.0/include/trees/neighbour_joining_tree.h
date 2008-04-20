/*********************************************************************************************************

 Library name     :   Biolib++ (Rutines for biological applications)
 File name        :   dtrees.C
 Defined Objects  :   neighbour_joining_tree  

 Programmer       :   Victor E. Bazterra 
 Directed by      :   Julio C. Facelli
                      
 Initialized      :   Wed Nov 12 10:59:47 MST 2003  
 Last updated     :   Wed Nov 12 10:59:47 MST 2003
   
 DESCRIPTION:
    This objects implement all the necesary tools to create binary trees given a set of distances.
 
**********************************************************************************************************/

#ifndef BIOLIB_TREES_NEIGHBOUR_JOINING_TREE_INCLUDED
#define BIOLIB_TREES_NEIGHBOUR_JOINING_TREE_INCLUDED

#include <vector>
#include <iostream>
#include <biolib.h>
#include <trees/binary_tree.h>

namespace biolib { namespace trees {

//! Binary tree construction using neighbour joing tree algorithm.
class neighbour_joining_tree : public binary_tree<int,float> {

// Functions and variables to construct an unrooted tree

 int point;

 int	node;
 int	left_node;
 int	right_node;
 float 	left_dist;
 float 	right_dist;   
 
 matrix distmat;

 std::vector<int> used_nodes;
 std::vector<int> unused_nodes;

 void neighbour_graph_itr();
 void neighbour_graph_remove(int);
 int  neighbour_graph_add(std::vector<float>&);
 void neighbour_graph_nodes(std::vector<int>&); 

 void make_tree();

// Function to calculate a root for the tree

 void step(vertex*, vertex*&, std::vector<bool>&, float&, bool&) const;
 void delta(vertex*, float&, float&, float&, float&, int&, int&, int&) const;
  
public:

  // Constructor by distance matrix.
  neighbour_joining_tree(matrix const &);

  // Creates a root node if the tree.
  void rooted();

};  



}


}
#endif
