#include <trees/neighbour_joining_tree.h>


namespace biolib { namespace trees {


neighbour_joining_tree::neighbour_joining_tree(matrix const & dists)
{
    point = 0;
    distmat   = dists;
    
    used_nodes.resize(dists.size());
    unused_nodes.resize(dists.size());

    for(std::size_t i=0; i<dists.size(); i++)
      used_nodes[i] = i;	

    make_tree();
}
        
void neighbour_joining_tree::neighbour_graph_remove(int node) {
    used_nodes[node] = -1;	
    unused_nodes[point] = node;
    point++;
}

void neighbour_joining_tree::neighbour_graph_nodes(std::vector<int>& nodes) 
{
    int idx;
    nodes.clear();
    for (int i=0; i<used_nodes.size(); i++) {
      idx = used_nodes[i]; 
      if(idx != -1) 
         nodes.push_back(idx);
    }
}

int neighbour_joining_tree::neighbour_graph_add(std::vector<float>& dists)
{
    int i, len, node; 
    
    point--; 

    if (point >= 0)
      node = unused_nodes[point];
    else 
      node = used_nodes.size();
     
    used_nodes[node] = node;	 

    std::vector<int> nodes;

    neighbour_graph_nodes(nodes);
  
    len = nodes.size();
  
    for(i=0; i<len; i++) 
      if (node != nodes[i])
        distmat[node][nodes[i]] = distmat[nodes[i]][node] = dists[nodes[i]];
      else
        distmat[node][node] = 0;
   
    return node;
}

void neighbour_joining_tree::neighbour_graph_itr() {

    int flag = false;
    int nodei, nodej;
    float dist, D, minD, sumd;

    std::vector<float>        r;     
    std::vector<int> 	  nodes;

    std::vector<float> 	     dk(distmat.size());     

    neighbour_graph_nodes(nodes);
 
    int len = nodes.size();
                
    r.resize(len);     

    for (int i=0; i<len; i++) {
      sumd = 0;
      for (int j=0; j<len; j++) 
        sumd += distmat[nodes[i]][nodes[j]];
      r[i] = sumd/static_cast<float>(len - 2);
    }  

    if(len > 2) {
      for (int i=0; i<len; i++)
        for (int j=i+1; j<len; j++) {
          dist = distmat[nodes[i]][nodes[j]];
    	  D = dist - r[i] - r[j];
      	  if (dist > 0 && flag == false) {
      	    minD = D;
      	    nodei = i; nodej = j;
      	    flag = true;	
      	  } else if (dist > 0 && minD > D) {
       	    minD = D;
     	    nodei = i; nodej = j;
          }
        }
      
      for (int i=0; i<len; i++) 
        dk[nodes[i]] = 0.5*(distmat[nodes[nodei]][nodes[i]] + 
		            distmat[nodes[nodej]][nodes[i]] -
	                    distmat[nodes[nodei]][nodes[nodej]]);

      left_node   = nodes[nodei];
      right_node  = nodes[nodej];
           
      left_dist   = 0.5*(distmat[nodes[nodei]][nodes[nodej]] + r[nodei] - r[nodej]);
      right_dist  = 0.5*(distmat[nodes[nodei]][nodes[nodej]] + r[nodej] - r[nodei]);

      neighbour_graph_remove(nodes[nodei]);
      neighbour_graph_remove(nodes[nodej]);

    } else if (len == 2) { 

      left_node   = nodes[0];
      right_node  = nodes[1];
      left_dist   = distmat[nodes[0]][nodes[1]];
      right_dist  = left_dist;

      neighbour_graph_remove(nodes[0]);
      neighbour_graph_remove(nodes[1]);
      
    }
    
    node = neighbour_graph_add(dk);
}

void neighbour_joining_tree::make_tree() 
{ 
    int len = distmat.size();
    binary_tree<int,float> tmptree;
    std::vector<binary_tree<int,float> > nodes;    

// allocate leaf nodes

    for (int i=0; i<len; i++)
      nodes.push_back(binary_tree<int,float>(i));
    
    int count = len;
    
// calculate the tree until the last two subtrees
     
    for (int i=0; i<len-2; i++) {
      neighbour_graph_itr();
      nodes[left_node ].edge(left_dist);
      nodes[right_node].edge(right_dist);
      tmptree.initialize(count);
      tmptree.left (nodes[left_node] );
      tmptree.right(nodes[right_node]);
      nodes[node].release(tmptree);
      count++;
    } 
 
    neighbour_graph_itr();

// rotation to the right of the left subtree

    nodes[left_node].right_rotation();

// union of the rotate left subtree and right subtree

    nodes[right_node].edge(right_dist);
    nodes[left_node].right(nodes[right_node]);

// release rotated left subtree        
    release(nodes[left_node]);
}

void neighbour_joining_tree::step(vertex* refnode, vertex*& itrnode, 
                                  std::vector<bool>& pass, float& dist, bool& flag) const
{
  float tmp;
  
  if (itrnode == NULL) return;
  
  if (itrnode == refnode) {
    dist = 0;
    flag = false;
    pass.clear();
    pass.resize(size(),false);
  }
        
  pass[itrnode->data()] = true;                         

  if (itrnode->left() != NULL) 
  {
    itrnode = itrnode->left();

    if(!pass[itrnode->data()]) {
      dist += itrnode->edge();  
      pass[itrnode->data()] = true;                         
    } else 
      dist -= itrnode->edge();          
  } 
  else if (itrnode->right() != NULL) 
  {
    itrnode = itrnode->right();
 
    if(!pass[itrnode->data()]) 
    {        
       dist += itrnode->edge();          
       pass[itrnode->data()] = true;           
    } 
    else 
       dist -= itrnode->edge();          
  } 
  else 
  {
    vertex* old;

    do {       
         	
      old = itrnode;
            
      if (itrnode->parent() == NULL) 
        return;
      
      tmp = itrnode->edge();          
      itrnode = itrnode->parent();

      if(!pass[itrnode->data()]) {        
         dist += tmp;          
         pass[itrnode->data()] = true;           
      } else 
         dist -= tmp;          

      if ( refnode->parent() == NULL)
        flag = false;
      else if (itrnode == refnode->parent())
        flag = true;

    } while( itrnode->right() == NULL || itrnode->right() == old );

    itrnode = itrnode->right();
 
    if(!pass[itrnode->data()]) {        
      dist += itrnode->edge();          
      pass[itrnode->data()] = true;           
    } else 
      dist -= itrnode->edge();                 
  }
  return;
}


void neighbour_joining_tree::delta(vertex* refnode , 
				   float& umaxdist , float& dmaxdist,
				   float& usumdist , float& dsumdist,
				   int&   uleaves  , int&   dleaves, int& nnodes   ) const
{
  float dist = 0;

  vertex* itrnode = refnode;

  nnodes   = 0;
  uleaves  = 0;
  dleaves  = 0;
  usumdist = 0;
  dsumdist = 0;
  umaxdist = 0;
  dmaxdist = 0;

  if(itrnode == NULL) return;

  bool up;

  float udelta = 0;
  float ddelta = 0;

  std::vector<bool> pass;

  step(refnode, itrnode, pass, dist, up);
  
  while(refnode != itrnode) {
    nnodes++;
    if (up) {
      if (itrnode->terminal()) {
        uleaves++;
        usumdist += dist;
        if (dist > umaxdist) umaxdist = dist;
      }
    } else {
      if (itrnode->terminal()) {
        dleaves++;
        dsumdist += dist;
        if (dist > dmaxdist) dmaxdist = dist;
      }
    }
    step(refnode, itrnode, pass, dist, up);
  }
} 


void neighbour_joining_tree::rooted()
{

  if (empty()) return;
  
  bool     first = true;
  float    minxroot = 0;
  float    minmaxdist = 0;
  
  float usum, dsum;
  float umax, dmax; 
  int   uleaves, dleaves, nnodes;

  float xroot, dist, maxdist;

  vertex* minnode;

  _current = _pointer;
  
  do {  
    
    if (_current->parent() == NULL) {
      preorder();
      continue;
    }
     
    delta(_current, umax, dmax, usum, dsum, uleaves, dleaves, nnodes);

    float dA, dB;

    if(dleaves == 0)
      dA = usum/uleaves;
    else
      dA = usum/uleaves - dsum/dleaves;

    if(dleaves == 0)
    {
      dB = usum/uleaves - 2 * _current->edge();
      dsum = _current->edge();
      dleaves = 1;
    }
    else
      dB = (usum - uleaves*_current->edge())/uleaves 
         - (dsum + dleaves*_current->edge())/dleaves;

    if ( dA * dB < 0) {
      xroot = 0.5 * (usum/uleaves - dsum/dleaves);    
      maxdist = std::max(xroot + dmax, umax - xroot);      
      if (first) {
        minxroot   = xroot;      	
      	minmaxdist = maxdist;
      	minnode    = _current;
        first = false;
      } else if (maxdist < minmaxdist) {
        minxroot   = xroot;      	
      	minmaxdist = maxdist;
      	minnode    = _current;
      }
    }   
          
    preorder();
  } while (_current != _pointer);

  float oldedge = minnode->edge();

  _current = minnode;
  root(nnodes+1);
  _pointer->left()->edge(minxroot);
  _pointer->right()->edge(oldedge - minxroot);
  
}



}


}
