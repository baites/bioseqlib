#include <trees/neighbour_joining_tree.h>


template<typename I, typename J>
std::size_t random_tree_spliter (biolib::trees::binary_tree<I,J> * tree, std::vector<I> & leaves)
{
  typename biolib::trees::binary_tree<I,J>::iterator itr (tree->begin());  
  typename biolib::trees::binary_tree<I,J>::iterator stop(tree->end());  

  std::srand(time(NULL));

  int rcount = (std::rand() % (tree->size()-3))+1;      
  leaves.clear();

  for (int i=0; i<rcount; i++) itr++;    
  while(itr->leaf()) itr++;
    
  biolib::trees::binary_tree<I,J> half_tree(itr.pointer());
  
  if (itr->parent()->left() == itr.pointer())
    itr->parent()->remove_left();
  else if (itr->parent()->right() == itr.pointer())
    itr->parent()->remove_right();
     
  half_tree.print();   
  tree->print();
     
  for(itr=half_tree.begin(); itr!=half_tree.end(); ++itr) 
    if(itr->leaf()) leaves.push_back(itr->data());
  if(itr->leaf()) leaves.push_back(itr->data());
    
  std::size_t pointer = leaves.size();
  
  for(itr=tree->begin(); itr!=tree->end(); ++itr) 
    if(itr->leaf()) leaves.push_back(itr->data());
  if(itr->leaf()) leaves.push_back(itr->data());

  return pointer;
}


main()
{
  std::vector<std::vector<float> > dists(6, std::vector<float>(6,0));
  
  dists[0][1] = dists[1][0] =  5;
  dists[0][2] = dists[2][0] =  6;
  dists[0][3] = dists[3][0] = 11;
  dists[0][4] = dists[4][0] =  9;
  dists[0][5] = dists[5][0] =  9;

  dists[1][2] = dists[2][1] =  5;
  dists[1][3] = dists[3][1] = 10;
  dists[1][4] = dists[4][1] =  8;
  dists[1][5] = dists[5][1] =  8;

  dists[2][3] = dists[3][2] =  7;
  dists[2][4] = dists[4][2] =  5;
  dists[2][5] = dists[5][2] =  5;

  dists[3][4] = dists[4][3] =  8;
  dists[3][5] = dists[5][3] =  8;

  dists[4][5] = dists[5][4] =  2;

  std::cout << "Distances: " << std::endl;
  for(int i=0; i<dists.size(); i++) {
    for(int j=0; j<dists.size(); j++)
      std::cout << std::setw(4) << dists[i][j];
    std::cout << std::endl; 
  }  

  std::cout << std::endl;
  
  biolib::trees::neighbour_joining_tree tree(dists);

  tree.rooted();

  tree.print();

  std::vector<int> leaves;

  std::cout << random_tree_spliter(&tree, leaves) << std::endl;

  for(std::size_t i=0; i<leaves.size(); i++)
    std::cout << leaves[i] << " ";
  
  std::cout << std::endl;
}
