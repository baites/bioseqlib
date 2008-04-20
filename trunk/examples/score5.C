/*********************************************************************************************************

 Library name     :   Biolib++ (Rutines for biological applications)
 File name        :   trees7.C
 Defined Objects  :   --  

 Programmer       :   Victor E. Bazterra 
 Directed by      :   Julio C. Facelli
                      
 Initialized      :   Wed Nov 12 10:59:47 MST 2003  
 Last updated     :   Wed Nov 12 10:59:47 MST 2003
   
 DESCRIPTION:
    A examples of weighting a tree using branch_proportional_weighting object.
 
**********************************************************************************************************/

#include <iostream>
#include <vector>
#include <iomanip>
#include <alignment.h>
#include <amino_acids.h>
#include <scores/alignment_score.h>
#include <scores/two_sequences_clusterw.h>
#include <trees/neighbour_joining_tree.h>
#include <trees/branch_proportional_weighting.h>

int main() 
{
  // Define an aminoacid alphabet.
  biolib::amino_acid_alphabet alphabet;

  // Define a alignment object that can hold the alphabet.
  biolib::alignment alignment(alphabet);
 
  // Define another aligment object to hold the initial sequences.
  biolib::alignment sequences(alphabet);

  // Add the initial sequences.
  sequences.add("VSP1_AGKCO","VIGGDECNINEHRFLALVYANGSLCGGTLINQEWVLTARHCDRGNMRIYLGMHNLKVLNKDALRRFPKEKYFCLNTRNDTIWDKDIMLIRLNRPVRNSAHIAPLSLPSNPPSVGSVCRIMGWGTITSPNATLPDVPHCANINILDYAVCQAAYKGLAATTLCAGILEGGKDTCKGDSGGPLICNGQFQGILSVGGNPCAQPRKPGIYTKVFDYTDWDYTDWIQSIISGNTDATCPP");

  sequences.add("ACRO_RAT  ","MVEMLPTVVALVLAVSVVAKDNTTCDGPCGLRFRQNPQAGIRIVGGQTSSRWAWPWMVSLQIFTSHNSRRYHACGGSLLNSHWVLTAAHCFDNKKKVYDWRLVFGAHEIEYGRNKPVKEPQQERYVQKIVIHEKYNAVTEGNDIALLKVTPPVTCGDFVGPGCLPHFKSGPPRIPHTCYVTGWGYIKDNAPRPSPVLMEARVDLIDLDLCNSTQWYNGRVTSTNVCAGYPEGKIDTCQGDSGGPLMCRDTRRQPFVIVGITSWGVGCARAKRPGVYTATWDYLDWIASKIGPTALHLIQPATPHPPTTQQPVISFHPPSTPPSLVLPTPVSSAALPTPPRPLLHQPSSVHTSSAPVIPLLSLLTPVQPVSFTLAAYHTRHHTTLSFASALQHLIEALKMRTYPIKYPSRYSGPVNYQHRFSTFEPLSNKPSEPLLHS");

  sequences.add("VSP1_AGKRH","VIGGDECNINEHRFLVAVYEGTNWTFICGGVLIHPEWVITAEHCARRRMNLVFGMHRKSEKFDDEQERYPKKRYFIRCNKTRTSWDEDIMLIRLNKPVNNSEHIAPLSLPSNPPIVGSDCRVMGWGSINRRIDVLSDEPRCANINLHNFTMCHGLFRKMPKKGRVLCAGDLRGRRDSCNSDSGGPLICNEELHGIVARGPNPCAQPNKPALYTSIYDYRDWVNNVIAGNATCSP");  

  // Add the alignment to be evaluated.
  alignment.add("VSP1_AGKCO","-----VIGGDECNINEHRFLALVYANGSLCGGTLINQEWVLTARHCDRGNMRIYLGMHNLKVLNKDALRRFPKEKYFCLNTRNDTIWDKDIMLIRLNRPVRNSAHIAPLSLPSNPPSVGSVCRIMGWGTITSPNATLPDVPHCANINILDYAVCQAAYKGLAATTLCAGILEGGKDTCKGDSGGPLICNGQFQGILSVGGNPCAQPRKPGIYTKVFDYTDWIQSIISGNTDATCPP---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");

  alignment.add("ACRO_RAT  ","------------------------------------------------MVEMLPTVVALVLAVSVVAKDNTTCDGPCGLRFRQNPQAGIRIVGGQTSSRWAWPWMVSLQIFTSHNSRRYHACGGSLLNSHWVLTAAHCFDNKKKVYDWRLVFGAHEIEYGRNKPVKEPQQERYVQKIVIHEKYNAVTEGNDIALLKVTPPVTCGDFVGPGCLPHFKSGPPRIPHTCYVTGWGYIKDNAPRPSPVLMEARVDLIDLDLCNSTQWYNGRVTSTNVCAGYPEGKIDTCQGDSGGPLMCRDTRRQPFVIVGITSWGVGCARAKRPGVYTATWDYLDWIASKIGPTALHLIQPATPHPPTTQQPVISFHPPSTPPSLVLPTPVSSAALPTPPRPLLHQPSSVHTSSAPVIPLLSLLTPVQPVSFTLAAYHTRHHTTLSFASALQHLIEALKMRTYPIKYPSRYSGPVNYQHRFSTFEPLSNKPSEPLLHS"); 

  alignment.add("VSP1_AGKRH","VIGGDECNINEHRFLVAVYEGTNWTFICGGVLIHPEWVITAEHCARRRMNLVFGMHRKSEKFDDEQERYPKKRYFIRCNKTRTSWDEDIMLIRLNKPVNNSEHIAPLSLPSNPPIVGSDCRVMGWGSINRRIDVLSDEPRCANINLHNFTMCHGLFRKMPKKGRVLCAGDLRGRRDSCNSDSGGPLICNEELHGIVARGPNPCAQPNKPALYTSIYDYRDWVNNVIAGNATCSP-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");

  // Print both alignmen objects.
  std::cout << "Sequences related to aligment:" << std::endl;
  std::cout << sequences;
  std::cout << "Aligment to be evaluated:" << std::endl;
  std::cout << alignment;

  // Define the similarity using BLOSUM 50 matrix.
  biolib::similarity similarity(alphabet, biolib::blosum50mt, "BLOSUM 50");

  // Calculate the distance matrix betweem sequences.
  biolib::scores::two_sequences_clusterw scoring(similarity, 12, 2);  
  biolib::matrix distances(scoring.evaluate_matrix(sequences));
  
  // Printing the distances.
  std::cout << "Distances: " << std::endl;
  for(int i=0; i<distances.size(); i++) {
    for(int j=0; j<distances[i].size(); j++)
      std::cout << std::setw(10) << distances[i][j];
    std::cout << std::endl; 
  }  

  // Binary tree created by Neighbour Joining Tree
  biolib::trees::neighbour_joining_tree tree(distances);

  // Creates a root node for the tree.
  tree.rooted();

  // Print the trees.
  tree.print();

  // Calculate the sequence weights from the trees.
  biolib::trees::branch_proportional_weighting<int> weights(&tree);
  weights.solve();
  
  // Print the weights.
  std::cout << "Sequence weights:" << std::endl;
  std::vector<float> const & w = weights.get_weights();
  for(int i=0; i<w.size(); i++) 
    std::cout << w[i] << "  ";
  std::cout << std::endl << std::endl;

  // Calculate aligment score.
  biolib::scores::alignment_score score(similarity, 12, 2);

  // If this line is comented, the score is not weighted.
  // score.set_weights(w);

  // Calculare the score.
  score.evaluate(alignment);

  return 0;  
}

