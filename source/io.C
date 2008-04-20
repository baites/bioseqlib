#include <io.h>


namespace biolib {


bool load(std::string s, biolib::alignment& align)
{

  seqio::SEQFILE *sfp;

  if ((sfp = seqio::seqfopen(const_cast<char*>(s.c_str()), "r", NULL)) == NULL) 
  {
	std::cout << "ERROR: Format file cannot be recognized !" << std::endl;
	return false;
  }

  int len;
  char *tmpseq;
  
  align.clear();
   
  while ((tmpseq = seqio::seqfgetseq(sfp, &len, 0)) != NULL)
     if (tmpseq[0] != '\0') 
         align.add(seqio::seqfidlist(sfp,0), tmpseq);

  seqio::seqfclose(sfp);

  return true;
}


bool loadraw(std::string s, biolib::alignment& align)
{

  seqio::SEQFILE *sfp;

  if ((sfp = seqio::seqfopen(const_cast<char*>(s.c_str()), "r", NULL)) == NULL) {

	std::cout << "ERROR: Format file cannot be recognized !" << std::endl;
 	return false;
  }

  int len;
  char *tmpseq;
   
  align.clear();

  while ((tmpseq = seqio::seqfgetrawseq(sfp, &len, 0)) != NULL)
     if (tmpseq[0] != '\0') 
         align.add(seqio::seqfidlist(sfp,0), tmpseq);

  seqio::seqfclose(sfp);

  return true;
}


}
