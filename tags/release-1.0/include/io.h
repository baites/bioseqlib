#ifndef BIOLIB_IO_INCLUDED
#define BIOLIB_IO_INCLUDED

#include <iostream>
#include <alignment.h>
#include <seqio.h>

namespace biolib {

//! loads a set of sequences from a file.
bool load(std::string s, biolib::alignment& align);

//! load a set of sequences from a file (including gaps).
bool loadraw(std::string s, biolib::alignment& align);

}

#endif

