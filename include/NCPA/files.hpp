#pragma once

#include <string>
#include <stdexcept>
#include <sstream>
#include <cstdio>
#include <fstream>

namespace NCPA::files {

	/**
	@brief Counts the rows in a file.
	@input filename The filename to count the rows from.
	@returns The number of rows in the file.
	*/
	int count_rows_arbcol(const std::string& filename) {
	  int answer,c;
	  FILE *f=fopen(filename.c_str(),"r");

	  if(f==NULL) {
	      std::ostringstream es;
	      es << "file << " << filename << " could not be opened.\n";
	      throw std::invalid_argument(es.str());
	  }
	  answer=0;
	  //read_header(f);
	  while((c=getc(f))!=EOF){
	      if(c=='\n') answer=answer+1;
	  }
	  fclose(f);
	  return answer;
	}

	/**
	Determines if a given filename represents a readable file by
	attempting to open it and determining if it is in a good state.
	@brief Determines if a given filename represents a readable file.
	@param filename The name of the file to attempt to open.
	@returns true if the file can be read on open, false otherwise.
	*/
	bool fexists( const char *filename ) {
		std::ifstream ifile( filename );
		bool tf = ifile.good();
		ifile.close();
		return tf;
	}

	/**
	@brief Prints two arbitrary columns to a stream with specified delimiter.
	@param out The stream to output to.
	@param nz The number of values in each array.
	@param col1 The first column.
	@param col2 The second column.
	@param del The delimiter.  Defaults to a single space.
	*/
	template<typename T, typename U>
	void print_2_columns( std::ostream &out, size_t nz,
		T *col1, U *col2, std::string del = " " ) {
		for (size_t i = 0; i < nz; i++) {
			out << col1[i] << del << col2[i] << std::endl;
		}
	}

}
