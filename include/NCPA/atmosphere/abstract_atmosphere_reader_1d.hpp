#pragma once

#ifndef FILE_SEPARATOR
#  ifdef _WIN32
#    define FILE_SEPARATOR '\\'
#  else
#    define FILE_SEPARATOR '/'
#  endif
#endif

#include "NCPA/atmosphere/Atmosphere1D.hpp"
// #include "NCPA/atmosphere/builders.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/files.hpp"

#include <fstream>
#include <iostream>
#include <memory>
#include <vector>

static void swap( NCPA::atmos::_abstract_atmosphere_reader_1d&,
                  NCPA::atmos::_abstract_atmosphere_reader_1d& ) noexcept;

namespace NCPA {
    namespace atmos {
        class _abstract_atmosphere_reader_1d {
            public:
                virtual ~_abstract_atmosphere_reader_1d() {}

                virtual std::unique_ptr<_abstract_atmosphere_reader_1d> clone()
                    const
                    = 0;

                virtual Atmosphere1D read( std::istream& in1 ) = 0;
                virtual Atmosphere1D read(
                    std::vector<std::istream *>& streams )
                    = 0;
                virtual Atmosphere1D read( const std::string& filename ) = 0;
                virtual Atmosphere1D read(
                    const std::vector<std::string>& filenames )
                    = 0;
        };
    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::_abstract_atmosphere_reader_1d& a,
                  NCPA::atmos::_abstract_atmosphere_reader_1d& b ) noexcept {}
