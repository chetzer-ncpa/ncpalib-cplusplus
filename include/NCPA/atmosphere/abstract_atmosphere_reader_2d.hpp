#pragma once

#ifndef FILE_SEPARATOR
#ifdef _WIN32
#define FILE_SEPARATOR '\\'
#else
#define FILE_SEPARATOR '/'
#endif
#endif

#include "NCPA/atmosphere/Atmosphere1D.hpp"
#include "NCPA/atmosphere/Atmosphere2D.hpp"
#include "NCPA/atmosphere/Atmosphere3D.hpp"
// #include "NCPA/atmosphere/builders.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/files.hpp"

#include <fstream>
#include <iostream>
#include <memory>
#include <vector>

static void swap( NCPA::atmos::_abstract_atmosphere_reader_2d&,
                  NCPA::atmos::_abstract_atmosphere_reader_2d& ) noexcept;

namespace NCPA {
    namespace atmos {
        class _abstract_atmosphere_reader_2d {
            public:
                virtual ~_abstract_atmosphere_reader_2d() {}

                virtual std::unique_ptr<_abstract_atmosphere_reader_2d> clone()
                    const
                    = 0;

                virtual Atmosphere2D read( const std::string& filename ) = 0;
                virtual Atmosphere2D read(
                    const std::vector<std::string>& filenames )
                    = 0;
                virtual _abstract_atmosphere_reader_2d& set_axis_units(
                    units_ptr_t u )
                    = 0;
        };
    }
}

static void swap( NCPA::atmos::_abstract_atmosphere_reader_2d& a,
                  NCPA::atmos::_abstract_atmosphere_reader_2d& b ) noexcept {}