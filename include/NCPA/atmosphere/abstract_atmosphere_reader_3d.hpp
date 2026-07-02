#pragma once

/**
 * @file
 * @brief 
 */

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
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/files.hpp"

#include <fstream>
#include <iostream>
#include <memory>
#include <vector>

/**
 * @brief 
 * @param  a 
 * @param  b 
 */
static void swap( NCPA::atmos::_abstract_atmosphere_reader_3d&,
                  NCPA::atmos::_abstract_atmosphere_reader_3d& ) noexcept;

namespace NCPA {
    namespace atmos {

        /**
         * @brief 
         */
        class _abstract_atmosphere_reader_3d {
            public:
                /**
                 * @brief 
                 */
                virtual ~_abstract_atmosphere_reader_3d() {}

                /**
                 * @brief 
                 * @return std::unique_ptr<_abstract_atmosphere_reader_3d> 
                 */
                virtual std::unique_ptr<_abstract_atmosphere_reader_3d> clone()
                    const
                    = 0;

                /**
                 * @brief 
                 * @param filename 
                 * @return Atmosphere3D 
                 */
                virtual Atmosphere3D read( const std::string& filename ) = 0;

                /**
                 * @brief 
                 * @param filenames 
                 * @return Atmosphere3D 
                 */
                virtual Atmosphere3D read(
                    const std::vector<std::string>& filenames )
                    = 0;

                /**
                 * @brief 
                 * @param n 
                 * @param u 
                 * @return _abstract_atmosphere_reader_3d& 
                 */
                virtual _abstract_atmosphere_reader_3d& set_axis_units(
                    size_t n, units_ptr_t u )
                    = 0;
        };
    }  // namespace atmos
}  // namespace NCPA

/**
 * @brief 
 * @param a 
 * @param b 
 */
static void swap( NCPA::atmos::_abstract_atmosphere_reader_3d& a,
                  NCPA::atmos::_abstract_atmosphere_reader_3d& b ) noexcept {}