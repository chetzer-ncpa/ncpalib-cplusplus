#pragma once

#ifndef FILE_SEPARATOR
#  ifdef _WIN32
#    define FILE_SEPARATOR '\\'
#  else
#    define FILE_SEPARATOR '/'
#  endif
#endif

#include "NCPA/atmosphere/Atmosphere1D.hpp"
#include "NCPA/atmosphere/Atmosphere2D.hpp"
#include "NCPA/atmosphere/Atmosphere3D.hpp"
// #include "NCPA/atmosphere/builders.hpp"
#include "NCPA/atmosphere/abstract_atmosphere_reader_2d.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/atmosphere/ncpaprop_atmosphere_reader_1d.hpp"
#include "NCPA/files.hpp"

#include <fstream>
#include <iostream>
#include <memory>
#include <vector>

static void swap( NCPA::atmos::ncpaprop_atmosphere_reader_2d&,
                  NCPA::atmos::ncpaprop_atmosphere_reader_2d& ) noexcept;

namespace NCPA {
    namespace atmos {
        class ncpaprop_atmosphere_reader_2d
            : public _abstract_atmosphere_reader_2d {
            public:
                ncpaprop_atmosphere_reader_2d() :
                    _abstract_atmosphere_reader_2d(),
                    _axis_units { &NCPA::units::KILOMETERS } {}

                ncpaprop_atmosphere_reader_2d(
                    const ncpaprop_atmosphere_reader_2d& other ) :
                    ncpaprop_atmosphere_reader_2d() {
                    _axis_units = other._axis_units;
                }

                DECLARE_BOILERPLATE_METHODS( ncpaprop_atmosphere_reader_2d,
                                             _abstract_atmosphere_reader_2d )

                virtual bool stratified() const { return false; }

                virtual Atmosphere2D read(
                    const std::string& filename ) override {
                    // auto const pos = filename.find_last_of(
                    // NCPA::files::filesep() );
                    std::string basedir = NCPA::files::pathname( filename );
                    // if (pos == filename.npos) {
                    //     basedir = ".";
                    // } else {
                    //     basedir = filename.substr( 0, pos );
                    // }
                    std::ifstream summary( filename, std::ios_base::in );
                    vector_u_t ranges( 0, "km" );
                    std::vector<std::string> fns;
                    std::string line;
                    do {
                        std::getline( summary, line );
                        line = NCPA::strings::deblank( line );
                        if (line.length() > 0 && line[ 0 ] != '#') {
                            auto tokens = NCPA::strings::split( line );
                            if (tokens.size() != 2) {
                                std::ostringstream oss;
                                oss << "Atmosphere2D.read(): Malformed line:"
                                    << std::endl
                                    << "'" << line << "'" << std::endl
                                    << "Must be '<range> <filename>'";
                                throw std::invalid_argument( oss.str() );
                            }
                            ranges.push_back( std::stod( tokens[ 0 ] ) );
                            fns.push_back( NCPA::files::fullfile(
                                basedir, tokens[ 1 ] ) );
                        }
                    } while (!summary.eof());
                    summary.close();

                    std::vector<Atmosphere1D> atms;
                    ncpaprop_atmosphere_reader_1d reader;
                    for (auto it = fns.begin(); it != fns.end(); ++it) {
                        atms.push_back( reader.read( *it ) );
                    }
                    Atmosphere2D atm2d;
                    if (atms.size() == 1) {
                        atm2d = Atmosphere2D(
                            _atm_2d_ptr_t( new stratified_atmosphere_2d() ) );
                        atm2d.set( ranges, atms );
                    } else {
                        if (this->stratified()) {
                            atm2d = Atmosphere2D( _atm_2d_ptr_t(
                                new piecewise_stratified_atmosphere_2d() ) );
                        } else {
                            atm2d = Atmosphere2D(
                                _atm_2d_ptr_t( new grid_atmosphere_2d() ) );
                        }
                        atm2d.set( ranges, atms );
                    }
                    return atm2d;
                }

                virtual Atmosphere2D read(
                    const std::vector<std::string>& filenames ) override {
                    if (filenames.size() != 2) {
                        throw std::range_error(
                            "Filename vector must have 2 elements: summary "
                            "file, header file" );
                    }
                    std::ifstream summary( filenames[ 0 ], std::ios_base::in );
                    vector_u_t ranges( 0, _axis_units );
                    std::vector<std::string> fns;
                    std::string basedir
                        = NCPA::files::pathname( filenames[ 1 ] );
                    std::string line;
                    do {
                        std::getline( summary, line );
                        line = NCPA::strings::deblank( line );
                        if (line.length() > 0 && line[ 0 ] != '#') {
                            auto tokens = NCPA::strings::split( line );
                            if (tokens.size() != 2) {
                                std::ostringstream oss;
                                oss << "ncpaprop_atmosphere_reader_2d.read(): "
                                       "Malformed line:"
                                    << std::endl
                                    << "'" << line << "'" << std::endl
                                    << "Must be '<range> <filename>'";
                                throw std::invalid_argument( oss.str() );
                            }
                            ranges.push_back( std::stod( tokens[ 0 ] ) );
                            fns.push_back( NCPA::files::fullfile(
                                basedir, tokens[ 1 ] ) );
                        }
                    } while (!summary.eof());
                    summary.close();

                    std::vector<Atmosphere1D> atms;
                    ncpaprop_atmosphere_reader_1d reader;
                    for (auto it = fns.cbegin(); it != fns.cend(); ++it) {
                        // std::vector<const std::string&> filepair { *it,
                        //                                      filenames[ 1 ]
                        //                                      };
                        atms.push_back( reader.read( std::vector<std::string> {
                            *it, filenames[ 1 ] } ) );
                    }
                    Atmosphere2D atm2d;
                    if (atms.size() == 1) {
                        atm2d = Atmosphere2D(
                            _atm_2d_ptr_t( new stratified_atmosphere_2d() ) );
                        atm2d.set( ranges, atms );
                    } else {
                        if (this->stratified()) {
                            atm2d = Atmosphere2D( _atm_2d_ptr_t(
                                new piecewise_stratified_atmosphere_2d() ) );
                        } else {
                            atm2d = Atmosphere2D(
                                _atm_2d_ptr_t( new grid_atmosphere_2d() ) );
                        }
                        atm2d.set( ranges, atms );
                    }
                    return atm2d;
                }

                virtual _abstract_atmosphere_reader_2d& set_axis_units(
                    units_ptr_t u ) override {
                    _axis_units = u;
                    return static_cast<_abstract_atmosphere_reader_2d&>(
                        *this );
                }

            protected:
                units_ptr_t _axis_units;
        };
    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::ncpaprop_atmosphere_reader_2d& a,
                  NCPA::atmos::ncpaprop_atmosphere_reader_2d& b ) noexcept {
    ::swap( dynamic_cast<NCPA::atmos::_abstract_atmosphere_reader_2d&>( a ),
            dynamic_cast<NCPA::atmos::_abstract_atmosphere_reader_2d&>( b ) );
}
