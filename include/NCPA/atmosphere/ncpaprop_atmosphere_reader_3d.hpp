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
#include "NCPA/atmosphere/abstract_atmosphere_reader_3d.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/atmosphere/ncpaprop_atmosphere_reader_1d.hpp"
#include "NCPA/files.hpp"

#include <fstream>
#include <iostream>
#include <memory>
#include <vector>

static void swap( NCPA::atmos::ncpaprop_atmosphere_reader_3d&,
                  NCPA::atmos::ncpaprop_atmosphere_reader_3d& ) noexcept;

namespace NCPA {
    namespace atmos {
        class ncpaprop_atmosphere_reader_3d
            : public _abstract_atmosphere_reader_3d {
            public:
                ncpaprop_atmosphere_reader_3d() :
                    _abstract_atmosphere_reader_3d(),
                    _axis_units {
                        { &NCPA::units::KILOMETERS, &NCPA::units::KILOMETERS }
                } {}

                ncpaprop_atmosphere_reader_3d(
                    const ncpaprop_atmosphere_reader_3d& other ) :
                    ncpaprop_atmosphere_reader_3d() {
                    _axis_units = other._axis_units;
                }

                ncpaprop_atmosphere_reader_3d(
                    ncpaprop_atmosphere_reader_3d&& source ) noexcept :
                    ncpaprop_atmosphere_reader_3d() {
                    ::swap( *this, source );
                }

                virtual ~ncpaprop_atmosphere_reader_3d() {}

                ncpaprop_atmosphere_reader_3d& operator=(
                    ncpaprop_atmosphere_reader_3d other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap(
                    ncpaprop_atmosphere_reader_3d& a,
                    ncpaprop_atmosphere_reader_3d& b ) noexcept;

                virtual std::unique_ptr<_abstract_atmosphere_reader_3d> clone()
                    const override {
                    return std::unique_ptr<_abstract_atmosphere_reader_3d>(
                        new ncpaprop_atmosphere_reader_3d( *this ) );
                }

                virtual bool stratified() const { return false; }

                virtual Atmosphere3D read(
                    const std::string& filename ) override {
                    vector_u_t x1( 0, _axis_units[ 0 ] ),
                        x2( 0, _axis_units[ 1 ] );
                    std::vector<size_t> x1inds, x2inds;
                    NCPA::arrays::ndvector<2, std::string> fns;
                    this->_parse_summary_file( filename, x1inds, x1, x2inds,
                                               x2, fns );
                    // std::cout << "Done parsing summary file" << std::endl;
                    NCPA::arrays::ndvector<2, Atmosphere1D> atms;
                    atms.reshape( fns.shape() );
                    // std::cout << "Created atmosphere 2d matrix of shape {"
                    //           << atms.shape()[ 0 ] << ", " << atms.shape()[
                    //           1 ]
                    //           << "}" << std::endl;
                    ncpaprop_atmosphere_reader_1d reader;
                    vector_u_t ux1 = x1, ux2 = x2;
                    ux1.resize( atms.shape()[ 0 ] );
                    ux2.resize( atms.shape()[ 1 ] );
                    for (size_t n = 0; n < x1.size(); ++n) {
                        Atmosphere1D atm1d
                            = reader.read( fns[ x1inds[ n ] ][ x2inds[ n ] ] );
                        atms[ x1inds[ n ] ][ x2inds[ n ] ] = atm1d;
                        ux1[ x1inds[ n ] ]                 = x1[ n ];
                        ux2[ x2inds[ n ] ]                 = x2[ n ];
                    }
                    Atmosphere3D atm3d;
                    if (x1.size() == 1 && x2.size() == 1) {
                        atm3d = Atmosphere3D(
                            _atm_3d_ptr_t( new stratified_atmosphere_3d() ) );
                        atm3d.set( atms[ 0 ][ 0 ] );
                    } else {
                        atm3d = Atmosphere3D(
                            _atm_3d_ptr_t( new grid_atmosphere_3d() ) );
                        atm3d.set( ux1, ux2, atms );
                    }
                    return atm3d;
                }

                virtual Atmosphere3D read(
                    const std::vector<std::string>& filenames ) override {
                    if (filenames.size() != 2) {
                        throw std::range_error(
                            "Filename vector must have 2 elements: summary "
                            "file, header file" );
                    }
                    vector_u_t x1( 0, _axis_units[ 0 ] ),
                        x2( 0, _axis_units[ 1 ] );
                    std::vector<size_t> x1inds, x2inds;
                    NCPA::arrays::ndvector<2, std::string> fns;
                    this->_parse_summary_file( filenames[ 0 ], x1inds, x1,
                                               x2inds, x2, fns );
                    NCPA::arrays::ndvector<2, Atmosphere1D> atms;
                    atms.reshape( fns.shape() );
                    ncpaprop_atmosphere_reader_1d reader;
                    vector_u_t ux1 = x1, ux2 = x2;
                    ux1.resize( atms.shape()[ 0 ] );
                    ux2.resize( atms.shape()[ 1 ] );
                    for (size_t n = 0; n < fns.size(); ++n) {
                        atms[ x1inds[ n ] ][ x2inds[ n ] ]
                            = reader.read( { fns[ x1inds[ n ] ][ x2inds[ n ] ],
                                             filenames[ 1 ] } );
                        ux1[ x1inds[ n ] ] = x1[ n ];
                        ux2[ x2inds[ n ] ] = x2[ n ];
                    }
                    Atmosphere3D atm3d;
                    if (x1.size() && x2.size() == 1) {
                        atm3d = Atmosphere3D(
                            _atm_3d_ptr_t( new stratified_atmosphere_3d() ) );
                        atm3d.set( atms[ 0 ][ 0 ] );
                    } else {
                        atm3d = Atmosphere3D(
                            _atm_3d_ptr_t( new grid_atmosphere_3d() ) );
                        atm3d.set( ux1, ux2, atms );
                    }
                    return atm3d;
                }

                virtual _abstract_atmosphere_reader_3d& set_axis_units(
                    size_t n, units_ptr_t u ) override {
                    _axis_units[ n ] = u;
                    return static_cast<_abstract_atmosphere_reader_3d&>(
                        *this );
                }

            protected:
                std::array<units_ptr_t, 2> _axis_units;

                void _parse_summary_file(
                    const std::string& filename, std::vector<size_t>& x1inds,
                    vector_u_t& x1, std::vector<size_t>& x2inds,
                    vector_u_t& x2,
                    NCPA::arrays::ndvector<2, std::string>& fns ) {
                    std::string basedir = NCPA::files::pathname( filename );
                    std::string line;
                    std::ifstream summary( filename, std::ios_base::in );
                    std::vector<std::string> tmpfilenames;
                    do {
                        std::getline( summary, line );
                        line = NCPA::strings::deblank( line );
                        if (line.length() > 0 && line[ 0 ] != '#') {
                            auto tokens = NCPA::strings::split( line );
                            // std::cout << "Parsed line into { ";
                            // NCPA::arrays::write( std::cout, tokens );
                            // std::cout << " }" << std::endl;
                            if (tokens.size() != 5) {
                                std::ostringstream oss;
                                oss << "ncpaprop_atmosphere_reader_3d.read(): "
                                       "Malformed line:"
                                    << std::endl
                                    << line << std::endl
                                    << "Must be '<x1 index> <x1> <x2 index> "
                                       "<x2> <filename>'";
                                throw std::invalid_argument( oss.str() );
                            }
                            x1inds.push_back( std::stoi( tokens[ 0 ] ) );
                            x1.push_back( std::stod( tokens[ 1 ] ) );
                            x2inds.push_back( std::stoi( tokens[ 2 ] ) );
                            x2.push_back( std::stod( tokens[ 3 ] ) );
                            tmpfilenames.push_back( NCPA::files::fullfile(
                                basedir, tokens[ 4 ] ) );
                        }
                    } while (!summary.eof());
                    summary.close();

                    std::vector<size_t> ux1inds = x1inds, ux2inds = x2inds;
                    std::sort( ux1inds.begin(), ux1inds.end() );
                    std::sort( ux2inds.begin(), ux2inds.end() );
                    size_t nx1 = std::distance(
                        ux1inds.begin(),
                        std::unique( ux1inds.begin(), ux1inds.end() ) );
                    size_t nx2 = std::distance(
                        ux2inds.begin(),
                        std::unique( ux2inds.begin(), ux2inds.end() ) );
                    fns.reshape( { nx1, nx2 } );
                    for (size_t i = 0; i < tmpfilenames.size(); ++i) {
                        fns[ x1inds[ i ] ][ x2inds[ i ] ] = tmpfilenames[ i ];
                    }
                }
        };
    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::ncpaprop_atmosphere_reader_3d& a,
                  NCPA::atmos::ncpaprop_atmosphere_reader_3d& b ) noexcept {
    ::swap( dynamic_cast<NCPA::atmos::_abstract_atmosphere_reader_3d&>( a ),
            dynamic_cast<NCPA::atmos::_abstract_atmosphere_reader_3d&>( b ) );
}
