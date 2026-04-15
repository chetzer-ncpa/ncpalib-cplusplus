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
#include "NCPA/atmosphere/abstract_atmosphere_reader_1d.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/files.hpp"

#include <fstream>
#include <iostream>
#include <memory>
#include <vector>

static void swap( NCPA::atmos::AtmosphereReader1D&,
                  NCPA::atmos::AtmosphereReader1D& ) noexcept;

namespace NCPA {
    namespace atmos {
        class AtmosphereReader1D {
            public:
                AtmosphereReader1D() {}

                AtmosphereReader1D(
                    std::unique_ptr<_abstract_atmosphere_reader_1d> ptr ) :
                    AtmosphereReader1D() {
                    _ptr = std::move( ptr );
                }

                // copy constructor
                AtmosphereReader1D( const AtmosphereReader1D& other ) :
                    AtmosphereReader1D() {
                    _ptr = std::move( other._ptr->clone() );
                }

                DECLARE_WRAPPER_BOILERPLATE_METHODS( AtmosphereReader1D )

                virtual Atmosphere1D read( std::istream& in1 ) {
                    _check_pointer();
                    return _ptr->read( in1 );
                }

                virtual Atmosphere1D read(
                    std::vector<std::istream *>& streams ) {
                    _check_pointer();
                    return _ptr->read( streams );
                }

                virtual Atmosphere1D read( const std::string& filename ) {
                    _check_pointer();
                    return _ptr->read( filename );
                }

                virtual Atmosphere1D read(
                    const std::vector<std::string>& filenames ) {
                    _check_pointer();
                    return _ptr->read( filenames );
                }

            protected:
                void _check_pointer() const {
                    if (!_ptr) {
                        throw std::logic_error(
                            "AtmosphereReader1D: No pointer set!" );
                    }
                }

                std::unique_ptr<_abstract_atmosphere_reader_1d> _ptr;
        };
    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::AtmosphereReader1D& a,
                  NCPA::atmos::AtmosphereReader1D& b ) noexcept {
    using std::swap;
    std::swap( a._ptr, b._ptr );
}
