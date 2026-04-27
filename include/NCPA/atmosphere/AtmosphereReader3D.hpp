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
#include "NCPA/files.hpp"

#include <fstream>
#include <iostream>
#include <memory>
#include <vector>

static void swap( NCPA::atmos::AtmosphereReader3D&,
                  NCPA::atmos::AtmosphereReader3D& ) noexcept;

namespace NCPA {
    namespace atmos {
        class AtmosphereReader3D {
            public:
                AtmosphereReader3D() {}

                AtmosphereReader3D(
                    std::unique_ptr<_abstract_atmosphere_reader_3d> ptr ) :
                    AtmosphereReader3D() {
                    _ptr = std::move( ptr );
                }

                // copy constructor
                AtmosphereReader3D( const AtmosphereReader3D& other ) :
                    AtmosphereReader3D() {
                    _ptr = std::move( other._ptr->clone() );
                }

                AtmosphereReader3D( AtmosphereReader3D&& source ) noexcept :
                    AtmosphereReader3D() {
                    ::swap( *this, source );
                }

                virtual ~AtmosphereReader3D() {}

                AtmosphereReader3D& operator=( AtmosphereReader3D other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap( AtmosphereReader3D& a,
                                    AtmosphereReader3D& b ) noexcept;

                virtual Atmosphere3D read( const std::string& filename,
                                           bool stratified = false ) {
                    _check_pointer();
                    return _ptr->read( filename );
                }

                virtual Atmosphere3D read(
                    const std::vector<std::string>& filenames,
                    bool stratified = false ) {
                    _check_pointer();
                    return _ptr->read( filenames );
                }

                virtual AtmosphereReader3D& set_axis_units( size_t n,
                                                            units_ptr_t u ) {
                    _check_pointer();
                    _ptr->set_axis_units( n, u );
                    return *this;
                }

            protected:
                void _check_pointer() const {
                    if (!_ptr) {
                        throw std::logic_error(
                            "AtmosphereReader3D: No pointer set!" );
                    }
                }

                std::unique_ptr<_abstract_atmosphere_reader_3d> _ptr;
        };
    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::AtmosphereReader3D& a,
                  NCPA::atmos::AtmosphereReader3D& b ) noexcept {
    using std::swap;
    std::swap( a._ptr, b._ptr );
}
