#pragma once

#include "NCPA/atmosphere/abstract_atmosphere_reader_2d.hpp"
#include "NCPA/atmosphere/Atmosphere1D.hpp"
#include "NCPA/atmosphere/Atmosphere2D.hpp"
#include "NCPA/atmosphere/Atmosphere3D.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/files.hpp"

#include <fstream>
#include <iostream>
#include <memory>
#include <vector>

namespace NCPA {
    namespace atmos {
        class AtmosphereReader2D {
            public:
                AtmosphereReader2D() {}

                AtmosphereReader2D(
                    std::unique_ptr<_abstract_atmosphere_reader_2d> ptr ) :
                    AtmosphereReader2D() {
                    _ptr = std::move( ptr );
                }

                // copy constructor
                AtmosphereReader2D( const AtmosphereReader2D& other ) :
                    AtmosphereReader2D() {
                    _ptr = std::move( other._ptr->clone() );
                }

                AtmosphereReader2D( AtmosphereReader2D&& source ) noexcept :
                    AtmosphereReader2D() {
                    swap( *this, source );
                }

                virtual ~AtmosphereReader2D() {}

                AtmosphereReader2D& operator=( AtmosphereReader2D other ) {
                    swap( *this, other );
                    return *this;
                }

                friend void swap( AtmosphereReader2D& a,
                                  AtmosphereReader2D& b ) noexcept {
                    using std::swap;
                    swap( a._ptr, b._ptr );
                }

                virtual Atmosphere2D read( const std::string& filename,
                                           bool stratified = false ) {
                    _check_pointer();
                    return _ptr->read( filename );
                }

                virtual Atmosphere2D read(
                    const std::vector<std::string>& filenames,
                    bool stratified = false ) {
                    _check_pointer();
                    return _ptr->read( filenames );
                }

            protected:
                void _check_pointer() const {
                    if (!_ptr) {
                        throw std::logic_error(
                            "AtmosphereReader2D: No pointer set!" );
                    }
                }

                std::unique_ptr<_abstract_atmosphere_reader_2d> _ptr;
        };
    }  // namespace atmos
}  // namespace NCPA
