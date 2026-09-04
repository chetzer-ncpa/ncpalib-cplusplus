#pragma once

#include "NCPA/atmosphere/Atmosphere1D.hpp"
#include "NCPA/atmosphere/Atmosphere2D.hpp"
#include "NCPA/atmosphere/Atmosphere3D.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/atmosphere/readers/abstract_atmosphere_reader_3d.hpp"
#include "NCPA/atmosphere/readers/ncpaprop_atmosphere_reader_3d.hpp"
#include "NCPA/files.hpp"

#include <fstream>
#include <iostream>
#include <memory>
#include <vector>

namespace NCPA {
    namespace atmos {
        class ncpaprop_atmosphere_reader_stratified_3d
            : public ncpaprop_atmosphere_reader_3d {
            public:
                ncpaprop_atmosphere_reader_stratified_3d() :
                    ncpaprop_atmosphere_reader_3d() {}

                ncpaprop_atmosphere_reader_stratified_3d(
                    const ncpaprop_atmosphere_reader_stratified_3d& other ) :
                    ncpaprop_atmosphere_reader_stratified_3d() {}

                ncpaprop_atmosphere_reader_stratified_3d(
                    ncpaprop_atmosphere_reader_stratified_3d&&
                        source ) noexcept :
                    ncpaprop_atmosphere_reader_stratified_3d() {
                    swap( *this, source );
                }

                virtual ~ncpaprop_atmosphere_reader_stratified_3d() {}

                ncpaprop_atmosphere_reader_stratified_3d& operator=(
                    ncpaprop_atmosphere_reader_stratified_3d other ) {
                    swap( *this, other );
                    return *this;
                }

                friend void swap(
                    ncpaprop_atmosphere_reader_stratified_3d& a,
                    ncpaprop_atmosphere_reader_stratified_3d& b ) noexcept {
                    using std::swap;
                    swap( static_cast<ncpaprop_atmosphere_reader_3d&>( a ),
                          static_cast<ncpaprop_atmosphere_reader_3d&>( b ) );
                }

                virtual std::unique_ptr<_abstract_atmosphere_reader_3d> clone()
                    const override {
                    return std::unique_ptr<_abstract_atmosphere_reader_3d>(
                        new ncpaprop_atmosphere_reader_stratified_3d(
                            *this ) );
                }

                virtual Atmosphere3D read(
                    const std::string& filename ) override {
                    ncpaprop_atmosphere_reader_1d reader;
                    Atmosphere1D atm1d = reader.read( filename );
                    Atmosphere3D atm3d = Atmosphere3D(
                        _atm_3d_ptr_t( new stratified_atmosphere_3d() ) );
                    atm3d.set( atm1d );
                    return atm3d;
                }

                virtual bool stratified() const override { return true; }
        };
    }  // namespace atmos
}  // namespace NCPA
