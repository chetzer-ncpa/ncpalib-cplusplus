#pragma once

#include "NCPA/atmosphere/Atmosphere1D.hpp"
#include "NCPA/atmosphere/Atmosphere2D.hpp"
#include "NCPA/atmosphere/Atmosphere3D.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/atmosphere/readers/ncpaprop_atmosphere_reader_2d.hpp"
#include "NCPA/files.hpp"

#include <fstream>
#include <iostream>
#include <memory>
#include <vector>

namespace NCPA {
    namespace atmos {
        class ncpaprop_atmosphere_reader_stratified_2d
            : public ncpaprop_atmosphere_reader_2d {
            public:
                ncpaprop_atmosphere_reader_stratified_2d() :
                    ncpaprop_atmosphere_reader_2d() {}

                ncpaprop_atmosphere_reader_stratified_2d(
                    const ncpaprop_atmosphere_reader_stratified_2d& other ) :
                    ncpaprop_atmosphere_reader_2d() {}

                ncpaprop_atmosphere_reader_stratified_2d(
                    ncpaprop_atmosphere_reader_stratified_2d&&
                        source ) noexcept :
                    ncpaprop_atmosphere_reader_stratified_2d() {
                    swap( *this, source );
                }

                virtual ~ncpaprop_atmosphere_reader_stratified_2d() {}

                ncpaprop_atmosphere_reader_stratified_2d& operator=(
                    ncpaprop_atmosphere_reader_stratified_2d other ) {
                    swap( *this, other );
                    return *this;
                }

                friend void swap(
                    ncpaprop_atmosphere_reader_stratified_2d& a,
                    ncpaprop_atmosphere_reader_stratified_2d& b ) noexcept {
                    using std::swap;
                    swap( static_cast<ncpaprop_atmosphere_reader_2d&>( a ),
                          static_cast<ncpaprop_atmosphere_reader_2d&>( b ) );
                }

                virtual std::unique_ptr<_abstract_atmosphere_reader_2d> clone()
                    const override {
                    return std::unique_ptr<_abstract_atmosphere_reader_2d>(
                        new ncpaprop_atmosphere_reader_stratified_2d(
                            *this ) );
                }

                virtual bool stratified() const override { return true; }
        };
    }  // namespace atmos
}  // namespace NCPA
