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
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/atmosphere/ncpaprop_atmosphere_reader_2d.hpp"
#include "NCPA/files.hpp"

#include <fstream>
#include <iostream>
#include <memory>
#include <vector>

static void swap(
    NCPA::atmos::ncpaprop_atmosphere_reader_stratified_2d&,
    NCPA::atmos::ncpaprop_atmosphere_reader_stratified_2d& ) noexcept;

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
                    ::swap( *this, source );
                }

                virtual ~ncpaprop_atmosphere_reader_stratified_2d() {}

                ncpaprop_atmosphere_reader_stratified_2d& operator=(
                    ncpaprop_atmosphere_reader_stratified_2d other ) {
                    ::swap( *this, other );
                    return *this;
                }

                friend void ::swap(
                    ncpaprop_atmosphere_reader_stratified_2d& a,
                    ncpaprop_atmosphere_reader_stratified_2d& b ) noexcept;

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

static void swap(
    NCPA::atmos::ncpaprop_atmosphere_reader_stratified_2d& a,
    NCPA::atmos::ncpaprop_atmosphere_reader_stratified_2d& b ) noexcept {
    ::swap( dynamic_cast<NCPA::atmos::ncpaprop_atmosphere_reader_2d&>( a ),
            dynamic_cast<NCPA::atmos::ncpaprop_atmosphere_reader_2d&>( b ) );
}
