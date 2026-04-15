#pragma once

#include "NCPA/atmosphere/abstract_atmosphere_1d.hpp"
#include "NCPA/atmosphere/abstract_atmosphere_2d.hpp"
#include "NCPA/atmosphere/abstract_atmosphere_3d.hpp"
#include "NCPA/atmosphere/Atmosphere1D.hpp"
#include "NCPA/atmosphere/Atmosphere2D.hpp"
#include "NCPA/atmosphere/Atmosphere3D.hpp"
#include "NCPA/atmosphere/AtmosphericProperty1D.hpp"
#include "NCPA/atmosphere/AtmosphericProperty2D.hpp"
#include "NCPA/atmosphere/AtmosphericProperty3D.hpp"
#include "NCPA/atmosphere/grid_atmosphere_2d.hpp"
#include "NCPA/atmosphere/grid_atmosphere_3d.hpp"
#include "NCPA/atmosphere/piecewise_stratified_atmosphere_2d.hpp"
#include "NCPA/atmosphere/readers.hpp"
#include "NCPA/atmosphere/stratified_atmosphere_2d.hpp"
#include "NCPA/atmosphere/stratified_atmosphere_3d.hpp"
#include "NCPA/strings.hpp"
#include "NCPA/units.hpp"

#include <fstream>
#include <iostream>
#include <string>

namespace NCPA {
    namespace atmos {
        class AtmosphereFactory {
            public:
                static bool can_build( atmospheric_property_1d_t proptype ) {
                    switch (proptype) {
                        case atmospheric_property_1d_t::TUPLE:
                            return true;
                            break;
                        default:
                            return false;
                    }
                }

                static AtmosphericProperty1D build(
                    atmospheric_property_1d_t proptype ) {
                    switch (proptype) {
                        case atmospheric_property_1d_t::TUPLE:
                            return AtmosphericProperty1D( _atm_prop_1d_ptr_t(
                                new tuple_atmospheric_property_1d() ) );
                            break;
                        default:
                            throw std::range_error(
                                "Requested 1-D atmospheric property type "
                                "unrecognized or not yet implemented" );
                    }
                }

                static bool can_build( atmospheric_property_2d_t proptype ) {
                    switch (proptype) {
                        case atmospheric_property_2d_t::STRATIFIED:
                        case atmospheric_property_2d_t::GRID:
                            return true;
                            break;
                        default:
                            return false;
                    }
                }

                static AtmosphericProperty2D build(
                    atmospheric_property_2d_t proptype ) {
                    switch (proptype) {
                        case atmospheric_property_2d_t::STRATIFIED:
                            return AtmosphericProperty2D( _atm_prop_2d_ptr_t(
                                new stratified_atmospheric_property_2d() ) );
                            break;
                        case atmospheric_property_2d_t::GRID:
                            return AtmosphericProperty2D( _atm_prop_2d_ptr_t(
                                new grid_atmospheric_property_2d() ) );
                            break;
                        default:
                            throw std::range_error(
                                "Requested 2-D atmospheric property type "
                                "unrecognized or not yet implemented" );
                    }
                }

                static bool can_build( atmospheric_property_3d_t proptype ) {
                    switch (proptype) {
                        case atmospheric_property_3d_t::STRATIFIED:
                        case atmospheric_property_3d_t::GRID:
                            return true;
                            break;
                        default:
                            return false;
                    }
                }

                static AtmosphericProperty3D build(
                    atmospheric_property_3d_t proptype ) {
                    switch (proptype) {
                        case atmospheric_property_3d_t::STRATIFIED:
                            return AtmosphericProperty3D( _atm_prop_3d_ptr_t(
                                new stratified_atmospheric_property_3d() ) );
                            break;
                        case atmospheric_property_3d_t::GRID:
                            return AtmosphericProperty3D( _atm_prop_3d_ptr_t(
                                new grid_atmospheric_property_3d() ) );
                            break;
                        default:
                            throw std::range_error(
                                "Requested 3-D atmospheric property type "
                                "unrecognized or not yet implemented" );
                    }
                }

                static bool can_build( atmosphere_1d_t proptype ) {
                    switch (proptype) {
                        case atmosphere_1d_t::TUPLE:
                            return true;
                            break;
                        default:
                            return false;
                    }
                }

                static Atmosphere1D build( atmosphere_1d_t atmostype ) {
                    Atmosphere1D atmos;
                    switch (atmostype) {
                        case atmosphere_1d_t::TUPLE:
                            return Atmosphere1D(
                                _atm_1d_ptr_t( new tuple_atmosphere_1d() ) );
                            break;
                        default:
                            throw std::range_error(
                                "Requested 1-D atmosphere type "
                                "unrecognized or not yet implemented" );
                    }
                }

                static bool can_build( atmosphere_2d_t proptype ) {
                    switch (proptype) {
                        case atmosphere_2d_t::GRID:
                        case atmosphere_2d_t::STRATIFIED:
                        case atmosphere_2d_t::PIECEWISE_STRATIFIED:
                            return true;
                            break;
                        default:
                            return false;
                    }
                }

                static Atmosphere2D build( atmosphere_2d_t atmostype ) {
                    Atmosphere1D atmos;
                    switch (atmostype) {
                        case atmosphere_2d_t::STRATIFIED:
                            return Atmosphere2D( _atm_2d_ptr_t(
                                new stratified_atmosphere_2d() ) );
                            break;
                        case atmosphere_2d_t::PIECEWISE_STRATIFIED:
                            return Atmosphere2D( _atm_2d_ptr_t(
                                new piecewise_stratified_atmosphere_2d() ) );
                            break;
                        case atmosphere_2d_t::GRID:
                            return Atmosphere2D(
                                _atm_2d_ptr_t( new grid_atmosphere_2d() ) );
                            break;
                        default:
                            throw std::range_error(
                                "Requested 2-D atmosphere type "
                                "unrecognized or not yet implemented" );
                    }
                }

                static bool can_build( atmosphere_3d_t proptype ) {
                    switch (proptype) {
                        case atmosphere_3d_t::GRID:
                        case atmosphere_3d_t::STRATIFIED:
                            return true;
                            break;
                        default:
                            return false;
                    }
                }

                static Atmosphere3D build( atmosphere_3d_t atmostype ) {
                    switch (atmostype) {
                        case atmosphere_3d_t::STRATIFIED:
                            return Atmosphere3D( _atm_3d_ptr_t(
                                new stratified_atmosphere_3d() ) );
                            break;
                        case atmosphere_3d_t::GRID:
                            return Atmosphere3D(
                                _atm_3d_ptr_t( new grid_atmosphere_3d() ) );
                            break;

                        default:
                            throw std::range_error(
                                "Requested 3-D atmosphere type "
                                "unrecognized or not yet implemented" );
                    }
                }

                static bool can_build( reader_1d_t proptype ) {
                    switch (proptype) {
                        case reader_1d_t::NCPAPROP:
                            return true;
                            break;
                        default:
                            return false;
                    }
                }

                static AtmosphereReader1D build( reader_1d_t atmostype ) {
                    switch (atmostype) {
                        case reader_1d_t::NCPAPROP:
                            return AtmosphereReader1D(
                                std::unique_ptr<
                                    _abstract_atmosphere_reader_1d>(
                                    new ncpaprop_atmosphere_reader_1d() ) );
                            break;
                        default:
                            throw std::range_error(
                                "Requested 1-D atmosphere type "
                                "unrecognized or not yet implemented" );
                    }
                }

                static bool can_build( reader_2d_t proptype ) {
                    switch (proptype) {
                        case reader_2d_t::NCPAPROP:
                            return true;
                            break;
                        default:
                            return false;
                    }
                }

                static AtmosphereReader2D build( reader_2d_t atmostype ) {
                    switch (atmostype) {
                        case reader_2d_t::NCPAPROP:
                            return AtmosphereReader2D(
                                std::unique_ptr<
                                    _abstract_atmosphere_reader_2d>(
                                    new ncpaprop_atmosphere_reader_2d() ) );
                            break;
                        case reader_2d_t::NCPAPROP_PIECEWISE_STRATIFIED:
                            return AtmosphereReader2D(
                                std::unique_ptr<
                                    _abstract_atmosphere_reader_2d>(
                                    new ncpaprop_atmosphere_reader_stratified_2d() ) );
                            break;
                        default:
                            throw std::range_error(
                                "Requested 2-D atmosphere type "
                                "unrecognized or not yet implemented" );
                    }
                }

                static bool can_build( reader_3d_t proptype ) {
                    switch (proptype) {
                        case reader_3d_t::NCPAPROP:
                        case reader_3d_t::NCPAPROP_STRATIFIED:
                            return true;
                            break;
                        default:
                            return false;
                    }
                }

                static AtmosphereReader3D build( reader_3d_t atmostype ) {
                    switch (atmostype) {
                        case reader_3d_t::NCPAPROP:
                            return AtmosphereReader3D(
                                std::unique_ptr<
                                    _abstract_atmosphere_reader_3d>(
                                    new ncpaprop_atmosphere_reader_3d() ) );
                            break;
                        case reader_3d_t::NCPAPROP_STRATIFIED:
                            return AtmosphereReader3D(
                                std::unique_ptr<
                                    _abstract_atmosphere_reader_3d>(
                                    new ncpaprop_atmosphere_reader_stratified_3d() ) );
                            break;
                        default:
                            throw std::range_error(
                                "Requested 3-D atmosphere type "
                                "unrecognized or not yet implemented" );
                    }
                }


            private:
                std::vector<std::string> headerlines;
        };
    }  // namespace atmos
}  // namespace NCPA
