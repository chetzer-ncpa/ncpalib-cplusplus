#pragma once

#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/units.hpp"

#define RETURN_THIS_AS_NCPA_ATMOS_ABSTRACT_ATMOSPHERIC_PROPERTY \
    return static_cast<NCPA::atmos::abstract::atmospheric_property&>( *this )

namespace NCPA {
    namespace atmos {
        namespace abstract {
            class atmospheric_property {
                public:
                    virtual ~atmospheric_property() {}

                    virtual size_t dimensions() const           = 0;
                    virtual const units_ptr_t get_units() const = 0;
                    virtual const units_ptr_t get_independent_units(
                        size_t dimension ) const
                        = 0;
                protected:
                    virtual atmospheric_property& _convert_units(
                        const NCPA::units::Unit& u )
                        = 0;
                    virtual atmospheric_property& _convert_independent_units(
                        size_t dimension, const NCPA::units::Unit& u )
                        = 0;
                    // virtual atmospheric_property& setup_spline() = 0;
                    // virtual double get( const std::vector<double>& coords )
                    //     = 0;
                    // virtual double get_first_derivative(
                    //     const std::vector<double>& coords,
                    //     std::vector<size_t>& rel )
                    //     = 0;
                    // virtual double get_second_derivative(
                    //     const std::vector<double>& coords,
                    //     std::vector<size_t>& rel )
                    //     = 0;

                    // virtual double get_first_derivative(
                    //     const std::vector<double>& coords ) {
                    //     std::vector<size_t> empty;
                    //     return this->get_first_derivative( coords, empty );
                    // }

                    // virtual double get_second_derivative(
                    //     const std::vector<double>& coords ) {
                    //     std::vector<size_t> empty;
                    //     return this->get_second_derivative( coords, empty );
                    // }
            };
        }  // namespace abstract
    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::abstract::atmospheric_property& a,
                  NCPA::atmos::abstract::atmospheric_property& b ) noexcept {}
