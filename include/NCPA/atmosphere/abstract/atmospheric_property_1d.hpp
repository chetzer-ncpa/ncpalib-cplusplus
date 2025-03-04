#pragma once

#include "NCPA/atmosphere/abstract/atmospheric_property.hpp"
#include "NCPA/interpolation.hpp"

#define RETURN_THIS_AS_NCPA_ATMOS_ABSTRACT_ATMOSPHERIC_PROPERTY_1D       \
    return static_cast<NCPA::atmos::abstract::atmospheric_property_1d&>( \
        *this )


static void swap( NCPA::atmos::abstract::atmospheric_property_1d& a,
                  NCPA::atmos::abstract::atmospheric_property_1d& b ) noexcept;

namespace NCPA {
    namespace atmos {
        namespace abstract {
            class atmospheric_property_1d : public atmospheric_property {
                public:
                    virtual ~atmospheric_property_1d() {}

                    friend void ::swap( atmospheric_property_1d& a,
                                        atmospheric_property_1d& b ) noexcept;

                    virtual size_t dimensions() const { return 1; }

                    virtual double get( double x1 )                   = 0;
                    virtual double get_first_derivative( double x1 )  = 0;
                    virtual double get_second_derivative( double x1 ) = 0;
                    virtual atmospheric_property_1d& set_interpolator(
                        NCPA::interpolation::interpolator_1d_type_t
                            interp_type )
                        = 0;
                    virtual size_t size() const = 0;
                    virtual atmospheric_property_1d& set( const vector_u_t& z,
                                                          const vector_u_t& p )
                        = 0;
                    virtual vector_u_t& vector()             = 0;
                    virtual const vector_u_t& vector() const = 0;
                    virtual atmospheric_property_1d& resample(
                        vector_u_t new_z )
                        = 0;
                    virtual NCPA::interpolation::Interpolator1D<double,
                                                                double>&
                        spline()
                        = 0;
            };
        }  // namespace abstract
    }  // namespace atmos
}  // namespace NCPA

static void swap(
    NCPA::atmos::abstract::atmospheric_property_1d& a,
    NCPA::atmos::abstract::atmospheric_property_1d& b ) noexcept {
    using std::swap;
    ::swap( dynamic_cast<NCPA::atmos::abstract::atmospheric_property&>( a ),
            dynamic_cast<NCPA::atmos::abstract::atmospheric_property&>( b ) );
    // std::swap( a._spline, b._spline );
    // swap( a._spline, b._spline );
}
