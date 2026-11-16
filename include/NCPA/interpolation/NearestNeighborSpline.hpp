#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/cloneable.hpp"
#include "NCPA/interpolation/abstract_spline.hpp"
#include "NCPA/interpolation/declarations.hpp"
#include "NCPA/interpolation/exceptions.hpp"
#include "NCPA/math.hpp"

namespace NCPA {
    namespace interpolation {
        template<typename FloatType>
        class NearestNeighborSpline
            : public Cloneable<
                  NearestNeighborSpline<FloatType>,
                  spline1d<FloatType>> {
                typedef Cloneable<
                  NearestNeighborSpline<FloatType>,
                  spline1d<FloatType>>
                    parent_t;

            public:
                using parent_t::eval_dddf;
                using parent_t::eval_ddf;
                using parent_t::eval_df;
                using parent_t::eval_f;
                using parent_t::init;

                NearestNeighborSpline() = default;
                NearestNeighborSpline(
                    const NearestNeighborSpline<FloatType>& other ) = default;
                NearestNeighborSpline(
                    NearestNeighborSpline<FloatType>&& other ) noexcept
                    = default;
                virtual ~NearestNeighborSpline() = default;

                NearestNeighborSpline<FloatType>& operator=(
                    NearestNeighborSpline<FloatType> other ) {
                    swap( *this, other );
                    return *this;
                }

                void clear() override { parent_t::clear(); }

                FloatType eval_f( const std::vector<FloatType>& x ) override {
                    this->_check_dimensions( x );
                    return this->f().at( NCPA::math::find_closest_index(
                        this->axes(0).begin(), this->axes(0).end(),
                        x.front() ) );
                }

                FloatType eval_df( const std::vector<FloatType>& x,
                                   size_t wrt ) override {
                    this->_check_dimensions( x );
                    this->_check_wrt( wrt );
                    auto upper
                        = std::lower_bound( this->axes(0).begin(),
                                            this->axes(0).end(), x.front() );
                    if (upper == this->axes(0).begin()
                        || upper == this->axes(0).end()) {
                        return 0.0;
                    }
                    auto lower = std::prev( upper );
                    size_t ind1
                        = std::distance( this->axes(0).begin(), upper );
                    size_t ind0 = ind1 - 1;
                    return ( this->f().at( ind1 ) - this->f().at( ind0 ) )
                         / ( *upper - *lower );
                }

                FloatType eval_ddf( const std::vector<FloatType>& x,
                                    size_t wrt1, size_t wrt2 ) override {
                    this->_check_dimensions( x );
                    this->_check_wrt( wrt1 );
                    this->_check_wrt( wrt2 );
                    return 0.0;
                }

                FloatType eval_dddf( const std::vector<FloatType>& x,
                                     size_t wrt1, size_t wrt2,
                                     size_t wrt3 ) override {
                    this->_check_dimensions( x );
                    this->_check_wrt( wrt1 );
                    this->_check_wrt( wrt2 );
                    this->_check_wrt( wrt3 );
                    return 0.0;
                }

                // void fill(
                //     const NCPA::arrays::SmartVector<FloatType>& x,
                //     const NCPA::arrays::SmartVector<FloatType>& f ) override
                //     { if (x.rank() != this->dimensions()
                //         || f.rank() != this->dimensions()) {
                //         throw DimensionsDoNotMatch(
                //             "fill() expected vectors of rank 1" );
                //     }
                //     if (x.size() != f.size()) {
                //         throw InputSizesDoNotMatch(
                //             "fill() expects input vectors of equal size" );
                //     }
                //     _axes.at(0). = x;
                //     _f = f;
                // }

                // void fill( const std::vector<FloatType>& x,
                //            const std::vector<FloatType>& f ) override {
                //     if (x.size() != f.size()) {
                //         throw InputSizesDoNotMatch(
                //             "fill() expects input vectors of equal size" );
                //     }
                //     _axes.at(0).resize( x.size() );
                //     _f_vals.resize( f.size() );
                //     _axes.at(0).set( x );
                //     _f_vals.set( f );
                // }

                // void init( std::vector<size_t> dims ) override {
                //     if (dims.size() != this->dimensions()) {
                //         throw DimensionsDoNotMatch(
                //             "init() expected one value for dimensions" );
                //     }
                //     this->clear();
                //     _axes.at(0).resize( dims.front() );
                //     _f_vals.resize( dims );
                // }

                interpolator_type_t interptype() const override {
                    return interpolator_type_t::NEAREST_NEIGHBOR;
                }

                // std::vector<std::pair<FloatType, FloatType>> limits()
                //     const override {
                //     return std::vector<std::pair<FloatType, FloatType>> {
                //         { _axes.at(0).front(), _axes.at(0).back() }
                //     };
                // }

                void ready() override {
                    if (this->axes(0).size() == 0 || this->f().size() == 0) {
                        throw InterpolatorNotReady(
                            "Interpolator has not been set up" );
                    }
                }

                // virtual const NCPA::arrays::SmartVector<double>
                // source_axes.at(0).()
                //     const = 0;
                // virtual const NCPA::arrays::SmartVector<double> source_f()
                //     const = 0;

            // protected:
                // NCPA::arrays::SmartVector<FloatType> _axes.at(0).;
                // NCPA::arrays::SmartVector<FloatType> _f;
        };
    }  // namespace interpolation
}  // namespace NCPA
