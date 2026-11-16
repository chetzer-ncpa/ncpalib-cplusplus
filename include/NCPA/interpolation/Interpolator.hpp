#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/interpolation/abstract_spline.hpp"
#include "NCPA/interpolation/builders.hpp"
#include "NCPA/interpolation/declarations.hpp"
#include "NCPA/interpolation/exceptions.hpp"
#include "NCPA/types.hpp"

#include <memory>
#include <vector>

namespace NCPA {
    namespace interpolation {
        template<typename FloatType>
        class Interpolator<
            FloatType, typename std::enable_if<
                           std::is_floating_point<FloatType>::value>::type> {
            friend class Interpolator<std::complex<FloatType>>;

            public:
                Interpolator() = default;

                Interpolator( interpolator_type_t interptype ) {
                    _spline
                        = InterpolatorFactory::build<FloatType>( interptype );
                }

                Interpolator(
                    std::unique_ptr<abstract_spline<FloatType>> spline ) {
                    _spline = std::move( spline );
                }

                Interpolator(
                    const std::unique_ptr<abstract_spline<FloatType>>&
                        spline ) {
                    _spline = spline->clone();
                }

                Interpolator(
                    const abstract_spline<FloatType>& spline ) {
                    _spline = spline.clone();
                }

                Interpolator( const Interpolator<FloatType>& other ) {
                    _spline = other._spline->clone();
                }

                Interpolator( Interpolator<FloatType>&& other ) noexcept
                    = default;

                Interpolator<FloatType>& operator=(
                    Interpolator<FloatType> other ) {
                    swap( *this, other );
                    return *this;
                }

                friend void swap( Interpolator<FloatType>& a,
                                  Interpolator<FloatType>& b ) noexcept {
                    using std::swap;
                    swap( a._spline, b._spline );
                }

                void clear() {
                    if (_spline) {
                        _spline->clear();
                    }
                }

                size_t dimensions() const {
                    if (_spline) {
                        return _spline->dimensions();
                    } else {
                        return 0;
                    }
                }

                // FloatType eval_f( const std::vector<FloatType>& x ) {
                //     _check_spline();
                //     return _spline->eval_f( x );
                // }

                // FloatType eval_df( const std::vector<FloatType>& x,
                //                    size_t wrt ) {
                //     _check_spline();
                //     return _spline->eval_df( x, wrt );
                // }

                // FloatType eval_ddf( const std::vector<FloatType>& x,
                //                     size_t wrt1, size_t wrt2 ) {
                //     _check_spline();
                //     return _spline->eval_ddf( x, wrt1, wrt2 );
                // }

                // FloatType eval_dddf( const std::vector<FloatType>& x,
                //                      size_t wrt1, size_t wrt2, size_t wrt3 ) {
                //     _check_spline();
                //     return _spline->eval_dddf( x, wrt1, wrt2, wrt3 );
                // }

                void fill( const std::vector<FloatType>& x,
                           const std::vector<FloatType>& f ) {
                    _check_spline();
                    _spline->fill( x, f );
                }

                void fill( const NCPA::arrays::SmartVector<FloatType>& x,
                           const NCPA::arrays::SmartVector<FloatType>& f ) {
                    _check_spline();
                    _spline->fill( x, f );
                }

                void init( std::vector<size_t> dims ) {
                    _check_spline();
                    _spline->init(dims);
                }

                interpolator_type_t interptype() const {
                    _check_spline();
                    return _spline->interptype();
                }

                std::vector<std::pair<FloatType, FloatType>> limits() const {
                    _check_spline();
                    return _spline->limits();
                }

                void ready() {
                    _check_spline();
                    _spline->ready();
                }

                template<typename... Coords>
                FloatType eval_f( Coords... coords ) {
                    _check_spline();
                    return _spline->eval_f( coords... );
                }

                template<typename... Coords>
                FloatType eval_df( Coords... coords ) {
                    _check_spline();
                    return _spline->eval_df( coords... );
                }

                template<typename... Coords>
                FloatType eval_ddf( Coords... coords) {
                    _check_spline();
                    return _spline->eval_ddf( coords... );
                }

                template<typename... Coords>
                FloatType eval_dddf( Coords... coords ) {
                    _check_spline();
                    return _spline->eval_dddf( coords... );
                }

                template<typename... Dims>
                void init( Dims... dims ) {
                    this->init( std::vector<size_t> {
                        static_cast<size_t>( dims )... } );
                }


            protected:
                std::unique_ptr<abstract_spline<FloatType>> _spline;

                void _check_spline() const {
                    if (!_spline) {
                        throw InterpolatorNotReady();
                    }
                }
        };

        template<typename ComplexType>
        class Interpolator<ComplexType,
                           typename std::enable_if<NCPA::types::is_complex<
                               ComplexType>::value>::type> {
                using FloatType = typename ComplexType::value_type;

            public:
                Interpolator() = default;

                Interpolator( interpolator_type_t interptype ) {
                    _rspline
                        = InterpolatorFactory::build<FloatType>( interptype );
                    _ispline
                        = InterpolatorFactory::build<FloatType>( interptype );
                }

                Interpolator(
                    std::unique_ptr<abstract_spline<FloatType>> rspline,
                    std::unique_ptr<abstract_spline<FloatType>> ispline ) {
                    _rspline = std::move( rspline );
                    _ispline = std::move( ispline );
                }

                Interpolator(
                    const std::unique_ptr<abstract_spline<FloatType>>& rspline,
                    const std::unique_ptr<abstract_spline<FloatType>>&
                        ispline ) {
                    _rspline = rspline->clone();
                    _ispline = ispline->clone();
                }

                Interpolator( const Interpolator<ComplexType>& other ) {
                    _rspline = other._rspline->clone();
                    _ispline = other._ispline->clone();
                }

                Interpolator( const Interpolator<FloatType>& other ) {
                    _rspline = other._spline->clone();
                    // _ispline = other._ispline->clone();
                }

                Interpolator( const Interpolator<FloatType>& other_r,
                const Interpolator<FloatType>& other_i ) {
                    _rspline = other_r._spline->clone();
                    _ispline = other_i._spline->clone();
                    // _ispline = other._ispline->clone();
                }

                Interpolator( Interpolator<ComplexType>&& other ) noexcept
                    = default;

                Interpolator<ComplexType>& operator=(
                    Interpolator<ComplexType> other ) {
                    swap( *this, other );
                    return *this;
                }

                friend void swap( Interpolator<ComplexType>& a,
                                  Interpolator<ComplexType>& b ) noexcept {
                    using std::swap;
                    swap( a._rspline, b._rspline );
                    swap( a._ispline, b._ispline );
                }

                void clear() {
                    if (_rspline) {
                        _rspline->clear();
                    }
                    if (_ispline) {
                        _ispline->clear();
                    }
                }

                size_t dimensions() const {
                    if (_rspline) {
                        return _rspline->dimensions();
                    } else {
                        return 0;
                    }
                }

                // ComplexType eval_f( const std::vector<FloatType>& x ) {
                //     _check_rspline();
                //     return ComplexType { _rspline->eval_f( x ),
                //                          _ispline ? _ispline->eval_f( x ),
                //                          0.0 };
                // }

                // ComplexType eval_df( const std::vector<FloatType>& x,
                //                    size_t wrt ) {
                //     _check_rspline();
                //     return ComplexType { _rspline->eval_df( x, wrt ),
                //                          _ispline
                //                          ? _ispline->eval_df( x, wrt ),
                //                          0.0 };
                // }

                // ComplexType eval_ddf( const std::vector<FloatType>& x,
                //                     size_t wrt1, size_t wrt2 ) {
                //     _check_rspline();
                //     return ComplexType { _rspline->eval_ddf( x, wrt1, wrt2 ),
                //                          _ispline
                //                          ? _ispline->eval_ddf( x, wrt1, wrt2 ),
                //                          0.0 };
                // }

                // ComplexType eval_dddf( const std::vector<FloatType>& x,
                //                      size_t wrt1, size_t wrt2, size_t wrt3 ) {
                //     _check_rspline();
                //     return ComplexType {
                //         _rspline->eval_dddf( x, wrt1, wrt2, wrt3 ),
                //         _ispline ? _ispline->eval_dddf( x, wrt1, wrt2, wrt3 ),
                //         0.0
                //     };
                // }

                void fill( const std::vector<FloatType>& x,
                           const std::vector<FloatType>& fr,
                        const std::vector<FloatType>& fi ) {
                    _check_rspline();
                    _check_ispline();
                    _rspline->fill( x, fr );
                    _ispline->fill( x, fi );
                }

                void fill( const std::vector<FloatType>& x,
                           const std::vector<ComplexType>& f ) {
                    _check_rspline();
                    _check_ispline();
                    std::vector<FloatType> rvec, ivec;
                    rvec.reserve(f.size());
                    ivec.reserve(f.size());
                    for (auto& val : f) {
                        rvec.push_back( f.real() );
                        ivec.push_vack( f.imag() );
                    }
                    _rspline->fill( x, rvec );
                    _ispline->fill( x, ivec );
                }

                void fill( const NCPA::arrays::SmartVector<FloatType>& x,
                           const NCPA::arrays::SmartVector<FloatType>& fr,
                            const NCPA::arrays::SmartVector<FloatType>& fi ) {
                    _check_rspline();
                    _check_ispline();
                    _rspline->fill( x, fr );
                    _ispline->fill( x, fi );
                }

                void fill( const NCPA::arrays::SmartVector<FloatType>& x,
                           const NCPA::arrays::SmartVector<ComplexType>& f ) {
                    _check_rspline();
                    _check_ispline();
                    std::vector<FloatType> rvec, ivec;
                    rvec.reserve(f.size());
                    ivec.reserve(f.size());
                    for (auto& val : f) {
                        rvec.push_back( f.real() );
                        ivec.push_vack( f.imag() );
                    }
                    _rspline->fill( x, rvec );
                    _ispline->fill( x, ivec );
                }

                void init( std::vector<size_t> dims ) {
                    _check_rspline();
                    _rspline->init(dims);
                    if (_ispline) {
                        _ispline->init(dims);
                    }
                }

                interpolator_type_t interptype() const {
                    _check_rspline();
                    return _rspline->interptype();
                }

                std::vector<std::pair<FloatType, FloatType>> limits() const {
                    _check_rspline();
                    return _rspline->limits();
                }

                void ready() {
                    _check_rspline();
                    _rspline->ready();
                    if (_ispline) {
                        _ispline->ready();
                    }
                }

                template<typename... Coords>
                ComplexType eval_f( Coords... coords ) {
                    _check_rspline();
                    return ComplexType{
                        _rspline->eval_f( coords... ),
                        _ispline ? _ispline->eval_f( coords... ) : 0.0
                    };
                }

                template<typename... Coords>
                ComplexType eval_df( Coords... coords ) {
                    _check_rspline();
                    return ComplexType{
                        _rspline->eval_df( coords... ),
                        _ispline ? _ispline->eval_df( coords... ) : 0.0
                    };
                }

                template<typename... Coords>
                ComplexType eval_ddf( Coords... coords ) {
                    _check_rspline();
                    return ComplexType{
                        _rspline->eval_ddf( coords... ),
                        _ispline ? _ispline->eval_ddf( coords... ) : 0.0
                    };
                }

                template<typename... Coords>
                ComplexType eval_dddf( Coords... coords ) {
                    _check_rspline();
                    return ComplexType{
                        _rspline->eval_dddf( coords... ),
                        _ispline ? _ispline->eval_dddf( coords... ) : 0.0
                    };
                }

                template<typename... Dims>
                void init( Dims... dims ) {
                    this->init( std::vector<size_t> {
                        static_cast<size_t>( dims )... } );
                }


            protected:
                std::unique_ptr<abstract_spline<FloatType>> _rspline, _ispline;

                void _check_rspline() const {
                    if (!_rspline) {
                        throw InterpolatorNotReady();
                    }
                }
                void _check_ispline() const {
                    if (!_ispline) {
                        throw InterpolatorNotReady();
                    }
                }
        };
    }  // namespace interpolation
}  // namespace NCPA
