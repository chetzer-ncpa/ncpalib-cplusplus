#pragma once

#pragma once

#include "gsl/gsl_spline.h"
#include "NCPA/arrays.hpp"
#include "NCPA/cloneable.hpp"
#include "NCPA/interpolation/abstract_spline.hpp"
#include "NCPA/interpolation/declarations.hpp"
#include "NCPA/interpolation/exceptions.hpp"
#include "NCPA/interpolation/gsl/functions.hpp"
#include "NCPA/math.hpp"

#include <vector>

namespace NCPA {
    namespace interpolation {
        namespace GSL {
            template<typename FloatType>
            class GSLSpline
                : public Cloneable<GSLSpline<FloatType>, spline1d<FloatType>> {
                    typedef spline1d<FloatType> clone_t;
                    typedef Cloneable<GSLSpline<FloatType>,
                                      spline1d<FloatType>>
                        parent_t;

                public:
                    using parent_t::eval_dddf;
                    using parent_t::eval_ddf;
                    using parent_t::eval_df;
                    using parent_t::eval_f;
                    using parent_t::init;
                    GSLSpline() = default;

                    GSLSpline( const gsl_interp_type *gsl_type ) :
                        parent_t(), _interp_type { gsl_type } {}

                    GSLSpline( const GSLSpline& other )     = default;
                    GSLSpline( GSLSpline&& other ) noexcept = default;

                    virtual ~GSLSpline() { this->clear(); }

                    GSLSpline& operator=( GSLSpline other ) {
                        swap( *this, other );
                        return *this;
                    }

                    friend void swap( GSLSpline<FloatType>& a,
                                      GSLSpline<FloatType>& b ) noexcept {
                        using std::swap;
                        swap( static_cast<parent_t&>( a ),
                              static_cast<parent_t&>( b ) );
                        swap( a._interp_type, b._interp_type );
                        swap( a._spline, b._spline );
                        swap( a._accel, b._accel );
                        swap( a._allocated, b._allocated );
                        // swap( a._filled, b._filled );
                        swap( a._ready, b._ready );
                        // swap( a._x, b._x );
                        // swap( a._f, b._f );
                    }

                    void clear() override {
                        _free_spline();
                        _allocated = false;
                        _ready     = false;
                        clone_t::clear();
                        // static_cast<parent_t*>(this)->clear();
                    }

                    FloatType eval_f(
                        const std::vector<FloatType>& x ) override {
                        this->_check_dimensions( x );
                        this->ready();
                        return static_cast<FloatType>(
                            gsl_spline_eval( _spline, x.front(), _accel ) );
                    }

                    FloatType eval_df( const std::vector<FloatType>& x,
                                       size_t wrt ) override {
                        this->_check_wrt( wrt );
                        this->_check_dimensions( x );
                        this->ready();
                        return static_cast<FloatType>( gsl_spline_eval_deriv(
                            _spline, x.front(), _accel ) );
                    }

                    double eval_ddf( const std::vector<FloatType>& x,
                                     size_t wrt1, size_t wrt2 ) override {
                        this->_check_wrt( wrt1 );
                        this->_check_wrt( wrt2 );
                        this->_check_dimensions( x );
                        this->ready();
                        return static_cast<FloatType>( gsl_spline_eval_deriv2(
                            _spline, x.front(), _accel ) );
                    }

                    double eval_dddf( const std::vector<double>& x,
                                      size_t wrt1, size_t wrt2,
                                      size_t wrt3 ) override {
                        throw UnsupportedDerivative(
                            "GSL splines do not support a third derivative" );
                    }

                    // void fill( const NCPA::arrays::SmartVector<FloatType>&
                    // x,
                    //            const NCPA::arrays::SmartVector<FloatType>& f
                    //            )
                    //     override {
                    //     if (x.rank() != this->dimensions()
                    //         || f.rank() != this->dimensions()) {
                    //         throw DimensionsDoNotMatch(
                    //             "fill() expected vectors of rank 1" );
                    //     }
                    //     this->fill( x.flatten(), f.flatten() );
                    // }

                    // void fill( const std::vector<FloatType>& x,
                    //            const std::vector<FloatType>& f ) override {
                    //     _fill( x, f );
                    // }

                    // void init( std::vector<size_t> dims ) override {
                    //     if (dims.size() != this->dimensions()) {
                    //         throw DimensionsDoNotMatch(
                    //             "init() expected one value for dimensions"
                    //             );
                    //     }
                    //     this->clear();
                    //     this->axes(0).resize( dims.front() );
                    //     this->f().resize( dims.front() );
                    //     _allocate_spline( dims.front() );
                    // }

                    interpolator_type_t interptype() const override {
                        return interpolator_type( _interp_type );
                    }

                    // std::vector<std::pair<FloatType, FloatType>> limits()
                    //     const override {
                    //     return std::vector<std::pair<FloatType, FloatType>>
                    //     {
                    //         { static_cast<FloatType>( this->axes(0).front()
                    //         ),
                    //          static_cast<FloatType>( this->axes(0).back() )
                    //          }
                    //     };
                    // }

                    void ready() override {
                        if (this->size() == 0) {
                            throw InterpolatorNotReady(
                                "Input data not set." );
                        }
                        if (!_allocated) {
                            _allocate_spline( this->size() );
                        }
                        if (!_ready) {
                            _finalize_spline();
                        }
                    }

                protected:
                    const gsl_interp_type *_interp_type = nullptr;
                    gsl_spline *_spline                 = nullptr;
                    gsl_interp_accel *_accel            = nullptr;
                    bool _allocated                     = false;
                    // bool _filled                        = false;
                    bool _ready                         = false;
                    // std::vector<double> _x;
                    // std::vector<double> _f;

                    /*
                        Workflow:
                        init | allocate(n):
                            set _spline;
                            set _accel;
                        fill(x,f):
                            if !allocated
                                allocate(size(x))
                            fill
                            spline_init

                    */

                    void _allocate_spline( size_t n ) {
                        _check_engine();
                        if (n < _interp_type->min_size) {
                            std::ostringstream oss;
                            oss << "GSL interpolation type "
                                << _interp_type->name << " requires at least "
                                << _interp_type->min_size << " points.";
                            throw BadInput( oss.str() );
                        }
                        _spline    = gsl_spline_alloc( _interp_type, n );
                        _accel     = gsl_interp_accel_alloc();
                        _allocated = true;
                    }

                    void _check_accel() const {
                        if (_accel == nullptr) {
                            throw std::logic_error(
                                "Spline accelerator has not been "
                                "initialized" );
                        }
                    }

                    void _check_engine() const {
                        if (_interp_type == nullptr) {
                            throw std::logic_error(
                                "GSL Spline type has not been set!" );
                        }
                    }

                    void _check_input_sizes() const {
                        if (this->axes( 0 ).empty() || this->f().empty()) {
                            throw InterpolatorNotReady(
                                "X and/or F vectors empty" );
                        }
                        if (this->axes( 0 ).size() != this->f().size()) {
                            throw InterpolatorNotReady(
                                "X and F vectors different sizes" );
                        }
                    }

                    // void _check_ready() {
                    //     if (!_ready) {
                    //         _finalize_spline();
                    //     }
                    // }

                    void _check_spline() const {
                        if (_spline == nullptr) {
                            throw std::logic_error(
                                "Spline has not been initialized!" );
                        }
                    }

                    // template<
                    //     typename T                         = FloatType,
                    //     typename std::enable_if<std::is_same<T,
                    //     double>::value,
                    //                             int>::type = 0>
                    // void _fill( const std::vector<FloatType>& x,
                    //             const std::vector<FloatType>& f ) {
                    //     if (x.size() != f.size()) {
                    //         throw InputSizesDoNotMatch(
                    //             "fill() expects input vectors of equal size"
                    //             );
                    //     }
                    //     if (!_allocated) {
                    //         _allocate_spline( x.size() );
                    //     }
                    //     _x      = x;
                    //     _f      = f;
                    //     _filled = true;
                    // }

                    // template<typename T = FloatType,
                    //          typename std::enable_if<
                    //              !std::is_same<T, double>::value, int>::type
                    //          = 0>
                    // void _fill( const std::vector<FloatType>& x,
                    //             const std::vector<FloatType>& f ) {
                    //     if (x.size() != f.size()) {
                    //         throw InputSizesDoNotMatch(
                    //             "fill() expects input vectors of equal size"
                    //             );
                    //     }
                    //     if (!_allocated) {
                    //         _allocate_spline( x.size() );
                    //     }
                    //     this->axes(0).clear();
                    //     this->axes(0).resize( x.size() );
                    //     this->f().clear();
                    //     this->f().resize( f.size() );
                    //     for (size_t i = 0; i < x.size(); ++i) {
                    //         _x[ i ] = x[ i ];
                    //         _f[ i ] = f[ i ];
                    //     }
                    //     _filled = true;
                    // }

                    void _finalize_spline() {
                        if (!( _allocated && this->size() > 0 )) {
                            throw InterpolatorNotReady(
                                "Spline has not been populated" );
                        }
                        if (!_ready) {
                            gsl_spline_init( _spline, this->axes( 0 ).data(),
                                             this->f().data(),
                                             this->axes( 0 ).size() );
                            _ready = true;
                        }
                    }

                    void _free_accel() {
                        if (_accel != nullptr) {
                            gsl_interp_accel_free( _accel );
                            _accel = nullptr;
                        }
                        _ready = false;
                    }

                    void _free_engine() { _interp_type = nullptr; }

                    void _free_spline() {
                        if (_spline != nullptr) {
                            gsl_spline_free( _spline );
                            _spline = nullptr;
                        }
                        _free_accel();
                        _ready     = false;
                        _allocated = false;
                    }
            };
        }  // namespace GSL
    }  // namespace interpolation
}  // namespace NCPA
