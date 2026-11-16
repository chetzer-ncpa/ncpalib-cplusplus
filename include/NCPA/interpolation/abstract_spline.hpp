#pragma once

#include "NCPA/arrays.hpp"
#include "NCPA/cloneable.hpp"
#include "NCPA/interpolation/declarations.hpp"
#include "NCPA/interpolation/exceptions.hpp"
#include "NCPA/types.hpp"

#include <vector>

namespace NCPA {
    namespace interpolation {

        template<typename FloatType>
        class abstract_spline : public CloneBase<abstract_spline<FloatType>> {
            public:
                abstract_spline()          = default;
                virtual ~abstract_spline() = default;

                virtual FloatType eval_f( const std::vector<FloatType>& x )
                    = 0;
                virtual FloatType eval_df( const std::vector<FloatType>& x,
                                           size_t wrt )                = 0;
                virtual FloatType eval_ddf( const std::vector<FloatType>& x,
                                            size_t wrt1, size_t wrt2 ) = 0;
                virtual FloatType eval_dddf( const std::vector<FloatType>& x,
                                             size_t wrt1, size_t wrt2,
                                             size_t wrt3 )             = 0;
                virtual interpolator_type_t interptype() const         = 0;
                virtual void ready()                                   = 0;
                virtual void set_f_value( std::vector<size_t> offsets,
                                          FloatType val )              = 0;
                // virtual void set_f_value( size_t offset1, size_t offset2,
                // FloatType val ) = 0; virtual void set_f_value( size_t
                // offset1, size_t offset2, size_t offset3, FloatType val ) =
                // 0;


                // these are implemented downstream
                virtual void clear()                          = 0;
                virtual size_t dimensions() const             = 0;
                virtual void init( std::vector<size_t> dims ) = 0;
                virtual std::vector<std::pair<FloatType, FloatType>> limits()
                    const                                    = 0;
                virtual void set_axis_value( size_t axnum, size_t index,
                                             FloatType val ) = 0;
                virtual size_t size() const                  = 0;

                template<typename... Coords>
                FloatType eval_f( Coords... coords ) {
                    return this->eval_f( std::vector<FloatType> {
                        static_cast<FloatType>( coords )... } );
                }

                template<
                    typename... Coords,
                    typename = typename std::enable_if<
                        NCPA::types::all_arithmetic<Coords...>::value>::type>
                FloatType eval_df( Coords... coords ) {
                    std::vector<FloatType> c { static_cast<FloatType>(
                        coords )... };
                    if (c.size() == this->dimensions() + 1) {
                        size_t wrt = static_cast<size_t>( c.back() + 0.1 );
                        c.pop_back();
                        return this->eval_df( c, wrt );
                    } else if (c.size() == this->dimensions()
                               && this->dimensions() == 1) {
                        return this->eval_df( c, 0 );
                    } else {
                        throw _ambiguous_coordinates( c, 1, "first" );
                    }
                }

                template<
                    typename... Coords,
                    typename = typename std::enable_if<
                        NCPA::types::all_arithmetic<Coords...>::value>::type>
                FloatType eval_ddf( Coords... coords ) {
                    std::vector<FloatType> c { static_cast<FloatType>(
                        coords )... };
                    if (c.size() == this->dimensions() + 2) {
                        size_t wrt2 = static_cast<size_t>( c.back() + 0.1 );
                        c.pop_back();
                        size_t wrt1 = static_cast<size_t>( c.back() + 0.1 );
                        c.pop_back();
                        return this->eval_ddf( c, wrt1, wrt2 );
                    } else if (c.size() == this->dimensions()
                               && this->dimensions() == 1) {
                        return this->eval_ddf( c, 0, 0 );
                    } else {
                        throw _ambiguous_coordinates( c, 2, "second" );
                    }
                }

                template<
                    typename... Coords,
                    typename = typename std::enable_if<
                        NCPA::types::all_arithmetic<Coords...>::value>::type>
                FloatType eval_dddf( Coords... coords ) {
                    std::vector<FloatType> c { static_cast<FloatType>(
                        coords )... };
                    if (c.size() == this->dimensions() + 3) {
                        size_t wrt3 = static_cast<size_t>( c.back() + 0.1 );
                        c.pop_back();
                        size_t wrt2 = static_cast<size_t>( c.back() + 0.1 );
                        c.pop_back();
                        size_t wrt1 = static_cast<size_t>( c.back() + 0.1 );
                        c.pop_back();
                        return this->eval_dddf( c, wrt1, wrt2, wrt3 );
                    } else if (c.size() == this->dimensions()
                               && this->dimensions() == 1) {
                        return this->eval_dddf( c, 0, 0, 0 );
                    } else {
                        throw _ambiguous_coordinates( c, 3, "third" );
                    }
                }

                template<typename Container1, typename Container2>
                void fill( const Container1& x, const Container2& f ) {
                    this->_check_dimensions( 1 );
                    if (x.size() != f.size()) {
                        throw InputSizesDoNotMatch(
                            "x and f vectors different sizes" );
                    }
                    this->init( x.size() );
                    for (size_t i = 0; i < x.size(); ++i) {
                        this->set_axis_value( 0, i, x.at( i ) );
                        this->set_f_value( i, f.at( i ) );
                    }
                }

                template<typename Container1, typename Container2>
                void fill( const Container1& x1, const Container1& x2,
                           const Container2& f ) {
                    this->_check_dimensions( 2 );
                    if (( x1.size() * x2.size() ) != f.size()) {
                        throw InputSizesDoNotMatch(
                            "x and f vectors different sizes" );
                    }
                    this->init( x1.size(), x2.size() );
                    for (size_t i = 0; i < x1.size(); ++i) {
                        this->set_axis_value( 0, i, x1.at( i ) );
                        for (size_t j = 0; j < x2.size(); ++j) {
                            if (i == 0) {
                                this->set_axis_value( 1, j, x2.at( j ) );
                            }
                            this->set_f_value(
                                std::vector<size_t> { i, j },
                                f.at( std::vector<size_t> { i, j } ) );
                        }
                    }
                }

                template<typename Container1, typename Container2>
                void fill( const Container1& x1, const Container1& x2,
                           const Container1& x3, const Container2& f ) {
                    if (( x1.size() * x2.size() * x3.size() ) != f.size()) {
                        throw InputSizesDoNotMatch(
                            "x and f vectors different sizes" );
                    }
                    this->init( x1.size(), x2.size(), x3.size() );
                    for (size_t i = 0; i < x1.size(); ++i) {
                        this->set_axis_value( 0, i, x1.at( i ) );
                        for (size_t j = 0; j < x2.size(); ++j) {
                            if (i == 0) {
                                this->set_axis_value( 1, j, x2.at( j ) );
                            }
                            for (size_t k = 0; k < x3.size(); ++k) {
                                if (i == 0 && j == 0) {
                                    this->set_axis_value( 2, k, x3.at( k ) );
                                }
                                this->set_f_value(
                                    std::vector<size_t> { i, j, k },
                                    f.at( std::vector<size_t> { i, j, k } ) );
                            }
                        }
                    }
                }

                template<typename... Dims>
                void init( Dims... dims ) {
                    this->init( std::vector<size_t> {
                        static_cast<size_t>( dims )... } );
                }

                virtual void set_f_value( size_t offset, FloatType val ) {
                    this->set_f_value( std::vector<size_t> { offset }, val );
                }

                virtual void set_f_value( size_t o1, size_t o2,
                                          FloatType val ) {
                    this->set_f_value( std::vector<size_t> { o1, o2 }, val );
                }

                virtual void set_f_value( size_t o1, size_t o2, size_t o3,
                                          FloatType val ) {
                    this->set_f_value( std::vector<size_t> { o1, o2, o3 },
                                       val );
                }

            protected:
                AmbiguousCoordinates _ambiguous_coordinates(
                    std::vector<FloatType> c, size_t derivnum,
                    std::string derivname ) const {
                    std::ostringstream oss;
                    oss << "Expected ";
                    if (this->dimensions() == 1) {
                        oss << this->dimensions() << " or ";
                    }
                    oss << this->dimensions() + derivnum << " arguments for "
                        << derivname << " derivative: ("
                        << _num_list( "x", this->dimensions() ) << ","
                        << _num_list( "wrt", 1 ) << ")";

                    if (this->dimensions() == 1) {
                        oss << " or (" << _num_list( "x", this->dimensions() )
                            << ") assuming wrt1 == 0";
                    }
                    oss << ".";
                    return AmbiguousCoordinates( oss.str() );
                }

                template<typename T>
                void _check_dimensions( const std::vector<T>& x ) const {
                    if (x.size() != this->dimensions()) {
                        throw DimensionsDoNotMatch(
                            "Expected " + std::to_string( this->dimensions() )
                            + " coordinates, got "
                            + std::to_string( x.size() ) );
                    }
                }

                void _check_dimensions( size_t n ) const {
                    if (n != this->dimensions()) {
                        throw DimensionsDoNotMatch(
                            "Dimension " + std::to_string( n )
                            + " invalid, only have "
                            + std::to_string( this->dimensions() )
                            + " dimensions." );
                    }
                }

                template<typename T>
                void _check_dimensions(
                    const NCPA::arrays::SmartVector<T>& x ) const {
                    if (x.rank() != this->dimensions()) {
                        throw DimensionsDoNotMatch(
                            "Expected " + std::to_string( this->dimensions() )
                            + " dimensions, got "
                            + std::to_string( x.rank() ) );
                    }
                }

                void _check_wrt( size_t n ) const {
                    if (n >= this->dimensions()) {
                        throw InvalidDimensions(
                            "Cannot differentiate by invalid dimension "
                            + std::to_string( n ) );
                    }
                }

                std::string _num_list( std::string base, size_t n ) const {
                    std::ostringstream oss;
                    for (size_t i = 0; i < n; ++i) {
                        if (i > 0) {
                            oss << ",";
                        }
                        oss << base << i + 1;
                    }
                    return oss.str();
                }
        };

        template<typename FloatType, typename AxesType, typename FType>
        class typed_spline : public abstract_spline<FloatType> {
            public:
                typed_spline()          = default;
                virtual ~typed_spline() = default;
                typed_spline(
                    const typed_spline<FloatType, AxesType, FType>& other )
                    = default;

                    

                AxesType& axes( size_t n ) {
                    if (n >= this->dimensions()) {
                        throw InvalidDimensions(
                            "Dimensions " + std::to_string( n )
                            + " too large for "
                            + std::to_string( this->dimensions() )
                            + "-D interpolator" );
                    }
                    this->init_axes();
                    return _axes.at( n );
                }

                const AxesType& axes( size_t n ) const {
                    if (n >= this->dimensions()) {
                        throw InvalidDimensions(
                            "Dimensions " + std::to_string( n )
                            + " too large for "
                            + std::to_string( this->dimensions() )
                            + "-D interpolator" );
                    }
                    this->init_axes();
                    return _axes.at( n );
                }

                virtual void clear() override {
                    for (auto& ax : _axes) {
                        ax.clear();
                    }
                    _f_vals.clear();
                }

                FType& f() { return _f_vals; }

                const FType& f() const { return _f_vals; }

                virtual void init_axes() {
                    if (_axes.size() < this->dimensions()) {
                        _axes.resize( this->dimensions() );
                    }
                }

                void set_axis_value( size_t axnum, size_t index,
                                     FloatType val ) override {
                    _axes.at( axnum ).at( index ) = val;
                }

                // void set_f_value( std::vector<size_t> offsets,
                //                   FloatType val ) override {
                //     if (offsets.size() == 1) {
                //         this->set_f_value( offsets.front(), val );
                //     } else {
                //         this->_check_dimensions
                //     }
                // }

                // void set_f_value( size_t offset, FloatType val ) override {
                //     this->_check_dimensions( 1 );
                //     _f_vals.at( offset ) = val;
                // }

                // void set_f_value( size_t offset1, size_t offset2,
                //                   FloatType val ) override {
                //     this->_check_dimensions( 2 );
                //     _f_vals.at( offset1, offset2 ) = val;
                // }

                // void set_f_value( size_t offset1, size_t offset2,
                //                   size_t offset3, FloatType val ) override {
                //     this->_check_dimensions( 3 );
                //     _f_vals.at( offset1, offset2, offset3 ) = val;
                // }

                // template<typename InputIterator, typename OutputIterator>
                // void fill( InputIterator x1_begin, InputIterator x1_end,
                // OutputIterator out_it ) {
                //     for (auto it = x1_begin; it != x1_end; ++it) {
                //         *(out_it++) = *it;
                //     }
                // }

                // template<typename InputIterator, typename OutputIterator>
                // void fill( InputIterator x1_begin, InputIterator x1_end,
                //     InputIterator x2_begin, InputIterator x2_end,
                //     OutputIterator out_it ) { for (auto it1 = x1_begin; it1
                //     != x1_end; ++it1) {
                //         for (auto it2 = x2_begin; it2 != x2.end(); ++it2) {
                //             *(out_it++) = *it;
                //         }
                //     }
                // }


                // virtual void fill( const std::vector<FloatType>& x,
                //                    const std::vector<FloatType>& f ) {
                //     throw InvalidDimensions(
                //         "1-D input not valid for "
                //         + std::to_string( this->dimensions() )
                //         + "-D interpolator" );
                // }

                // virtual void fill(
                //     const NCPA::arrays::SmartVector<FloatType>& x,
                //     const NCPA::arrays::SmartVector<FloatType>& f ) {
                //     throw InvalidDimensions(
                //         "1-D input not valid for "
                //         + std::to_string( this->dimensions() )
                //         + "-D interpolator" );
                // }

                std::vector<std::pair<FloatType, FloatType>> limits()
                    const override {
                    std::vector<std::pair<FloatType, FloatType>> lim;
                    for (auto ax : _axes) {
                        lim.push_back(
                            std::make_pair( ax.front(), ax.back() ) );
                    }
                    return lim;
                }

                // // new methods to override as needed
                // virtual void set_basis( size_t x1, FloatType f1 ) {
                //     throw InvalidDimensions(
                //         "1-D input not valid for "
                //         + std::to_string( this->dimensions() )
                //         + "-D interpolator" );
                // }

                // virtual void set_basis( size_t x1, size_t x2, FloatType f1 )
                // {
                //     throw InvalidDimensions(
                //         "2-D input not valid for "
                //         + std::to_string( this->dimensions() )
                //         + "-D interpolator" );
                // }

                // virtual void set_basis( size_t x1, size_t x2, size_t x3,
                //                         FloatType f1 ) {
                //     throw InvalidDimensions(
                //         "3-D input not valid for "
                //         + std::to_string( this->dimensions() )
                //         + "-D interpolator" );
                // }

                // template<typename InputIterator, typename OutputIterator>
                // void fill( InputIterator x1_begin, InputIterator x1_end,
                //            OutputIterator out_it ) {
                //     this->init( std::distance( x1_begin, x1_end ) );
                //     for (auto it = x1_begin; it != x1_end; ++it) {
                //         *( out_it++ ) = *it;
                //     }
                // }

                // template<typename InputIterator, typename OutputIterator>
                // void fill( InputIterator x1_begin, InputIterator x1_end,
                //            InputIterator x2_begin, InputIterator x2_end,
                //            OutputIterator out_it ) {
                //     for (auto it1 = x1_begin; it1 != x1_end; ++it1) {
                //         for () *( out_it++ ) = *it;
                //     }
                // }

                // template<typename... Coords>
                // FloatType eval_f( Coords... coords ) {
                //     return this->eval_f( std::vector<FloatType> {
                //         static_cast<FloatType>( coords )... } );
                // }

                // template<
                //     typename... Coords,
                //     typename = typename std::enable_if<
                //         NCPA::types::all_arithmetic<Coords...>::value>::type>
                // FloatType eval_df( Coords... coords ) {
                //     std::vector<FloatType> c { static_cast<FloatType>(
                //         coords )... };
                //     if (c.size() == this->dimensions() + 1) {
                //         size_t wrt = static_cast<size_t>( c.back() + 0.1 );
                //         c.pop_back();
                //         return this->eval_df( c, wrt );
                //     } else if (c.size() == this->dimensions()
                //                && this->dimensions() == 1) {
                //         return this->eval_df( c, 0 );
                //     } else {
                //         throw _ambiguous_coordinates( c, 1, "first" );
                //     }
                // }

                // template<
                //     typename... Coords,
                //     typename = typename std::enable_if<
                //         NCPA::types::all_arithmetic<Coords...>::value>::type>
                // FloatType eval_ddf( Coords... coords ) {
                //     std::vector<FloatType> c { static_cast<FloatType>(
                //         coords )... };
                //     if (c.size() == this->dimensions() + 2) {
                //         size_t wrt2 = static_cast<size_t>( c.back() + 0.1 );
                //         c.pop_back();
                //         size_t wrt1 = static_cast<size_t>( c.back() + 0.1 );
                //         c.pop_back();
                //         return this->eval_ddf( c, wrt1, wrt2 );
                //     } else if (c.size() == this->dimensions()
                //                && this->dimensions() == 1) {
                //         return this->eval_ddf( c, 0, 0 );
                //     } else {
                //         throw _ambiguous_coordinates( c, 2, "second" );
                //     }
                // }

                // template<
                //     typename... Coords,
                //     typename = typename std::enable_if<
                //         NCPA::types::all_arithmetic<Coords...>::value>::type>
                // FloatType eval_dddf( Coords... coords ) {
                //     std::vector<FloatType> c { static_cast<FloatType>(
                //         coords )... };
                //     if (c.size() == this->dimensions() + 3) {
                //         size_t wrt3 = static_cast<size_t>( c.back() + 0.1 );
                //         c.pop_back();
                //         size_t wrt2 = static_cast<size_t>( c.back() + 0.1 );
                //         c.pop_back();
                //         size_t wrt1 = static_cast<size_t>( c.back() + 0.1 );
                //         c.pop_back();
                //         return this->eval_dddf( c, wrt1, wrt2, wrt3 );
                //     } else if (c.size() == this->dimensions()
                //                && this->dimensions() == 1) {
                //         return this->eval_dddf( c, 0, 0, 0 );
                //     } else {
                //         throw _ambiguous_coordinates( c, 3, "third" );
                //     }
                // }


                size_t size() const override {
                    size_t n = 1;
                    for (const auto& ax : _axes) {
                        n *= ax.size();
                    }
                    return n;
                }


            protected:
                std::vector<AxesType> _axes;
                FType _f_vals;
        };

        template<typename FloatType>
        class spline1d : public typed_spline<FloatType, std::vector<FloatType>,
                                             std::vector<FloatType>> {
            public:
                spline1d()          = default;
                virtual ~spline1d() = default;

                size_t dimensions() const override { return 1; }

                void init( std::vector<size_t> dims ) override {
                    this->_check_dimensions( dims );
                    this->axes(0).resize(dims.front());
                    
                    this->f().resize( dims.front() );
                }

                 
                

                void set_f_value( std::vector<size_t> offsets,
                                  FloatType val ) override {
                    this->_check_dimensions( offsets );
                    this->f().at( offsets.front() ) = val;
                }
        };

        template<typename FloatType>
        class spline2d
            : public typed_spline<FloatType, std::vector<FloatType>,
                                  NCPA::arrays::SmartArray<2, FloatType>> {
            public:
                spline2d()          = default;
                virtual ~spline2d() = default;

                size_t dimensions() const override { return 2; }

                void init( std::vector<size_t> dims ) override {
                    this->_check_dimensions( dims );
                    for (size_t i = 0; i < this->dimensions(); ++i) {
                        this->axes(i).resize( dims.at(i));
                    }
                    this->f().resize( dims );
                }

                void set_f_value( std::vector<size_t> offsets,
                                  FloatType val ) override {
                    this->_check_dimensions( offsets );
                    this->f().at( offsets.at( 0 ), offsets.at( 1 ) ) = val;
                }
        };

        template<typename FloatType>
        class spline3d
            : public typed_spline<FloatType, std::vector<FloatType>,
                                  NCPA::arrays::SmartArray<3, FloatType>> {
            public:
                spline3d()          = default;
                virtual ~spline3d() = default;

                size_t dimensions() const override { return 3; }

                void init( std::vector<size_t> dims ) override {
                    this->_check_dimensions( dims );
                    for (size_t i = 0; i < this->dimensions(); ++i) {
                        this->axes(i).resize( dims.at(i));
                    }
                    this->f().resize( dims );
                }

                void set_f_value( std::vector<size_t> offsets,
                                  FloatType val ) override {
                    this->_check_dimensions( offsets );
                    this->f().at( offsets.at( 0 ), offsets.at( 1 ),
                                  offsets.at( 2 ) ) = val;
                }
        };
    }  // namespace interpolation
}  // namespace NCPA
