#pragma once
#include "NCPA/constants.hpp"

#include <stdexcept>
#include <vector>

namespace NCPA {
    namespace dsp {
        template<typename T>
        class _window;
        template<typename T>
        class CosineWindow;
        template<typename T>
        class HannWindow;
        template<typename T>
        class HammingWindow;
    }  // namespace dsp
}  // namespace NCPA

template<typename T>
void swap( NCPA::dsp::_window<T>& a, NCPA::dsp::_window<T>& b ) noexcept;
template<typename T>
void swap( NCPA::dsp::CosineWindow<T>& a,
           NCPA::dsp::CosineWindow<T>& b ) noexcept;
template<typename T>
void swap( NCPA::dsp::HannWindow<T>& a, NCPA::dsp::HannWindow<T>& b ) noexcept;
template<typename T>
void swap( NCPA::dsp::HammingWindow<T>& a,
           NCPA::dsp::HammingWindow<T>& b ) noexcept;

namespace NCPA {
    namespace dsp {
        template<typename T>
        class _window : public std::vector<T> {
            public:
                _window() : _window( 0 ) {}

                _window( size_t n, const T& val = 1.0 ) { _init( n, val ); }

                _window( const _window<T>& other ) :
                    std::vector<T> { other } {}

                _window( _window<T>&& other ) noexcept : _window<T>() {
                    ::swap( *this, other );
                }

                virtual ~_window() {}

                friend void ::swap<>( _window<T>& a, _window<T>& b ) noexcept;

                virtual T value( size_t winsize, size_t n ) const = 0;

                virtual const _window& apply( size_t winsize, size_t n,
                                              T& point ) const {
                    point *= this->value( winsize, n );
                    return *this;
                }

                virtual const _window& apply( size_t n, T& point ) const {
                    if (this->empty()) {
                        throw std::out_of_range( "No window size set!" );
                    }
                    if (this->size() <= n) {
                        throw std::out_of_range(
                            "Point " + std::to_string( n )
                            + " out of range for window of size "
                            + std::to_string( this->size() ) );
                    } else {
                        point *= this->at( n );
                        return *this;
                    }
                }

                virtual const _window& apply( std::vector<T>& signal ) const {
                    bool size_matches = ( signal.size() == this->size() );
                    for (size_t n = 0; n < signal.size(); ++n) {
                        if (size_matches) {
                            signal[ n ] *= this->at( n );
                        } else {
                            this->apply( signal.size(), n, signal[ n ] );
                        }
                    }
                    return *this;
                }

                virtual _window& build( size_t npts ) {
                    this->resize( npts, 1.0 );
                    for (size_t i = 0; i < npts; ++i) {
                        this->at( i ) *= this->value( npts, i );
                    }
                    return *this;
                }

            protected:
                void _init( size_t npts, const T& val ) {
                    this->resize( npts, static_cast<T>( val ) );
                }
        };

        template<typename T>
        class CosineWindow : public _window<T> {
            public:
                CosineWindow() : _window<T>(), _a {} {}

                CosineWindow( size_t N, const std::vector<T>& a ) :
                    _window<T>( N ), _a { a } {
                    this->build( N );
                }

                CosineWindow( size_t N, std::initializer_list<T>& a ) :
                    CosineWindow<T>( N, std::vector<T> { a } ) {}

                CosineWindow( const CosineWindow<T>& other ) :
                    _window<T>( other ), _a { other._a } {}

                CosineWindow( CosineWindow<T>&& other ) noexcept :
                    CosineWindow<T>() {
                    ::swap( *this, other );
                }

                virtual ~CosineWindow() {}

                friend void ::swap<>( CosineWindow<T>& a,
                                      CosineWindow<T>& b ) noexcept;

                virtual T value( size_t winsize, size_t n ) const override {
                    T w       = 0;
                    T kfactor = 1.0;
                    for (size_t k = 0; k < _a.size(); ++k) {
                        w += kfactor * _a[ k ]
                           * std::cos( 2.0 * NCPA::constants::PI * k * static_cast<T>( n )
                                       / static_cast<T>( this->size() - 1 ) );
                        kfactor = -kfactor;
                    }
                    return w;
                }

            private:
                std::vector<T> _a;
        };

        template<typename T>
        class HannWindow : public CosineWindow<T> {
            public:
                HannWindow() : CosineWindow<T>() {}

                HannWindow( size_t n ) : CosineWindow<T>( n, { 0.5, 0.5 } ) {}

                HannWindow( const HannWindow<T>& other ) :
                    CosineWindow<T>( other ) {}

                HannWindow( HannWindow<T>&& other ) noexcept :
                    HannWindow<T>() {
                    ::swap( *this, other );
                }

                virtual ~HannWindow() {}

                friend void ::swap<>( HannWindow<T>& a,
                                      HannWindow<T>& b ) noexcept;
        };

        template<typename T>
        using HanningWindow = HannWindow<T>;

        template<typename T>
        class HammingWindow : public CosineWindow<T> {
            public:
                HammingWindow() : CosineWindow<T>() {}

                HammingWindow( size_t n ) :
                    CosineWindow<T>( n, { 0.54, 0.46 } ) {}

                HammingWindow( const HammingWindow<T>& other ) :
                    CosineWindow<T>( other ) {}

                HammingWindow( HammingWindow<T>&& other ) noexcept :
                    HammingWindow<T>() {
                    ::swap( *this, other );
                }

                virtual ~HammingWindow() {}

                friend void ::swap<>( HammingWindow<T>& a,
                                      HammingWindow<T>& b ) noexcept;
        };
    }  // namespace dsp
}  // namespace NCPA

template<typename T>
void swap( NCPA::dsp::_window<T>& a, NCPA::dsp::_window<T>& b ) noexcept {
    using std::swap;
    swap( static_cast<std::vector<T>&>( a ),
          static_cast<std::vector<T>&>( b ) );
}

template<typename T>
void swap( NCPA::dsp::CosineWindow<T>& a,
           NCPA::dsp::CosineWindow<T>& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::dsp::_window<T>&>( a ),
            static_cast<NCPA::dsp::_window<T>&>( b ) );
    swap( a._a, b._a );
}

template<typename T>
void swap( NCPA::dsp::HannWindow<T>& a,
           NCPA::dsp::HannWindow<T>& b ) noexcept {
    ::swap( static_cast<NCPA::dsp::CosineWindow<T>&>( a ),
            static_cast<NCPA::dsp::CosineWindow<T>&>( b ) );
}

template<typename T>
void swap( NCPA::dsp::HammingWindow<T>& a,
           NCPA::dsp::HammingWindow<T>& b ) noexcept {
    ::swap( static_cast<NCPA::dsp::CosineWindow<T>&>( a ),
            static_cast<NCPA::dsp::CosineWindow<T>&>( b ) );
}
