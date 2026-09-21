#pragma once

#if __has_include( "fftw3.h" )
#  ifndef NCPA_HAVE_FFTW3
#    define NCPA_HAVE_FFTW3
#  endif
#  include "fftw3.h"
#endif

#if __has_include( "pocketfft_hdronly.h" )
#  ifndef NCPA_HAVE_POCKETFFT
#    define NCPA_HAVE_POCKETFFT
#  endif
#  include "pocketfft_hdronly.h"
#endif

#include "NCPA/cloneable.hpp"
#include "NCPA/constants.hpp"

#include <complex>
#include <memory>
#include <stdexcept>
#include <vector>

namespace NCPA {
    namespace dsp {

        enum class fft_sign_t { NEGATIVE, POSITIVE };
        enum class fft_t { FFTW, FFTW3, POCKETFFT };
        enum class fft_scaling_t { NONE, HALF, FULL };

        class FFT : public virtual CloneBase<FFT> {
            public:
                virtual void compute( std::complex<double> *in,
                                      std::complex<double> *out, size_t NFFT,
                                      fft_sign_t sign,
                                      fft_scaling_t scaling
                                      = fft_scaling_t::NONE ) = 0;

                virtual void compute( double *in, std::complex<double> *out,
                                      size_t NFFT, fft_sign_t sign,
                                      fft_scaling_t scaling
                                      = fft_scaling_t::NONE ) = 0;

            public:
                FFT() {}

                FFT( const std::vector<std::complex<double>>& input ) {
                    this->set( input );
                }

                FFT( const std::vector<double>& input ) { this->set( input ); }

                FFT( const FFT& other ) { _input = other._input; }

                FFT( FFT&& other ) noexcept { swap( *this, other ); }

                virtual ~FFT() {}

                friend void swap( FFT& a, FFT& b ) noexcept {
                    using std::swap;
                    swap( a._input, b._input );
                }

                virtual FFT& set(
                    const std::vector<std::complex<double>>& in ) {
                    _input = in;
                    return *this;
                }

                virtual FFT& set( const std::vector<double>& in ) {
                    _input.resize( in.size() );
                    for (size_t i = 0; i < in.size(); ++i) {
                        _input.at( i ).real( in.at( i ) );
                    }
                    return *this;
                }

                virtual void compute( std::vector<std::complex<double>>& in,
                                      std::vector<std::complex<double>>& out,
                                      fft_sign_t sign,
                                      fft_scaling_t scaling
                                      = fft_scaling_t::NONE ) {
                    out.resize( in.size() );
                    this->compute( in.data(), out.data(), in.size(), sign,
                                   scaling );
                }

                virtual void compute( std::vector<double>& in,
                                      std::vector<std::complex<double>>& out,
                                      fft_sign_t sign,
                                      fft_scaling_t scaling
                                      = fft_scaling_t::NONE ) {
                    out.resize( in.size() );
                    this->compute( in.data(), out.data(), in.size(), sign,
                                   scaling );
                }

                virtual std::vector<std::complex<double>> compute(
                    fft_sign_t sign,
                    fft_scaling_t scaling = fft_scaling_t::NONE ) {
                    std::vector<std::complex<double>> out( _input.size() );
                    this->compute( _input.data(), out.data(), _input.size(),
                                   sign, scaling );
                    return out;
                }

                virtual std::vector<std::complex<double>> compute(
                    std::vector<std::complex<double>>& in, fft_sign_t sign,
                    fft_scaling_t scaling = fft_scaling_t::NONE ) {
                    std::vector<std::complex<double>> out( in.size() );
                    this->compute( in.data(), out.data(), in.size(), sign,
                                   scaling );
                    return out;
                }

                virtual std::vector<std::complex<double>> compute(
                    std::vector<double>& in, fft_sign_t sign,
                    fft_scaling_t scaling = fft_scaling_t::NONE ) {
                    std::vector<std::complex<double>> out( in.size() );
                    this->compute( in.data(), out.data(), in.size(), sign,
                                   scaling );
                    return out;
                }

                static double scaling_factor( fft_scaling_t scaling_type,
                                              size_t NFFT ) {
                    if (scaling_type == fft_scaling_t::NONE) {
                        return 1.0;
                    } else if (scaling_type == fft_scaling_t::FULL) {
                        return 1.0 / (double)NFFT;
                    } else if (scaling_type == fft_scaling_t::HALF) {
                        return 2.0 / (double)NFFT;
                    } else {
                        throw std::out_of_range( "Unrecognized or unsupported "
                                                 "scaling type requested" );
                    }
                }

                static void mirror( std::vector<std::complex<double>>& in ) {
                    mirror( in.data(), in.size() );
                }

                static void mirror( std::complex<double> *in, size_t nfft ) {
                    for (size_t i = 1; i < nfft / 2; ++i) {
                        in[ nfft - i ] = std::conj( in[ i ] );
                    }
                }

            protected:
                std::vector<std::complex<double>> _input;
                // std::vector<std::complex<double>> _output;
        };

        typedef std::unique_ptr<FFT> fft_ptr_t;

#ifdef NCPA_HAVE_POCKETFFT
        class PocketFFT : public FFT,
                          public Cloneable<PocketFFT, FFT> {
            public:
                using FFT::compute;

                PocketFFT() : FFT() {}

                PocketFFT( const std::vector<std::complex<double>>& input ) :
                    FFT( input ) {}

                PocketFFT( const std::vector<double>& input ) : FFT( input ) {}

                PocketFFT( const PocketFFT& other ) : FFT( other ) {}

                PocketFFT( PocketFFT&& other ) noexcept {
                    swap( *this, other );
                }

                virtual ~PocketFFT() {}

                friend void swap( PocketFFT& a, PocketFFT& b ) noexcept {
                    using std::swap;
                    swap( static_cast<FFT&>( a ), static_cast<FFT&>( b ) );
                }

                // NCPA_CLONE_METHOD( PocketFFT, FFT )

                void compute( std::complex<double> *in,
                              std::complex<double> *out, size_t NFFT,
                              fft_sign_t sign,
                              fft_scaling_t scaling
                              = fft_scaling_t::NONE ) override {
                    using namespace pocketfft;
                    bool forward = ( sign == fft_sign_t::NEGATIVE );
                    shape_t shape { NFFT };
                    stride_t stride { sizeof( std::complex<double> ) };
                    shape_t axes { 0 };
                    c2c( shape, stride, stride, axes, forward, in, out,
                         FFT::scaling_factor( scaling, NFFT ), 1 );
                }

                void compute( double *in, std::complex<double> *out,
                              size_t NFFT, fft_sign_t sign,
                              fft_scaling_t scaling
                              = fft_scaling_t::NONE ) override {
                    using namespace pocketfft;
                    bool forward = ( sign == fft_sign_t::NEGATIVE );
                    shape_t shape { NFFT };
                    stride_t stride_in { sizeof( double ) };
                    stride_t stride_out { sizeof( std::complex<double> ) };
                    shape_t axes { 0 };
                    r2c( shape, stride_in, stride_out, axes, forward, in, out,
                         FFT::scaling_factor( scaling, NFFT ), 1 );
                    FFT::mirror( out, NFFT );
                }
        };
#endif

#ifdef NCPA_HAVE_FFTW3
        class FFTW : public FFT,
                     public Cloneable<FFTW, FFT> {
            public:
                using FFT::compute;

                FFTW() : FFT() {}

                FFTW( const std::vector<std::complex<double>>& input ) :
                    FFT( input ) {}

                FFTW( const std::vector<double>& input ) : FFT( input ) {}

                FFTW( const FFTW& other ) : FFT( other ) {}

                FFTW( FFTW&& other ) noexcept { swap( *this, other ); }

                virtual ~FFTW() {}

                friend void swap( FFTW& a, FFTW& b ) noexcept {
                    using std::swap;
                    swap( static_cast<FFT&>( a ), static_cast<FFT&>( b ) );
                }

                void compute( std::complex<double> *in,
                              std::complex<double> *out, size_t NFFT,
                              fft_sign_t sign,
                              fft_scaling_t scaling
                              = fft_scaling_t::NONE ) override {
                    // FORWARD in FFTW is negative sign
                    int fftwsign = ( sign == fft_sign_t::NEGATIVE ? -1 : 1 );
                    fftw_plan p  = fftw_plan_dft_1d(
                        (int)NFFT, reinterpret_cast<fftw_complex *>( in ),
                        reinterpret_cast<fftw_complex *>( out ), fftwsign,
                        FFTW_ESTIMATE );
                    fftw_execute( p );
                    fftw_destroy_plan( p );
                    if (scaling != fft_scaling_t::NONE) {
                        double scale = FFT::scaling_factor( scaling, NFFT );
                        for (size_t i = 0; i < NFFT; ++i) {
                            out[ i ] *= scale;
                        }
                    }
                }

                virtual void compute( double *in, std::complex<double> *out,
                                      size_t NFFT, fft_sign_t sign,
                                      fft_scaling_t scaling
                                      = fft_scaling_t::NONE ) override {
                    fftw_plan p = fftw_plan_dft_r2c_1d(
                        (int)NFFT, in, reinterpret_cast<fftw_complex *>( out ),
                        FFTW_ESTIMATE );
                    fftw_execute( p );

                    // for (size_t i = 1; i < NFFT / 2; ++i) {
                    //     out[ NFFT - i ] = std::conj( out[ i ] );
                    // }
                    fftw_destroy_plan( p );
                    if (scaling != fft_scaling_t::NONE) {
                        double scale = FFT::scaling_factor( scaling, NFFT );
                        for (size_t i = 0; i <= NFFT / 2; ++i) {
                            out[ i ] *= scale;
                        }
                    }
                    FFT::mirror( out, NFFT );
                }
        };
#endif

        class FFTFactory {
            public:
                static fft_ptr_t build( fft_t fft_type ) {
                    switch (fft_type) {
#ifdef NCPA_HAVE_FFTW3
                        case fft_t::FFTW:
                        case fft_t::FFTW3:
                            return fft_ptr_t( new FFTW() );
                            break;
#endif
#ifdef NCPA_HAVE_POCKETFFT
                        case fft_t::POCKETFFT:
                            return fft_ptr_t( new PocketFFT() );
                            break;
#endif
                        default:
                            throw std::out_of_range(
                                "Unsupported or unavailable FFT type "
                                "requested" );
                    }
                }

                static fft_ptr_t build(
                    fft_t fft_type,
                    const std::vector<std::complex<double>>& input ) {
                    fft_ptr_t fft = build( fft_type );
                    fft->set( input );
                    return fft;
                }

                static fft_ptr_t build( fft_t fft_type,
                                        const std::vector<double>& input ) {
                    fft_ptr_t fft = build( fft_type );
                    fft->set( input );
                    return fft;
                }
        };

        template<typename T>
        class _window;

        template<typename T>
        std::vector<T>& operator*=( std::vector<T>& v, const _window<T>& w );

        template<typename T>
        class _window {
            public:
                _window() : _window( 0 ) {}

                _window( size_t n, const T& val = 1.0 ) { _init( n, val ); }

                _window( const _window<T>& other ) { _vals = other._vals; }

                _window( _window<T>&& other ) noexcept : _window<T>() {
                    swap( *this, other );
                }

                virtual ~_window() {}

                friend void swap( _window<T>& a, _window<T>& b ) noexcept {
                    using std::swap;
                    swap( a._vals, b._vals );
                }

                friend std::vector<T> operator*( const std::vector<T>& v,
                                                 const _window<T>& w ) {
                    if (v.size() != w.size()) {
                        throw std::range_error( "Incompatible vector sizes "
                                                "for window multiplication" );
                    }
                    std::vector<T> x  = v;
                    x                *= w;
                    return x;
                }

                friend std::vector<T> operator*( const _window<T>& w,
                                                 const std::vector<T>& v ) {
                    return v * w;
                }

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
                        point *= _vals.at( n );
                        return *this;
                    }
                }

                virtual const _window& apply( std::vector<T>& signal ) const {
                    bool size_matches = ( signal.size() == this->size() );
                    for (size_t n = 0; n < signal.size(); ++n) {
                        if (size_matches) {
                            signal[ n ] *= _vals.at( n );
                        } else {
                            this->apply( signal.size(), n, signal[ n ] );
                        }
                    }
                    return *this;
                }

                virtual _window& build( size_t npts ) {
                    _vals.resize( npts, 1.0 );
                    for (size_t i = 0; i < npts; ++i) {
                        _vals.at( i ) *= this->value( npts, i );
                    }
                    return *this;
                }

                virtual bool empty() const { return _vals.empty(); }

                virtual size_t size() const { return _vals.size(); }

                T& operator[]( size_t n ) { return _vals[ n ]; }

                const T& operator[]( size_t n ) const { return _vals[ n ]; }

            protected:
                void _init( size_t npts, const T& val ) {
                    _vals.resize( npts, static_cast<T>( val ) );
                }

            private:
                std::vector<T> _vals;
        };

        template<typename T>
        std::vector<T>& operator*=( std::vector<T>& v, const _window<T>& w ) {
            if (v.size() != w.size()) {
                throw std::range_error(
                    "Incompatible vector sizes for window multiplication" );
            }
            for (size_t i = 0; i < v.size(); ++i) {
                v[ i ] *= w[ i ];
            }
            return v;
        }

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
                    swap( *this, other );
                }

                virtual ~CosineWindow() {}

                friend void swap( CosineWindow<T>& a,
                                  CosineWindow<T>& b ) noexcept {
                    using std::swap;
                    swap( static_cast<_window<T>&>( a ),
                          static_cast<_window<T>&>( b ) );
                    swap( a._a, b._a );
                }

                virtual T value( size_t winsize, size_t n ) const override {
                    T w       = 0;
                    T kfactor = 1.0;
                    for (size_t k = 0; k < _a.size(); ++k) {
                        w += kfactor * _a[ k ]
                           * std::cos( 2.0 * NCPA::constants::PI * k
                                       * static_cast<T>( n )
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
                    swap( *this, other );
                }

                virtual ~HannWindow() {}

                friend void swap( HannWindow<T>& a,
                                  HannWindow<T>& b ) noexcept {
                    swap( static_cast<CosineWindow<T>&>( a ),
                          static_cast<CosineWindow<T>&>( b ) );
                }
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
                    swap( *this, other );
                }

                virtual ~HammingWindow() {}

                friend void swap( HammingWindow<T>& a,
                                  HammingWindow<T>& b ) noexcept {
                    swap( static_cast<CosineWindow<T>&>( a ),
                          static_cast<CosineWindow<T>&>( b ) );
                }
        };
    }  // namespace dsp
}  // namespace NCPA
