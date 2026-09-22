#include "NCPA/arrays.hpp"
#include "NCPA/constants.hpp"
#include "NCPA/dsp.hpp"
#include "NCPA/gtest.hpp"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <limits>
#include <numbers>
#include <utility>

using namespace std;
using namespace NCPA::dsp;
using namespace NCPA::arrays;
using namespace testing;

TEST( NCPADSPLibraryTest, HannWindowIsCorrect ) {
    size_t nwin = 11;
    HannWindow<double> hann( nwin );
    for (size_t i = 0; i < nwin; ++i) {
        EXPECT_DOUBLE_EQ( hann[ i ],
                          ( 0.5
                            - 0.5
                                  * cos( 2.0 * NCPA::constants::PI * (double)i
                                         / (double)( nwin - 1 ) ) ) );
    }
}

TEST( NCPADSPLibraryTest, HammingWindowIsCorrect ) {
    size_t nwin = 11;
    HammingWindow<double> hamm( nwin );
    for (size_t i = 0; i < nwin; ++i) {
        EXPECT_DOUBLE_EQ( hamm[ i ],
                          ( 0.54
                            - 0.46
                                  * cos( 2.0 * NCPA::constants::PI * (double)i
                                         / (double)( nwin - 1 ) ) ) );
    }
}

TEST( NCPADSPLibraryTest, WindowScalesVectorCorrectly ) {
    size_t nwin = 11;
    HammingWindow<double> hamm( nwin );
    std::vector<double> v( nwin, 1.0 );
    std::vector<double> w = v * hamm;
    for (size_t i = 0; i < nwin; ++i) {
        EXPECT_DOUBLE_EQ( hamm[ i ], w[ i ] );
    }
    v *= hamm;
    for (size_t i = 0; i < nwin; ++i) {
        EXPECT_DOUBLE_EQ( hamm[ i ], v[ i ] );
    }
}

#ifdef NCPA_HAVE_FFTW3
TEST( NCPADSPLibraryTest, FFTWGeneratesExpectedForwardTransform ) {
    size_t nfft    = 128;
    size_t ncycles = 8;
    double df
        = ( 2 * NCPA::constants::PI / ( (double)nfft / (double)ncycles ) );
    vector<complex<double>> input( nfft );
    for (size_t i = 0; i < nfft; ++i) {
        input[ i ].real( std::sin( df * i ) );
    }
    fft_ptr_t f1 = FFTFactory::build( fft_t::FFTW3 );
    vector<complex<double>> output;
    f1->compute( input, output, fft_sign_t::POSITIVE );

    // with these settings, bin 8 should be (0,NFFT/2),
    // bin 120 should be (0,-NFFT/2), and zero otherwise
    size_t posbin = nfft / 2 / ncycles;
    for (size_t i = 0; i < nfft; ++i) {
        if (i == posbin) {
            EXPECT_NEAR( output.at( i ).real(), 0.0, 1e-12 );
            EXPECT_DOUBLE_EQ( ( output.at( i ).imag() ),
                              (double)( nfft / 2 ) );
        } else if (i == ( nfft - posbin )) {
            EXPECT_NEAR( output.at( i ).real(), 0.0, 1e-12 );
            EXPECT_DOUBLE_EQ( ( output.at( i ).imag() ),
                              -(double)( nfft / 2 ) );
        } else {
            EXPECT_NEAR( output.at( i ).real(), 0.0, 1e-12 );
            EXPECT_NEAR( output.at( i ).imag(), 0.0, 1e-12 );
        }
    }

    f1->compute( input, output, fft_sign_t::NEGATIVE );
    for (size_t i = 0; i < nfft; ++i) {
        if (i == posbin) {
            EXPECT_NEAR( output.at( i ).real(), 0.0, 1e-12 );
            EXPECT_DOUBLE_EQ( ( output.at( i ).imag() ),
                              -(double)( nfft / 2 ) );
        } else if (i == ( nfft - posbin )) {
            EXPECT_NEAR( output.at( i ).real(), 0.0, 1e-12 );
            EXPECT_DOUBLE_EQ( ( output.at( i ).imag() ),
                              (double)( nfft / 2 ) );
        } else {
            EXPECT_NEAR( output.at( i ).real(), 0.0, 1e-12 );
            EXPECT_NEAR( output.at( i ).imag(), 0.0, 1e-12 );
        }
    }
}

TEST( NCPADSPLibraryTest, FFTWGeneratesExpectedInverseTransform ) {
    size_t nfft    = 128;
    size_t ncycles = 8;
    double df
        = ( 2 * NCPA::constants::PI / ( (double)nfft / (double)ncycles ) );
    vector<complex<double>> input( nfft );
    for (size_t i = 0; i < nfft; ++i) {
        input[ i ].real( std::sin( df * i ) );
    }
    fft_ptr_t f1 = FFTFactory::build( fft_t::FFTW3 );
    vector<complex<double>> output;
    f1->compute( input, output, fft_sign_t::POSITIVE );
    vector<complex<double>> rebuilt;
    f1->compute( output, rebuilt, fft_sign_t::NEGATIVE );

    for (size_t i = 0; i < input.size(); ++i) {
        EXPECT_NEAR( (double)nfft * input[ i ].real(), rebuilt[ i ].real(),
                     1e-12 );
        EXPECT_NEAR( (double)nfft * input[ i ].imag(), rebuilt[ i ].imag(),
                     1e-12 );
    }

    f1->compute( input, output, fft_sign_t::NEGATIVE );
    f1->compute( output, rebuilt, fft_sign_t::POSITIVE );

    for (size_t i = 0; i < input.size(); ++i) {
        EXPECT_NEAR( (double)nfft * input[ i ].real(), rebuilt[ i ].real(),
                     1e-12 );
        EXPECT_NEAR( (double)nfft * input[ i ].imag(), rebuilt[ i ].imag(),
                     1e-12 );
    }
}

TEST( NCPADSPLibraryTest, FFTWGeneratesExpectedForwardRealTransform ) {
    size_t nfft    = 128;
    size_t ncycles = 8;
    double df
        = ( 2 * NCPA::constants::PI / ( (double)nfft / (double)ncycles ) );
    vector<double> input( nfft );
    for (size_t i = 0; i < nfft; ++i) {
        input[ i ] = std::sin( df * i );
    }
    fft_ptr_t f1 = FFTFactory::build( fft_t::FFTW3 );
    vector<complex<double>> output;
    f1->compute( input, output, fft_sign_t::POSITIVE );

    // with these settings, bin posbin should be (0,NFFT/2),
    // bin NFFT-posbin should be (0,-NFFT/2), and zero otherwise.
    // sign is flipped from previous test because FFTW does real-to-complex
    // transforms as forward (sign-positive) only
    size_t posbin = nfft / 2 / ncycles;
    for (size_t i = 0; i < nfft; ++i) {
        if (i == posbin) {
            EXPECT_NEAR( output.at( i ).real(), 0.0, 1e-12 );
            EXPECT_DOUBLE_EQ( ( output.at( i ).imag() ),
                              -(double)( nfft / 2 ) );
        } else if (i == ( nfft - posbin )) {
            EXPECT_NEAR( output.at( i ).real(), 0.0, 1e-12 );
            EXPECT_DOUBLE_EQ( ( output.at( i ).imag() ),
                              (double)( nfft / 2 ) );
        } else {
            EXPECT_NEAR( output.at( i ).real(), 0.0, 1e-12 );
            EXPECT_NEAR( output.at( i ).imag(), 0.0, 1e-12 );
        }
    }
}
#endif

#ifdef NCPA_HAVE_POCKETFFT
TEST( NCPADSPLibraryTest, PocketFFTGeneratesExpectedForwardTransform ) {
    size_t nfft    = 128;
    size_t ncycles = 8;
    double df
        = ( 2 * NCPA::constants::PI / ( (double)nfft / (double)ncycles ) );
    vector<complex<double>> input( nfft );
    for (size_t i = 0; i < nfft; ++i) {
        input[ i ].real( std::sin( df * i ) );
    }
    fft_ptr_t f1 = FFTFactory::build( fft_t::POCKETFFT );
    vector<complex<double>> output;
    f1->compute( input, output, fft_sign_t::POSITIVE );

    // with these settings, bin 8 should be (0,NFFT/2),
    // bin 120 should be (0,-NFFT/2), and zero otherwise
    size_t posbin = nfft / 2 / ncycles;
    for (size_t i = 0; i < nfft; ++i) {
        if (i == posbin) {
            EXPECT_NEAR( output.at( i ).real(), 0.0, 1e-12 );
            EXPECT_DOUBLE_EQ( ( output.at( i ).imag() ),
                              (double)( nfft / 2 ) );
        } else if (i == ( nfft - posbin )) {
            EXPECT_NEAR( output.at( i ).real(), 0.0, 1e-12 );
            EXPECT_DOUBLE_EQ( ( output.at( i ).imag() ),
                              -(double)( nfft / 2 ) );
        } else {
            EXPECT_NEAR( output.at( i ).real(), 0.0, 1e-12 );
            EXPECT_NEAR( output.at( i ).imag(), 0.0, 1e-12 );
        }
    }

    f1->compute( input, output, fft_sign_t::NEGATIVE );
    for (size_t i = 0; i < nfft; ++i) {
        if (i == posbin) {
            EXPECT_NEAR( output.at( i ).real(), 0.0, 1e-12 );
            EXPECT_DOUBLE_EQ( ( output.at( i ).imag() ),
                              -(double)( nfft / 2 ) );
        } else if (i == ( nfft - posbin )) {
            EXPECT_NEAR( output.at( i ).real(), 0.0, 1e-12 );
            EXPECT_DOUBLE_EQ( ( output.at( i ).imag() ),
                              (double)( nfft / 2 ) );
        } else {
            EXPECT_NEAR( output.at( i ).real(), 0.0, 1e-12 );
            EXPECT_NEAR( output.at( i ).imag(), 0.0, 1e-12 );
        }
    }
}

TEST( NCPADSPLibraryTest, PocketFFTGeneratesExpectedInverseTransform ) {
    size_t nfft    = 128;
    size_t ncycles = 8;
    double df
        = ( 2 * NCPA::constants::PI / ( (double)nfft / (double)ncycles ) );
    vector<complex<double>> input( nfft );
    for (size_t i = 0; i < nfft; ++i) {
        input[ i ].real( std::sin( df * i ) );
    }
    fft_ptr_t f1 = FFTFactory::build( fft_t::POCKETFFT );
    vector<complex<double>> output;
    f1->compute( input, output, fft_sign_t::POSITIVE );
    vector<complex<double>> rebuilt;
    f1->compute( output, rebuilt, fft_sign_t::NEGATIVE );

    for (size_t i = 0; i < input.size(); ++i) {
        EXPECT_NEAR( (double)nfft * input[ i ].real(), rebuilt[ i ].real(),
                     1e-12 );
        EXPECT_NEAR( (double)nfft * input[ i ].imag(), rebuilt[ i ].imag(),
                     1e-12 );
    }

    f1->compute( input, output, fft_sign_t::NEGATIVE );
    f1->compute( output, rebuilt, fft_sign_t::POSITIVE );

    for (size_t i = 0; i < input.size(); ++i) {
        EXPECT_NEAR( (double)nfft * input[ i ].real(), rebuilt[ i ].real(),
                     1e-12 );
        EXPECT_NEAR( (double)nfft * input[ i ].imag(), rebuilt[ i ].imag(),
                     1e-12 );
    }
}

TEST( NCPADSPLibraryTest, PocketFFTGeneratesExpectedForwardRealTransform ) {
    size_t nfft    = 128;
    size_t ncycles = 8;
    double df
        = ( 2 * NCPA::constants::PI / ( (double)nfft / (double)ncycles ) );
    vector<double> input( nfft );
    for (size_t i = 0; i < nfft; ++i) {
        input[ i ] = std::sin( df * i );
    }
    fft_ptr_t f1 = FFTFactory::build( fft_t::POCKETFFT );
    vector<complex<double>> output;
    f1->compute( input, output, fft_sign_t::POSITIVE );

    // with these settings, bin posbin should be (0,NFFT/2),
    // bin NFFT-posbin should be (0,-NFFT/2), and zero otherwise.
    size_t posbin = nfft / 2 / ncycles;
    for (size_t i = 0; i < nfft; ++i) {
        if (i == posbin) {
            EXPECT_NEAR( output.at( i ).real(), 0.0, 1e-12 );
            EXPECT_DOUBLE_EQ( ( output.at( i ).imag() ),
                              (double)( nfft / 2 ) );
        } else if (i == ( nfft - posbin )) {
            EXPECT_NEAR( output.at( i ).real(), 0.0, 1e-12 );
            EXPECT_DOUBLE_EQ( ( output.at( i ).imag() ),
                              -(double)( nfft / 2 ) );
        } else {
            EXPECT_NEAR( output.at( i ).real(), 0.0, 1e-12 );
            EXPECT_NEAR( output.at( i ).imag(), 0.0, 1e-12 );
        }
    }
}
#endif
