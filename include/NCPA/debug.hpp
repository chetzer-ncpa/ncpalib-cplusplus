#pragma once

#include "NCPA/defines.hpp"
#include "NCPA/linearalgebra.hpp"
#include "NCPA/ndvector.hpp"

#include <complex>
#include <fstream>
#include <iostream>
#include <sstream>

namespace NCPA {
    namespace debug {

        template<typename T, typename U>
        void print_debug_matrix( const std::string& base, size_t step,
                                 const std::vector<T>& x1,
                                 const std::vector<T>& x2,
                                 const NCPA::arrays::ndvector<2, U>& var,
                                 ENABLE_FUNCTION_IF_ARITHMETIC( U ) ) {
            std::ostringstream filenamestr;
            filenamestr << base << "_" << step << ".debug";
            std::ofstream wcstr( filenamestr.str(), std::ios::trunc );
            for (size_t ii = 0; ii < x1.size(); ++ii) {
                for (size_t jj = 0; jj < x2.size(); ++jj) {
                    wcstr << x1[ ii ] << " " << x2[ jj ] << " "
                          << var[ ii ][ jj ] << std::endl;
                }
            }
            wcstr.close();
        }

        template<typename U>
        void print_debug_matrix( const std::string& base, size_t step,
                                 //  const std::vector<T>& x1,
                                 //  const std::vector<T>& x2,
                                 NCPA::linear::Matrix<U>& var,
                                 ENABLE_FUNCTION_IF_ARITHMETIC( U ) ) {
            std::ostringstream filenamestr;
            filenamestr << base << "_" << step << ".debug";
            std::ofstream wcstr( filenamestr.str(), std::ios::trunc );
            var.print_nonzero( wcstr );
            wcstr.close();
        }

        template<typename T, typename U>
        void print_debug_matrix( const std::string& base, size_t step,
                                 const std::vector<T>& x1,
                                 const std::vector<T>& x2,
                                 const NCPA::arrays::ndvector<2, U>& var,
                                 ENABLE_FUNCTION_IF_COMPLEX( U ) ) {
            std::ostringstream filenamestr;
            filenamestr << base << "_" << step << ".debug";
            std::ofstream wcstr( filenamestr.str(), std::ios::trunc );
            for (size_t ii = 0; ii < x1.size(); ++ii) {
                for (size_t jj = 0; jj < x2.size(); ++jj) {
                    wcstr << x1[ ii ] << " " << x2[ jj ] << " "
                          << var[ ii ][ jj ].real() << " "
                          << var[ ii ][ jj ].imag() << std::endl;
                }
            }
            wcstr.close();
        }

        template<typename U>
        void print_debug_matrix( const std::string& base, size_t step,
                                 NCPA::linear::Matrix<U>& var,
                                 ENABLE_FUNCTION_IF_COMPLEX( U ) ) {
            std::ostringstream filenamestr;
            filenamestr << base << "_" << step << ".debug";
            std::ofstream wcstr( filenamestr.str(), std::ios::trunc );
            var.print_nonzero( wcstr );
            wcstr.close();
        }

        template<typename U>
        void print_debug_vector( const std::string& base, size_t step,
                                 NCPA::linear::Vector<U>& var,
                                 ENABLE_FUNCTION_IF_ARITHMETIC( U ) ) {
            std::ostringstream filenamestr;
            filenamestr << base << "_" << step << ".debug";
            std::ofstream wcstr( filenamestr.str(), std::ios::trunc );
            var.print_nonzero( wcstr );
            wcstr.close();
        }

        template<typename U>
        void print_debug_vector( const std::string& base, size_t step,
                                 NCPA::linear::Vector<U>& var,
                                 ENABLE_FUNCTION_IF_COMPLEX( U ) ) {
            std::ostringstream filenamestr;
            filenamestr << base << "_" << step << ".debug";
            std::ofstream wcstr( filenamestr.str(), std::ios::trunc );
            var.print_nonzero( wcstr );
            wcstr.close();
        }
    }  // namespace debug
}  // namespace NCPA
