#pragma once


#include "NCPA/atmosphere/Atmosphere1D.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/cloneable.hpp"
#include "NCPA/files.hpp"

#include <fstream>
#include <iostream>
#include <memory>
#include <vector>

namespace NCPA {
    namespace atmos {

        /**
         * @brief
         */
        class _abstract_atmosphere_reader_1d
            : public Cloneable<_abstract_atmosphere_reader_1d> {
            public:
                /**
                 * @brief
                 */
                virtual ~_abstract_atmosphere_reader_1d() {}

                friend void swap(
                    _abstract_atmosphere_reader_1d& a,
                    _abstract_atmosphere_reader_1d& b ) noexcept {}

                /**
                 * @brief
                 * @param in1
                 * @return Atmosphere1D
                 */
                virtual Atmosphere1D read( std::istream& in1 ) = 0;

                /**
                 * @brief
                 * @param streams
                 * @return Atmosphere1D
                 */
                virtual Atmosphere1D read(
                    std::vector<std::istream *>& streams ) = 0;

                /**
                 * @brief
                 * @param filename
                 * @return Atmosphere1D
                 */
                virtual Atmosphere1D read( const std::string& filename ) = 0;

                /**
                 * @brief
                 * @param filenames
                 * @return Atmosphere1D
                 */
                virtual Atmosphere1D read(
                    const std::vector<std::string>& filenames ) = 0;
        };
    }  // namespace atmos
}  // namespace NCPA
