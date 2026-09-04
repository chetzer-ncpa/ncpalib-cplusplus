#pragma once

/**
 * @file
 * @brief
 */

#include "NCPA/atmosphere/declarations.hpp"

#include <memory>
#include <stdexcept>

namespace NCPA {
    namespace atmos {

        /**
         * @brief
         */
        class abstract_atmospheric_property {
            public:
                /**
                 * @brief
                 */
                virtual ~abstract_atmospheric_property() {}

                friend void swap( abstract_atmospheric_property& a,
                                  abstract_atmospheric_property& b ) noexcept {
                }

                /**
                 * @brief
                 * @return size_t
                 */
                virtual size_t dimensions() const = 0;

                /**
                 * @brief
                 * @return std::unique_ptr<abstract_atmospheric_property>
                 */
                virtual std::unique_ptr<abstract_atmospheric_property> clone()
                    const = 0;

                /**
                 * @brief
                 * @param source
                 * @return abstract_atmospheric_property&
                 */
                virtual abstract_atmospheric_property& copy(
                    const abstract_atmospheric_property& source ) = 0;

                /**
                 * @brief
                 * @return true
                 * @return false
                 */
                virtual bool strict() const { return _strict; }

                /**
                 * @brief
                 * @param tf
                 * @return abstract_atmospheric_property&
                 */
                virtual abstract_atmospheric_property& set_strict( bool tf ) {
                    _strict = tf;
                    return *this;
                }

            protected:
                /**
                 * @brief
                 */
                bool _strict = false;

                /**
                 * @brief
                 * @param n
                 * @throws std::range_error
                 */
                virtual void _validate_axis( size_t n ) const {
                    if (n >= this->dimensions()) {
                        std::ostringstream oss;
                        oss << "Invalid axis " << n;
                        throw std::range_error( oss.str() );
                    }
                }
        };
    }  // namespace atmos
}  // namespace NCPA
