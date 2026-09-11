#pragma once

#include "NCPA/cloneable.hpp"
#include "NCPA/configuration/declarations.hpp"

#include <functional>
#include <string>
#include <utility>

namespace NCPA {
    namespace config {
        class Validation : public Cloneable<Validation> {
            public:
                Validation(const std::string& failmsg = "") : _failmessage{ failmsg } {}

                virtual ~Validation() {}

                Validation( const Validation& other ) {
                    _failmessage = other._failmessage;
                }

                Validation( Validation&& other ) noexcept {
                    swap( *this, other );
                }

                friend void swap( Validation& a, Validation& b ) noexcept {
                    using std::swap;
                    swap( a._failmessage, b._failmessage );
                }

                const std::string& failure_message() const {
                    return _failmessage;
                }

            private:
                std::string _failmessage;
        };

        

    }  // namespace config
}  // namespace NCPA
