#pragma once

#include "NCPA/processing/AbstractDataWrapper.hpp"
#include "NCPA/processing/DataWrapper.hpp"
#include "NCPA/processing/declarations.hpp"

#include <memory>
#include <string>

void swap( NCPA::processing::Packet& a, NCPA::processing::Packet& b ) noexcept;

namespace NCPA::processing {
    class Packet {
        public:
            Packet() {}

            Packet( const std::string& tag ) : _tag{ tag } {}

            Packet( const Packet& other ) : _tag { other._tag } {}

            Packet( Packet&& other ) noexcept {
                ::swap( *this, other );
            }

            virtual ~Packet() {}

            friend void ::swap( Packet& a, Packet& b ) noexcept;


            virtual Packet& set_tag( const std::string& tag ) {
                _tag = tag;
                return *this;
            }

            virtual std::string tag() const { return _tag; }

        private:
            std::string _tag;
    };
}  // namespace NCPA::processing

void swap( NCPA::processing::Packet& a,
           NCPA::processing::Packet& b ) noexcept {
    using std::swap;
    swap( a._tag, b._tag );
}
