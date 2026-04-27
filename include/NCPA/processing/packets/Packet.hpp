#pragma once

#include "NCPA/processing/AbstractDataWrapper.hpp"
#include "NCPA/processing/DataWrapper.hpp"
#include "NCPA/processing/declarations.hpp"

#include <algorithm>
#include <memory>
#include <string>

// void swap( NCPA::processing::Packet& a, NCPA::processing::Packet& b )
// noexcept;

namespace NCPA {
    namespace processing {
        class Packet {
            public:
                Packet() : Packet( "" ) {}

                Packet( const std::string& tag ) :
                    _tag { tag },
                    _packet_time { std::chrono::system_clock::now() } {}

                Packet( const Packet& other ) :
                    _tag { other._tag },
                    _packet_time { std::chrono::system_clock::now() } {}

                Packet( Packet&& other ) noexcept { swap( *this, other ); }

                virtual ~Packet() {}

                friend void swap( Packet& a, Packet& b ) noexcept;

                virtual time_point_t packet_time() const {
                    return _packet_time;
                }

                virtual Packet& set_tag( const std::string& tag ) {
                    _tag = tag;
                    return *this;
                }

                virtual std::string tag() const { return _tag; }

            private:
                std::string _tag;
                time_point_t _packet_time;
        };

        inline void swap( NCPA::processing::Packet& a,
                   NCPA::processing::Packet& b ) noexcept {
            using std::swap;
            swap( a._tag, b._tag );
            swap( a._packet_time, b._packet_time );
        }
    }  // namespace processing
}  // namespace NCPA
