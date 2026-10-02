#pragma once

#include "NCPA/processing/DataWrapper.hpp"
#include "NCPA/processing/declarations.hpp"
#include "NCPA/processing/packets/Packet.hpp"
#include "NCPA/strings.hpp"

#include <memory>
#include <sstream>
#include <vector>

// void swap( NCPA::processing::ResponsePacket& a,
//            NCPA::processing::ResponsePacket& b ) noexcept;

namespace NCPA::processing {
    class ResponsePacket : public Packet {
        public:
            ResponsePacket( response_id_t id ) : Packet(), _ID { id } {}

            ResponsePacket( response_id_t id, const std::string& tag ) :
                Packet( tag ), _ID { id }, _messages { "" } {}

            ResponsePacket( response_id_t id, const std::string& tag,
                            const std::vector<std::string>& message ) :
                Packet( tag ), _ID { id }, _messages { message } {}

            ResponsePacket( const ResponsePacket& other ) :
                Packet( other ),
                _ID { other._ID },
                _messages { other._messages },
                _flags { other._flags } {}

            virtual ~ResponsePacket() {}

            friend void swap( ResponsePacket& a, ResponsePacket& b ) noexcept {
                using std::swap;
                swap( dynamic_cast<NCPA::processing::Packet&>( a ),
                      dynamic_cast<NCPA::processing::Packet&>( b ) );
                swap( a._ID, b._ID );
                swap( a._messages, b._messages );
                swap( a._flags, b._flags );
            }

            response_id_t& ID() { return _ID; }

            const response_id_t& ID() const { return _ID; }

            std::string message() const {
                std::ostringstream oss;
                oss << NCPA::strings::join( _messages, "\n" );
                if (!_flags.empty()) {
                    oss << "\nFlags:\n" << NCPA::strings::join( _flags, "\n" );
                }
                return oss.str();
            }

            ResponsePacket& append_flag( const std::string& message ) {
                _flags.push_back( message );
                return *this;
            }

            ResponsePacket& prepend_flag( const std::string& message ) {
                _flags.insert( _flags.begin(), message );
                return *this;
            }

            const std::vector<std::string>& flags() const { return _flags; }

            static response_ptr_t build( response_id_t ID ) {
                return response_ptr_t( new ResponsePacket( ID ) );
            }

            static response_ptr_t build( response_id_t ID,
                                         const std::string& tag ) {
                return response_ptr_t( new ResponsePacket( ID, tag ) );
            }

            static response_ptr_t build( response_id_t ID,
                                         const std::string& tag,
                                         const std::vector<std::string>& message ) {
                return response_ptr_t(
                    new ResponsePacket( ID, tag, message ) );
            }

        private:
            response_id_t _ID;
            std::vector<std::string> _messages;
            std::vector<std::string> _flags;
    };
}  // namespace NCPA::processing
