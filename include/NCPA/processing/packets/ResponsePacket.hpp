#pragma once

#include "NCPA/processing/DataWrapper.hpp"
#include "NCPA/processing/declarations.hpp"
#include "NCPA/processing/packets/Packet.hpp"

#include <memory>
#include <vector>

void swap( NCPA::processing::ResponsePacket& a,
           NCPA::processing::ResponsePacket& b ) noexcept;

namespace NCPA::processing {
    class ResponsePacket : public Packet {
        public:
            ResponsePacket( response_id_t id ) : Packet(), _ID { id } {}

            ResponsePacket( response_id_t id, const std::string& tag ) :
                Packet( tag ), _ID { id }, _message { "" } {}

            ResponsePacket( response_id_t id, const std::string& tag,
                            const std::string& message ) :
                Packet( tag ), _ID { id }, _message { message } {}

            ResponsePacket( const ResponsePacket& other ) :
                Packet( other ),
                _ID { other._ID },
                _message { other._message },
                _flags { other._flags } {}

            virtual ~ResponsePacket() {}

            friend void ::swap( ResponsePacket& a,
                                ResponsePacket& b ) noexcept;

            response_id_t ID() const { return _ID; }

            std::string message() const {
                std::ostringstream oss;
                oss << _message << "\n" << "Flags:" << "\n";
                for (auto f : _flags) {
                    oss << f << std::endl;
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

            // static response_ptr_t build( response_id_t resptype,
            //                              const std::string& tag,
            //                              const std::string& message ) {
            //     return response_ptr_t(
            //         new ResponsePacket( resptype, tag, message ) );
            // }

        private:
            response_id_t _ID;
            std::string _message;
            std::vector<std::string> _flags;
    };
}  // namespace NCPA::processing

inline void swap( NCPA::processing::ResponsePacket& a,
                  NCPA::processing::ResponsePacket& b ) noexcept {
    using std::swap;
    ::swap( dynamic_cast<NCPA::processing::Packet&>( a ),
            dynamic_cast<NCPA::processing::Packet&>( b ) );
    swap( a._ID, b._ID );
    swap( a._message, b._message );
    swap( a._flags, b._flags );
}
