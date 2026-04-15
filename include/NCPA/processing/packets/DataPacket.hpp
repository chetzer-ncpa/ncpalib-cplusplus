#pragma once

#include "NCPA/processing/packets/InputPacket.hpp"

#include <memory>

template<typename T>
void swap( NCPA::processing::DataPacket<T>& a,
           NCPA::processing::DataPacket<T>& b ) noexcept;

namespace NCPA::processing {

    template<typename T>
    class DataPacket : public InputPacket,
                       public std::enable_shared_from_this<DataPacket<T>> {
        public:
            DataPacket() : InputPacket( input_id_t::DATA ) {
                // _internal = std::make_shared<T>();
            }

            DataPacket( const T& in ) : InputPacket( input_id_t::DATA ) {
                // _internal = std::make_shared<T>( in );
                _internal.set( in );
            }

            DataPacket( const T& in, const std::string& tag ) :
                InputPacket( input_id_t::DATA, tag ) {
                // _internal = std::make_shared<T>( in );
                _internal.set( in );
            }

            DataPacket( std::shared_ptr<T> in ) :
                InputPacket( input_id_t::DATA ) {
                // _internal = in;
                _internal.set( in );
            }

            DataPacket( std::shared_ptr<T> in, const std::string& tag ) :
                InputPacket( input_id_t::DATA, tag ) {
                // _internal = in;
                _internal.set( in );
            }

            // DataPacket( const DataPacket<T>& input ) : DataPacket<T>() {
            //     _internal = input.ptr();
            // }

            DataPacket( DataPacket<T>&& input ) noexcept : DataPacket<T>() {
                ::swap( *this, input );
            }

            virtual ~DataPacket() {}

            friend void swap<>( DataPacket<T>& a, DataPacket<T>& b ) noexcept;

            // DataPacket<T>& operator=( DataPacket<T> other ) {
            //     ::swap( *this, other );
            //     return *this;
            // }

            DataPacket<T>& set( const T& input ) {
                // _internal = std::make_unique<T>( input );
                _internal.set( input );
                return *this;
            }

            DataPacket<T>& set( std::shared_ptr<T> input ) {
                // _internal = std::make_unique<T>( input );
                _internal = input;
                return *this;
            }

            T& get() {
                if (_internal) {
                    return *_internal;
                } else {
                    throw std::logic_error( "DataPacket: Nothing has been "
                                            "assigned to internal pointer" );
                }
            }

            const T& get() const {
                if (_internal) {
                    return *_internal;
                } else {
                    throw std::logic_error( "DataPacket: Nothing has been "
                                            "assigned to internal pointer" );
                }
            }

            std::shared_ptr<T> ptr() {
                return _internal.ptr();
                // return std::shared_ptr<T>( this->shared_from_this(),
                //                            _internal.get() );
            }

            static std::unique_ptr<InputPacket> build() {
                return std::unique_ptr<InputPacket>( new DataPacket<T>() );
            }

            static std::unique_ptr<InputPacket> build( const T in ) {
                return std::unique_ptr<InputPacket>( new DataPacket<T>( in ) );
            }

            static std::unique_ptr<InputPacket> build(
                const T in, const std::string& tag ) {
                return std::unique_ptr<InputPacket>(
                    new DataPacket<T>( in, tag ) );
            }

            static std::unique_ptr<InputPacket> build(
                const DataPacket<T>& input ) {
                return std::unique_ptr<InputPacket>(
                    new DataPacket<T>( input ) );
            }

        private:
            // std::shared_ptr<T> _internal;
            DataWrapper<T> _internal;
    };
}  // namespace NCPA::processing

template<typename T>
void ::swap( NCPA::processing::DataPacket<T>& a,
             NCPA::processing::DataPacket<T>& b ) noexcept {
    using std::swap;
    ::swap( dynamic_cast<NCPA::processing::InputPacket&>( a ),
            dynamic_cast<NCPA::processing::InputPacket&>( b ) );
    swap( a._internal, b._internal );
}
