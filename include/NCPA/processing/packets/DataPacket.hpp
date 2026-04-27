#pragma once

#include "NCPA/processing/packets/InputPacket.hpp"

#include <chrono>
#include <memory>

template<typename T>
void swap( NCPA::processing::DataPacket<T>& a,
           NCPA::processing::DataPacket<T>& b ) noexcept;

namespace NCPA::processing {


    // typedef std::chrono::duration<double, std::ratio<1>> duration_t;

    template<typename T>
    class DataPacket : public InputPacket,
                       public std::enable_shared_from_this<DataPacket<T>> {
        public:
            DataPacket() : InputPacket( input_id_t::DATA ) {}

            DataPacket( const T& in ) : InputPacket( input_id_t::DATA ) {
                _internal.set( in );
            }

            DataPacket( std::shared_ptr<T> in ) :
                InputPacket( input_id_t::DATA ) {
                _internal.set( in );
            }

            DataPacket( DataWrapper<T>& in ) :
                InputPacket( input_id_t::DATA ) {
                _internal.set( in.ptr() );
            }

            DataPacket( const T& in, const std::string& tag ) :
                InputPacket( input_id_t::DATA, tag ) {
                _internal.set( in );
            }

            DataPacket( std::shared_ptr<T> in, const std::string& tag ) :
                InputPacket( input_id_t::DATA, tag ) {
                // _internal = in;
                _internal.set( in );
            }

            DataPacket( DataWrapper<T>& in, const std::string& tag ) :
                InputPacket( input_id_t::DATA, tag ) {
                // _internal = in;
                _internal.set( in.ptr() );
            }

            DataPacket( const T& in, const time_point_t& start_time,
                        std::chrono::nanoseconds duration
                        = duration_t::zero() ) :
                InputPacket( input_id_t::DATA ),
                _data_time { start_time,
                             std::chrono::duration_cast<std::chrono::seconds>(
                                 duration ) } {
                _internal.set( in );
            }

            DataPacket( std::shared_ptr<T> in, const time_point_t& start_time,
                        std::chrono::nanoseconds duration
                        = duration_t::zero() ) :
                InputPacket( input_id_t::DATA ),
                _data_time { start_time,
                             std::chrono::duration_cast<std::chrono::seconds>(
                                 duration ) } {
                _internal.set( in );
            }

            DataPacket( DataWrapper<T>& in, const time_point_t& start_time,
                        std::chrono::nanoseconds duration
                        = duration_t::zero() ) :
                InputPacket( input_id_t::DATA ),
                _data_time { start_time,
                             std::chrono::duration_cast<std::chrono::seconds>(
                                 duration ) } {
                _internal.set( in.ptr() );
            }

            DataPacket( const T& in, const std::string& tag,
                        const time_point_t& start_time,
                        std::chrono::nanoseconds duration
                        = duration_t::zero() ) :
                InputPacket( input_id_t::DATA, tag ),
                _data_time { start_time,
                             std::chrono::duration_cast<std::chrono::seconds>(
                                 duration ) } {
                _internal.set( in );
            }

            DataPacket( std::shared_ptr<T> in, const std::string& tag,
                        const time_point_t& start_time,
                        std::chrono::nanoseconds duration
                        = duration_t::zero() ) :
                InputPacket( in, tag ),
                _data_time { start_time,
                             std::chrono::duration_cast<std::chrono::seconds>(
                                 duration ) } {
                _internal.set( in );
            }

            DataPacket( DataWrapper<T>& in, const std::string& tag,
                        const time_point_t& start_time,
                        std::chrono::nanoseconds duration
                        = duration_t::zero() ) :
                InputPacket( in, tag ),
                _data_time { start_time,
                             std::chrono::duration_cast<std::chrono::seconds>(
                                 duration ) } {
                _internal.set( in.ptr() );
            }

            // DataPacket( const DataPacket<T>& input ) : DataPacket<T>() {
            //     _internal = input.ptr();
            // }

            DataPacket( const DataPacket<T>& other ) : InputPacket( other ) {
                _internal  = other._internal;
                _data_time = other._data_time;
                // _duration  = other._duration;
            }

            DataPacket( DataPacket<T>&& input ) noexcept : DataPacket<T>() {
                ::swap( *this, input );
            }

            virtual ~DataPacket() {}

            friend void swap<>( DataPacket<T>& a, DataPacket<T>& b ) noexcept;

            DataPacket<T>& operator=( DataPacket<T> other ) {
                ::swap( *this, other );
                return *this;
            }

            const duration_t& duration() const { return _data_time.duration; }

            const time_interval_t& interval() const { return _data_time; }

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
                    return _internal.get();
                } else {
                    throw std::logic_error( "DataPacket: Nothing has been "
                                            "assigned to internal pointer" );
                }
            }

            const T& get() const {
                if (_internal) {
                    return _internal.get();
                } else {
                    throw std::logic_error( "DataPacket: Nothing has been "
                                            "assigned to internal pointer" );
                }
            }

            std::shared_ptr<T> ptr() { return _internal.ptr(); }

            const time_point_t& time() const { return _data_time.time; }

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
            time_interval_t _data_time;
            // duration_t _duration = duration_t::zero();
            // time_point_t _data_time;
    };
}  // namespace NCPA::processing

template<typename T>
void ::swap( NCPA::processing::DataPacket<T>& a,
             NCPA::processing::DataPacket<T>& b ) noexcept {
    using std::swap;
    ::swap( dynamic_cast<NCPA::processing::InputPacket&>( a ),
            dynamic_cast<NCPA::processing::InputPacket&>( b ) );
    swap( a._internal, b._internal );
    swap( a._data_time, b._data_time );
    // swap( a._duration, b._duration );
}
