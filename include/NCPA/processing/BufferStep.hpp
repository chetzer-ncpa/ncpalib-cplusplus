#pragma once

#include "NCPA/DAMES/declarations.hpp"

#include <deque>

namespace NCPA::DAMES {
    using namespace NCPA::processing;

    template<typename contents_t>
    class BufferStep
        : public StatefulStep<contents_t, std::deque<contents_t>> {
        public:
            BufferStep() : BufferStep( "" ) {}

            BufferStep( const std::string& tag ) :
                StatefulStep<contents_t, std::deque<contents_t>>( tag ) {
                _define_parameters();
            }

            virtual ~BufferStep()                                 = default;
            BufferStep( const BufferStep<contents_t>& other )     = default;
            BufferStep( BufferStep<contents_t>&& other ) noexcept = default;

            bool apply_configuration(
                std::vector<std::string>& message ) override {
                this->set_capacity( this->parameter( "buffer_size" ).asInt() );
                return true;
            }

            int capacity() const { return _capacity; }

            BufferStep<contents_t>& clear() {
                this->state().clear();
                return *this;
            }

            std::deque<contents_t>& contents() { return this->state(); }

            const std::deque<contents_t>& contents() const {
                return this->state();
            }

            bool run_tests( std::vector<std::string>& message ) override {
                bool passed = true;
                int n       = 5;
                this->clear().set_capacity( n );
                this->input_data() = contents_t();
                for (int i = 0; i < n; ++i) {
                    this->_process_input();
                    if (this->state().size() != i + 1) {
                        message.push_back(
                            "Fed " + std::to_string( i + 1 )
                            + " frames but buffer reports size "
                            + std::to_string( this->state().size() ) + "." );
                        passed = false;
                    }
                }

                for (int i = 0; i < n; ++i) {
                    this->_process_input();
                    if (this->state().size() != n) {
                        message.push_back(
                            "Buffer set to max size " + std::to_string( n )
                            + " frames but buffer reports size "
                            + std::to_string( this->state().size() ) + "." );
                        passed = false;
                    }
                }
                return passed;
            }

            BufferStep<contents_t>& set_capacity( int i ) {
                _capacity = i;
                if (this->state().size() > i) {
                    this->state().resize( i );
                }
                return *this;
            }

            int size() const { return this->state().size(); }

        protected:
            int _capacity = 1;

            virtual void _define_parameters() override {
                BooleanParameter verbose( "verbose", false,
                                          "Use verbose output" );
                IntegerParameter buffer_size( "buffer_size", 1,
                                              "Number of items to keep" );

                this->add_parameter( verbose ).add_parameter( buffer_size );
            }

            virtual packet_processing_result_t _process_input() override {
                this->_flags.clear();
                this->state().push_back( this->input_data() );
                if (this->size() > this->capacity()) {
                    this->state().pop_front();
                }
                this->_product.set( this->input_data() );
                return packet_processing_result_t::SUCCESS_PRODUCT;
            }
    };


}  // namespace NCPA::DAMES
