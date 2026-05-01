#pragma once

/**
     * Examples using armadillo classes

    class Double2ComplexDoubleConverter
        : public ProcessingStep<arma::Col<double>,
                                arma::Col<std::complex<double>>> {
        public:
            Double2ComplexDoubleConverter() {}

            virtual ~Double2ComplexDoubleConverter() {}

        protected:
            virtual bool _process_internal() {
                _product.set(
                    arma::conv_to<arma::Mat<std::complex<double>>>::from(
                        this->_input.get() ) );
                return true;
            }
    };

    class ComplexDoubleInverter
        : public ProcessingStep<arma::Col<std::complex<double>>,
                                arma::Col<std::complex<double>>> {
        public:
            ComplexDoubleInverter() {}

            virtual ~ComplexDoubleInverter() {}

        protected:
            virtual bool _process_internal() {
                _product = _input;
                _product.get().transform( []( std::complex<double> val ) {
                    return std::complex( val.imag(), val.real() );
                } );
                return true;
            }
    };
    */

#include "NCPA/processing/AbstractDataWrapper.hpp"
#include "NCPA/processing/AbstractProcessingStep.hpp"
#include "NCPA/processing/declarations.hpp"
#include "NCPA/processing/packets.hpp"

#include <unordered_map>

template<typename intype, typename outtype>
void swap( NCPA::processing::ProcessingStep<intype, outtype>& a,
           NCPA::processing::ProcessingStep<intype, outtype>& b ) noexcept;

namespace NCPA {
    namespace processing {
        template<typename intype, typename outtype>
        class ProcessingStep : public AbstractProcessingStep {
            public:
                ProcessingStep( const std::string& tag,
                                bool shortcircuit = false ) :
                    AbstractProcessingStep( tag ),
                    _short_circuit { shortcircuit } {}

                ProcessingStep() :
                    AbstractProcessingStep(), _short_circuit { false } {}

                virtual ~ProcessingStep() {}

                ProcessingStep(
                    const ProcessingStep<intype, outtype>& other ) :
                    ProcessingStep<intype, outtype>() {
                    _short_circuit   = other._short_circuit;
                    _input           = other._input;
                    _product         = other._product;
                    _parameters      = other._parameters;
                    _input_data_time = other._input_data_time;
                }

                ProcessingStep(
                    ProcessingStep<intype, outtype>&& other ) noexcept :
                    ProcessingStep<intype, outtype>() {
                    ::swap( *this, other );
                }

                friend void ::swap<>(
                    ProcessingStep<intype, outtype>& a,
                    ProcessingStep<intype, outtype>& b ) noexcept;

                ProcessingStep& operator=(
                    ProcessingStep<intype, outtype> other ) {
                    ::swap( *this, other );
                    return *this;
                }

                virtual AbstractDataWrapper& product() override {
                    return _product;
                }

                // virtual response_ptr_t pass_to_next(
                //     InputPacket& packet,
                //     response_id_t response_type_if_no_next,
                //     const std::string& error_msg_if_no_next ) override {
                //     return this->has_next()
                //              ? this->next()->process( packet )
                //              : this->_build_product_packet();
                // }

                // virtual response_ptr_t pass_on( InputPacket& packet,
                //                                 const std::string& error_msg
                //                                 = "" ) override {
                //     return this->has_next()
                //              ? this->next()->process( packet )
                //              : this->_build_product_packet();
                // }

                virtual bool product_available() const override {
                    return ( _product );
                }


            protected:
                virtual response_ptr_t _build_product_packet() const override {
                    return response_ptr_t(
                        new ProductPacket<outtype>( _product.ptr() ) );
                }

                virtual input_ptr_t _build_next_input_packet() const override {
                    return input_ptr_t(
                        new DataPacket<outtype>( _product.ptr() ) );
                }

                virtual packet_processing_result_t _process_data_packet(
                    InputPacket& packet, std::string& message ) override {
                    if (auto packet_ptr
                        = dynamic_cast<DataPacket<intype> *>( &packet )) {
                        _input.set( packet_ptr->ptr() );
                        _input_data_time = packet_ptr->interval();
                        if (this->_configuration_changed) {
                            if (!this->apply_configuration()) {}
                            this->_configuration_changed = false;
                        }
                        return this->_process_input();
                    } else {
                        return packet_processing_result_t::PACKET_INVALID;
                    }
                }

                bool _short_circuit;
                DataWrapper<intype> _input;
                DataWrapper<outtype> _product;
                parameter_tree_t _parameters;
                time_interval_t _input_data_time;
        };
    }  // namespace processing
}  // namespace NCPA

template<typename intype, typename outtype>
void swap( NCPA::processing::ProcessingStep<intype, outtype>& a,
           NCPA::processing::ProcessingStep<intype, outtype>& b ) noexcept {
    using std::swap;
    ::swap( static_cast<NCPA::processing::AbstractProcessingStep&>( a ),
            static_cast<NCPA::processing::AbstractProcessingStep&>( b ) );
    swap( a._short_circuit, b._short_circuit );
    swap( a._input, b._input );
    swap( a._product, b._product );
    swap( a._parameters, b._parameters );
    swap( a._input_data_time, b._input_data_time );
}
