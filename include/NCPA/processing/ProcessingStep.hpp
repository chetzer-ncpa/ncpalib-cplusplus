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

namespace NCPA {
    namespace processing {
        template<typename intype, typename outtype>
        class ProcessingStep : public virtual AbstractProcessingStep {
            public:
                using input_t   = intype;
                using output_t  = outtype;
                using product_t = outtype;

                ProcessingStep( const std::string& tag,
                                bool shortcircuit = false ) :
                    AbstractProcessingStep( tag ),
                    _short_circuit { shortcircuit } {}

                ProcessingStep() : _short_circuit { false } {}

                virtual ~ProcessingStep() {}

                ProcessingStep(
                    const ProcessingStep<intype, outtype>& other ) :
                    AbstractProcessingStep( other ) {
                    _short_circuit   = other._short_circuit;
                    _input           = other._input;
                    _product         = other._product;
                    _input_data_time = other._input_data_time;

                    _parameters.clear();
                    for (auto it = other._parameters.begin();
                         it != other._parameters.end(); ++it) {
                        const std::vector<parameter_ptr_t>& vec = it->second;
                        std::vector<parameter_ptr_t> new_vec;
                        new_vec.reserve( vec.size() );
                        for (auto vit = vec.begin(); vit != vec.end(); ++vit) {
                            if (*vit) {
                                new_vec.push_back( ( *vit )->clone() );
                            }
                        }
                        _parameters.insert( std::make_pair(
                            it->first, std::move( new_vec ) ) );
                    }
                }

                ProcessingStep(
                    ProcessingStep<intype, outtype>&& other ) noexcept {
                    swap( *this, other );
                }

                friend void swap(
                    ProcessingStep<intype, outtype>& a,
                    ProcessingStep<intype, outtype>& b ) noexcept {
                    using std::swap;
                    swap( static_cast<AbstractProcessingStep&>( a ),
                          static_cast<AbstractProcessingStep&>( b ) );
                    swap( a._short_circuit, b._short_circuit );
                    swap( a._input, b._input );
                    swap( a._product, b._product );
                    swap( a._parameters, b._parameters );
                    swap( a._input_data_time, b._input_data_time );
                }

                virtual DataWrapper<intype>& input() override {
                    return _input;
                }

                virtual intype& input_data() { return _input.contents(); }

                virtual const intype& input_data() const {
                    return _input.contents();
                }

                virtual DataWrapper<outtype>& product() override {
                    return _product;
                }

                virtual outtype& product_data() { return _product.contents(); }

                virtual const outtype& product_data() const {
                    return _product.contents();
                }

                virtual bool product_available() const override {
                    return (bool)_product;
                }

                virtual ProcessingStep<intype, outtype>& reset() override {
                    return *this;
                }

            protected:
                virtual response_ptr_t _build_product_packet() const override {
                    return response_ptr_t(
                        new ProductPacket<outtype>( _product.ptr() ) );
                }

                virtual response_ptr_t _build_product_packet_from_input()
                    const override {
                    return response_ptr_t(
                        new ProductPacket<intype>( _input.ptr() ) );
                }

                virtual input_ptr_t _build_next_input_packet() const override {
                    return input_ptr_t(
                        new DataPacket<outtype>( _product.ptr() ) );
                }

                virtual input_ptr_t _ditto_input_packet() const override {
                    return input_ptr_t(
                        new DataPacket<intype>( _input.ptr() ) );
                }

                virtual packet_processing_result_t _process_data_packet(
                    InputPacket& packet,
                    std::vector<std::string>& message ) override {
                    if (auto packet_ptr
                        = dynamic_cast<DataPacket<intype> *>( &packet )) {
                        _input.set( packet_ptr->ptr() );
                        _input_data_time = packet_ptr->interval();
                        if (this->_configuration_changed) {
                            if (!this->apply_configuration( message )) {}
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
