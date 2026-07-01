#pragma once

#include "NCPA/processing/declarations.hpp"
#include "NCPA/processing/parameters/declarations.hpp"
#include "NCPA/processing/parameters/Parameter.hpp"
#include "NCPA/processing/parameters/VectorParameter.hpp"

#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#if HAVE_NLOHMANN_JSON_HPP
#  include "nlohmann/json.hpp"
#endif

namespace NCPA {
    namespace processing {
        class BooleanVectorParameter : public VectorParameter<bool> {
            public:
                BooleanVectorParameter( const std::string& key,
                                        const std::vector<bool>& value,
                                        const std::string& description = "" ) :
                    VectorParameter<bool>( key, value, description ) {}

                BooleanVectorParameter( const std::string& key,
                                        std::initializer_list<bool>& value,
                                        const std::string& description = "" ) :
                    VectorParameter<bool>( key, value, description ) {}

                virtual ~BooleanVectorParameter() {}

                friend void swap(
                    NCPA::processing::BooleanVectorParameter& a,
                    NCPA::processing::BooleanVectorParameter& b ) noexcept {
                    using std::swap;
                    swap(
                        dynamic_cast<NCPA::processing::VectorParameter<bool>&>(
                            a ),
                        dynamic_cast<NCPA::processing::VectorParameter<bool>&>(
                            b ) );
                }

                virtual std::unique_ptr<Parameter> clone() const override {
                    return std::unique_ptr<Parameter>(
                        new BooleanVectorParameter( this->key(), this->value(),
                                                    this->description() ) );
                }

                virtual parameter_type_t data_type() const override {
                    return parameter_type_t::BOOLEAN;
                }

                virtual std::vector<bool> asBoolVector() const override {
                    return this->value();
                }

                virtual std::vector<int> asIntVector() const override {
                    std::vector<int> vec;
                    for (size_t i = 0; i < this->value().size(); ++i) {
                        vec[ i ] = this->value().at( i ) ? 1 : 0;
                    }
                    return vec;
                }

                virtual std::vector<double> asDoubleVector() const override {
                    std::vector<double> vec;
                    for (size_t i = 0; i < this->value().size(); ++i) {
                        vec[ i ] = this->value().at( i ) ? 1.0 : 0.0;
                    }
                    return vec;
                }

                virtual std::vector<std::string> asStringVector()
                    const override {
                    std::vector<std::string> vec;
                    for (size_t i = 0; i < this->value().size(); ++i) {
                        vec[ i ] = this->value().at( i ) ? "true" : "false";
                    }
                    return vec;
                }

                virtual std::string data_type_str() const override {
                    return "boolean";
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key, const std::vector<bool>& value,
                    const std::string& description ) {
                    return std::unique_ptr<Parameter>(
                        new BooleanVectorParameter( key, value,
                                                    description ) );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key, const std::vector<bool>& value ) {
                    return std::unique_ptr<Parameter>(
                        new BooleanVectorParameter( key, value ) );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key,
                    const std::initializer_list<bool>& value,
                    const std::string& description ) {
                    return std::unique_ptr<Parameter>(
                        new BooleanVectorParameter( key, value,
                                                    description ) );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key,
                    const std::initializer_list<bool>& value ) {
                    return std::unique_ptr<Parameter>(
                        new BooleanVectorParameter( key, value ) );
                }
        };
    }  // namespace processing
}  // namespace NCPA
