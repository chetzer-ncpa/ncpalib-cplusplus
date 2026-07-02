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
        class IntegerVectorParameter : public VectorParameter<int> {
            public:
                IntegerVectorParameter( const std::string& key,
                                        const std::vector<int>& value,
                                        const std::string& description = "" ) :
                    VectorParameter<int>( key, value, description ) {}

                IntegerVectorParameter( const std::string& key,
                                        std::initializer_list<int>& value,
                                        const std::string& description = "" ) :
                    VectorParameter<int>( key, value, description ) {}

                virtual ~IntegerVectorParameter() {}

                friend void swap(
                    NCPA::processing::IntegerVectorParameter& a,
                    NCPA::processing::IntegerVectorParameter& b ) noexcept {
                    using std::swap;
                    swap(
                        dynamic_cast<NCPA::processing::VectorParameter<int>&>(
                            a ),
                        dynamic_cast<NCPA::processing::VectorParameter<int>&>(
                            b ) );
                }

                virtual std::unique_ptr<Parameter> clone() const override {
                    return std::unique_ptr<Parameter>(
                        new IntegerVectorParameter( this->key(), this->value(),
                                                    this->description() ) );
                }

                virtual parameter_type_t data_type() const override {
                    return parameter_type_t::INTEGER;
                }

                virtual std::vector<bool> asBoolVector() const override {
                    std::vector<bool> vec;
                    for (size_t i = 0; i < this->value().size(); ++i) {
                        vec[ i ] = ( this->value().at( i ) != 0 );
                    }
                    return vec;
                }

                virtual std::vector<int> asIntVector() const override {
                    return this->value();
                }

                virtual std::vector<double> asDoubleVector() const override {
                    std::vector<double> vec;
                    for (size_t i = 0; i < this->value().size(); ++i) {
                        vec[ i ]
                            = static_cast<double>( this->value().at( i ) );
                    }
                    return vec;
                }

                virtual std::vector<std::string> asStringVector()
                    const override {
                    std::vector<std::string> vec;
                    std::ostringstream oss;
                    for (size_t i = 0; i < this->value().size(); ++i) {
                        oss << this->value().at( i );
                        vec[ i ] = oss.str();
                        oss.str( "" );
                    }
                    return vec;
                }

                virtual std::string data_type_str() const override {
                    return "integer";
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key, const std::vector<int>& value,
                    const std::string& description ) {
                    return std::unique_ptr<Parameter>(
                        new IntegerVectorParameter( key, value,
                                                    description ) );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key, const std::vector<int>& value ) {
                    return std::unique_ptr<Parameter>(
                        new IntegerVectorParameter( key, value ) );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key,
                    const std::initializer_list<int>& value,
                    const std::string& description ) {
                    return std::unique_ptr<Parameter>(
                        new IntegerVectorParameter( key, value,
                                                    description ) );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key,
                    const std::initializer_list<int>& value ) {
                    return std::unique_ptr<Parameter>(
                        new IntegerVectorParameter( key, value ) );
                }
        };
    }  // namespace processing
}  // namespace NCPA
