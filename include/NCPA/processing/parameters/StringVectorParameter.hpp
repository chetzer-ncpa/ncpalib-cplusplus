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
        class StringVectorParameter : public VectorParameter<std::string> {
            public:
                StringVectorParameter( const std::string& key,
                                       const std::vector<std::string>& value,
                                       const std::string& description = "" ) :
                    VectorParameter<std::string>( key, value, description ) {}

                StringVectorParameter(
                    const std::string& key,
                    std::initializer_list<std::string>& value,
                    const std::string& description = "" ) :
                    VectorParameter<std::string>( key, value, description ) {}

                virtual ~StringVectorParameter() {}

                friend void swap(
                    NCPA::processing::StringVectorParameter& a,
                    NCPA::processing::StringVectorParameter& b ) noexcept {
                    using std::swap;
                    swap( dynamic_cast<
                              NCPA::processing::VectorParameter<std::string>&>(
                              a ),
                          dynamic_cast<
                              NCPA::processing::VectorParameter<std::string>&>(
                              b ) );
                }

                virtual std::unique_ptr<Parameter> clone() const override {
                    return std::unique_ptr<Parameter>(
                        new StringVectorParameter( this->key(), this->value(),
                                                   this->description() ) );
                }

                virtual parameter_type_t data_type() const override {
                    return parameter_type_t::STRING;
                }

                virtual std::vector<bool> asBoolVector() const override {
                    std::vector<bool> vec;
                    for (size_t i = 0; i < this->value().size(); ++i) {
                        vec[ i ] = ( this->value().at( i ) != "" );
                    }
                    return vec;
                }

                virtual std::vector<int> asIntVector() const override {
                    std::vector<int> vec;
                    for (size_t i = 0; i < this->value().size(); ++i) {
                        vec[ i ] = std::stoi( this->value().at( i ) );
                    }
                    return vec;
                }

                virtual std::vector<double> asDoubleVector() const override {
                    std::vector<double> vec;
                    for (size_t i = 0; i < this->value().size(); ++i) {
                        vec[ i ] = std::stod( this->value().at( i ) );
                    }
                    return vec;
                }

                virtual std::vector<std::string> asStringVector()
                    const override {
                    return this->value();
                }

                virtual std::string data_type_str() const override {
                    return "string";
                }

                virtual std::string json_value_str() const override {
                    std::ostringstream oss;
                    oss << "[ ";
                    for (auto it = this->value().begin();
                         it != this->value().end(); ++it) {
                        if (it != this->value().begin()) {
                            oss << ", ";
                        }
                        oss << "\"" << *it << "\"";
                    }
                    oss << " ]";
                    return oss.str();
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key,
                    const std::vector<std::string>& value,
                    const std::string& description ) {
                    return std::unique_ptr<Parameter>(
                        new StringVectorParameter( key, value, description ) );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key,
                    const std::vector<std::string>& value ) {
                    return std::unique_ptr<Parameter>(
                        new StringVectorParameter( key, value ) );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key,
                    const std::initializer_list<std::string>& value,
                    const std::string& description ) {
                    return std::unique_ptr<Parameter>(
                        new StringVectorParameter( key, value, description ) );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key,
                    const std::initializer_list<std::string>& value ) {
                    return std::unique_ptr<Parameter>(
                        new StringVectorParameter( key, value ) );
                }
        };
    }  // namespace processing
}  // namespace NCPA
