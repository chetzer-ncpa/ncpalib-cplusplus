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
        class DoubleVectorParameter : public VectorParameter<double> {
            public:
                DoubleVectorParameter( const std::string& key,
                                       const std::vector<double>& value,
                                       const std::string& description = "" ) :
                    VectorParameter<double>( key, value, description ) {}

                DoubleVectorParameter( const std::string& key,
                                       std::initializer_list<double>& value,
                                       const std::string& description = "" ) :
                    VectorParameter<double>( key, value, description ) {}

                virtual ~DoubleVectorParameter() {}

                friend void swap(
                    NCPA::processing::DoubleVectorParameter& a,
                    NCPA::processing::DoubleVectorParameter& b ) noexcept {
                    using std::swap;
                    swap(
                        dynamic_cast<
                            NCPA::processing::VectorParameter<double>&>( a ),
                        dynamic_cast<
                            NCPA::processing::VectorParameter<double>&>( b ) );
                }

                virtual std::unique_ptr<Parameter> clone() const override {
                    return std::unique_ptr<Parameter>(
                        new DoubleVectorParameter( this->key(), this->value(),
                                                   this->description() ) );
                }

                virtual parameter_type_t data_type() const override {
                    return parameter_type_t::FLOAT;
                }

                virtual std::vector<bool> asBoolVector() const override {
                    std::vector<bool> vec;
                    for (size_t i = 0; i < this->value().size(); ++i) {
                        vec[ i ] = ( this->value().at( i ) != 0.0 );
                    }
                    return vec;
                }

                virtual std::vector<int> asIntVector() const override {
                    std::vector<int> vec;
                    for (size_t i = 0; i < this->value().size(); ++i) {
                        vec[ i ] = static_cast<int>(
                            this->value().at( i )
                            + std::numeric_limits<double>::epsilon() );
                    }
                    return vec;
                }

                virtual std::vector<double> asDoubleVector() const override {
                    return this->value();
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
                    return "float";
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key, const std::vector<double>& value,
                    const std::string& description ) {
                    return std::unique_ptr<Parameter>(
                        new DoubleVectorParameter( key, value, description ) );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key,
                    const std::vector<double>& value ) {
                    return std::unique_ptr<Parameter>(
                        new DoubleVectorParameter( key, value ) );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key,
                    const std::initializer_list<double>& value,
                    const std::string& description ) {
                    return std::unique_ptr<Parameter>(
                        new DoubleVectorParameter( key, value, description ) );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key,
                    const std::initializer_list<double>& value ) {
                    return std::unique_ptr<Parameter>(
                        new DoubleVectorParameter( key, value ) );
                }
        };
    }  // namespace processing
}  // namespace NCPA
