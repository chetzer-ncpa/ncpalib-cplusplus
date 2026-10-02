#pragma once

#include "NCPA.hpp"
#include "NCPA/processing/declarations.hpp"
#include "NCPA/processing/parameters/declarations.hpp"
#include "NCPA/processing/parameters/Parameter.hpp"
#include "NCPA/processing/parameters/ScalarParameter.hpp"

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
        class StringParameter
            : public Cloneable<StringParameter, ScalarParameter<std::string>> {
            public:
                StringParameter() = default;

                StringParameter( const std::string& key, std::string value,
                                 const std::string& description = "" ) :
                    cloneable_t( key, value, description ) {}

                virtual ~StringParameter() {}

                friend void swap( StringParameter& a,
                                  StringParameter& b ) noexcept {
                    using std::swap;
                    swap( static_cast<cloneable_t&>( a ),
                          static_cast<cloneable_t&>( b ) );
                }

                virtual Parameter& clear() override {
                    this->set( "" );
                    return *this;
                }

                // virtual std::unique_ptr<Parameter> clone() const override {
                //     return std::unique_ptr<Parameter>(
                //         new StringParameter( this->key(), this->value() ) );
                // }

                virtual parameter_type_t data_type() const override {
                    return parameter_type_t::STRING;
                }

                virtual bool asBool() const override {
                    return ( this->value() != "" );
                }

                virtual int asInt() const override {
                    return std::stoi( this->value() );
                }

                virtual double asDouble() const override {
                    return std::stod( this->value() );
                }

                virtual std::string asString() const override {
                    return this->value();
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key, const std::string& value,
                    const std::string& description ) {
                    return std::unique_ptr<Parameter>(
                        new StringParameter( key, value, description ) );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key, const std::string& value ) {
                    return std::unique_ptr<Parameter>(
                        new StringParameter( key, value ) );
                }

                virtual std::string data_type_str() const override {
                    return "string";
                }

                virtual std::string json_value_str() const override {
                    return "\"" + this->asString() + "\"";
                }

#if HAVE_NLOHMANN_JSON_HPP
                virtual void from_json( nlohmann::json& json ) override {
                    this->set( json.at( "value" ).get<std::string>() );
                }
#endif
        };
    }  // namespace processing
}  // namespace NCPA
