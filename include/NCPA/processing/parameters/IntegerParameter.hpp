#pragma once

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
        class IntegerParameter
            : public Cloneable<IntegerParameter, ScalarParameter<int>> {
            public:
                IntegerParameter() = default;

                IntegerParameter( const std::string& key, int value,
                                  const std::string& description = "" ) :
                    cloneable_t( key, value, description ) {}

                virtual ~IntegerParameter() {}

                friend void swap( IntegerParameter& a,
                                  IntegerParameter& b ) noexcept {
                    using std::swap;
                    swap( static_cast<cloneable_t&>( a ),
                          static_cast<cloneable_t&>( b ) );
                }

                virtual Parameter& clear() override {
                    this->set( 0 );
                    return *this;
                }

                // virtual std::unique_ptr<Parameter> clone() const override {
                //     return std::unique_ptr<Parameter>(
                //         new IntegerParameter( this->key(), this->value() )
                //         );
                // }

                virtual parameter_type_t data_type() const override {
                    return parameter_type_t::INTEGER;
                }

                virtual bool asBool() const override {
                    return ( this->value() != 0 );
                }

                virtual int asInt() const override { return this->value(); }

                virtual double asDouble() const override {
                    return static_cast<int>( this->value() );
                }

                virtual std::string asString() const override {
                    std::ostringstream oss;
                    oss << this->value();
                    return oss.str();
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key, int value,
                    const std::string& description ) {
                    return std::unique_ptr<Parameter>(
                        new IntegerParameter( key, value, description ) );
                }

                static std::unique_ptr<Parameter> build(
                    const std::string& key, int value ) {
                    return std::unique_ptr<Parameter>(
                        new IntegerParameter( key, value ) );
                }

                virtual std::string data_type_str() const override {
                    return "integer";
                }

#if HAVE_NLOHMANN_JSON_HPP
                virtual void from_json( nlohmann::json& json ) override {
                    this->set( json.at( "value" ).get<int>() );
                }
#endif
        };
    }  // namespace processing
}  // namespace NCPA
