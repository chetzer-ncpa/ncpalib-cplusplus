#pragma once

#include "NCPA/processing/declarations.hpp"
#include "NCPA/processing/parameters/declarations.hpp"
#include "NCPA/processing/parameters/Parameter.hpp"

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
        template<typename T>
        class ScalarParameter : public Parameter {
            public:
                ScalarParameter( const std::string& key, const T& value,
                                 const std::string& description = "" ) :
                    Parameter( key, description ), _value { value } {}

                virtual ~ScalarParameter() {}

                friend void swap(
                    NCPA::processing::ScalarParameter<T>& a,
                    NCPA::processing::ScalarParameter<T>& b ) noexcept {
                    using std::swap;
                    swap( dynamic_cast<NCPA::processing::Parameter&>( a ),
                          dynamic_cast<NCPA::processing::Parameter&>( b ) );
                }

                virtual const T& value() const { return _value; }

                virtual ScalarParameter& set( const T& val ) {
                    _value = val;
                    return *this;
                }

#if HAVE_NLOHMANN_JSON_HPP
                virtual void add_value_to_json(
                    nlohmann::json& json ) const override {
                    json[ "value" ] = this->value();
                }
#endif

                virtual bool is_scalar() const override { return true; }

                virtual std::vector<bool> asBoolVector() const override {
                    return std::vector<bool> { this->asBool() };
                }

                virtual std::vector<int> asIntVector() const override {
                    return std::vector<int> { this->asInt() };
                }

                virtual std::vector<size_t> asSizeTVector() const override {
                    return std::vector<size_t> { this->asSizeT() };
                }

                virtual std::vector<double> asDoubleVector() const override {
                    return std::vector<double> { this->asDouble() };
                }

                virtual std::vector<std::string> asStringVector()
                    const override {
                    return std::vector<std::string> { this->asString() };
                }

                virtual parameter_form_t data_form() const override {
                    return parameter_form_t::SCALAR;
                }

                virtual std::string data_form_str() const override {
                    return "scalar";
                }

                virtual std::string json_value_str() const override {
                    return this->asString();
                }

            private:
                T _value;
        };

        template<typename T,
                 typename
                 = typename std::enable_if<std::is_scalar<T>::value>::type>
        static std::unique_ptr<Parameter> make_parameter(
            const std::string& key, const T& val ) {
            return std::unique_ptr<Parameter>(
                new ScalarParameter<T>( key, val ) );
        }
    }  // namespace processing
}  // namespace NCPA
