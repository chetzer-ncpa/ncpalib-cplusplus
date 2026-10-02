#pragma once

#ifdef HAVE_NCPA_JSON
#  undef HAVE_NCPA_JSON
#endif

#if __has_include( "nlohmann/json.hpp" )
#  define HAVE_NCPA_JSON true

#  include "NCPA/cloneable.hpp"
#  include "nlohmann/json.hpp"

#  include <map>
#  include <string>

namespace NCPA {
    namespace json {
        using JSON = nlohmann::json;

        enum class json_form_t { SCALAR = 0, VECTOR };

        enum class json_type_t { INTEGER, FLOAT, STRING, BOOLEAN };

        std::string to_string( json_form_t t ) {
            switch (t) {
                case json_form_t::SCALAR:
                    return "scalar";
                case json_form_t::VECTOR:
                    return "vector";
                default:
                    throw std::out_of_range( "unsupported JSON form" );
            }
        }

        std::string to_string( json_type_t t ) {
            switch (t) {
                case json_type_t::INTEGER:
                    return "integer";
                case json_type_t::FLOAT:
                    return "float";
                case json_type_t::STRING:
                    return "string";
                case json_type_t::BOOLEAN:
                    return "boolean";
                default:
                    throw std::out_of_range( "unsupported JSON type" );
            }
        }

        class JSONItem : public nlohmann::json {
            public:
                template<typename T>
                JSONItem( const std::string& key, const T& val,
                          json_form_t valform, json_type_t valtype,
                          const std::string& comment = "" ) :
                    nlohmann::json() {
                    ( *this )[ key ][ "form" ]    = to_string( valform );
                    ( *this )[ key ][ "type" ]    = to_string( valtype );
                    ( *this )[ key ][ "value" ]   = val;
                    ( *this )[ key ][ "comment" ] = comment;
                }

                virtual ~JSONItem() = default;

                JSONItem( const JSONItem& other ) = default;

                JSONItem( JSONItem&& other ) noexcept = default;
        };
    }  // namespace json
}  // namespace NCPA
#else
#  define HAVE_NCPA_JSON false
#endif
