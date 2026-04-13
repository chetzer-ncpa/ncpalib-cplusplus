#pragma once

#include "NCPA/configuration/declarations.hpp"

#include <string>
#include <vector>

namespace NCPA {
    namespace config {
        class Argument {
            public:
                Argument() : Argument( "", "" ) {}

                Argument( const std::string& tag ) : Argument( tag, "" ) {}

                Argument( const std::string& tag,
                          const std::string helptext ) :
                    _tag { tag }, _help_text { helptext } {}

                virtual std::string help_text() const { return _help_text; }

                virtual Argument& set_help_text( const std::string& text ) {
                    _help_text = text;
                    return *this;
                }

                virtual Argument& set_tag( const std::string& tag ) {
                    _tag = tag;
                    return *this;
                }

                virtual std::string tag() const { return _tag; }

                virtual bool expects_value() const = 0;
                virtual parse_status_t parse( std::vector<std::string>& args )
                    = 0;

            protected:
                std::string _tag;
                std::string _help_text;
        };
    }  // namespace config
}  // namespace NCPA
