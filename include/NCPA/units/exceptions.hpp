#pragma once

#include "NCPA/units/declarations.hpp"
#include "NCPA/exceptions.hpp"
#include <stdexcept>

namespace NCPA {
    namespace units {
        NCPA_DEFINE_EXCEPTION( InvalidConversion, std::out_of_range )
        NCPA_DEFINE_EXCEPTION( UnitsNotSet, std::logic_error )
        // class invalid_conversion : public std::out_of_range {
        //     public:
        //         invalid_conversion() :
        //             std::out_of_range( "Invalid conversion" ) {}

        //         invalid_conversion( const std::string& msg ) :
        //             std::out_of_range( msg ) {}

        //             invalid_conversion( const Unit *from, const Unit *to ) :
        //             std::out_of_range( "Invalid conversion from " + from->name()
        //                                + " to " + to->name() ) {

        //             }

        //         invalid_conversion( const Unit& from, const Unit& to ) :
        //             std::out_of_range( "Invalid conversion from " + from.name()
        //                                + " to " + to.name() ) {}

        //         invalid_conversion( const Unit *from, const Unit *to ) :
        //             invalid_conversion( *from, *to ) {}
        // };
    }
}