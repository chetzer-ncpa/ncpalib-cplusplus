#pragma once

#include "NCPA/atmosphere/Atmosphere1D.hpp"
#include "NCPA/atmosphere/abstract_atmosphere_1d.hpp"
#include "NCPA/atmosphere/stratified_atmosphere_1d.hpp"
#include "NCPA/strings.hpp"
#include "NCPA/units.hpp"

#include <fstream>
#include <iostream>
#include <string>

namespace NCPA {
    namespace atmos {
        class AtmosphereFactory {
            public:
                AtmosphereFactory() {}

                // Atmosphere1D build( const std::string& filename ) {
                //     return build( filename, "" );
                // }

                Atmosphere1D build( const std::string& filename,
                                    const std::string& headerfilename = "" ) {
                    std::string hfile( headerfilename );
                    if ( headerfilename.size() == 0 ) {
                        hfile = filename;
                    }
                    std::ifstream header_in( hfile );
                    read_header_from_stream( header_in );
                    header_in.close();

                    std::ifstream in( filename );
                    Atmosphere1D atm = read_values_from_stream( in );
                    in.close();
                    _headerlines.clear();
                    return atm;
                }

                void read_header_from_stream( std::istream& in ) {
                    std::string line;

                    std::getline( in, line );
                    while ( in.good() ) {
                        // lines will either be comments (# ), field
                        // descriptions (#% ), or field contents
                        line = NCPA::strings::deblank( line );
                        if ( line[ 0 ] == '#' ) {
                            // check second character
                            if ( line.size() > 1 && line[ 1 ] == '%' ) {
                                _headerlines.push_back( line.substr( 2 ) );
                            }  // otherwise it's a regular comment and can be
                               // ignored
                        }

                        std::getline( in, line );
                    }
                }

                Atmosphere1D read_values_from_stream( std::istream& in ) {
                    if ( !in.good() ) {
                        throw std::runtime_error(
                            "Atmosphere1D - Input stream not in good state" );
                    }

                    std::string line;
                    std::vector<std::string> atmlines, scalarlines;
                    std::ostringstream oss;  // for exceptions
                    size_t i;                // repeated index variable
                    details::stratified_atmosphere_1d atm;

                    std::getline( in, line );
                    while ( in.good() ) {
                        // lines will either be comments (# ), field
                        // descriptions (#% ), or field contents
                        line = NCPA::strings::deblank( line );
                        if ( line[ 0 ] == '#' ) {
                            // check second character
                            // if (line.size() > 1 && line[ 1 ] == '%') {
                            // 	_headerlines.push_back( line.substr( 2 ) );
                            // } // otherwise it's a regular comment and can be
                            // ignored

                        } else if ( line.size() == 0 ) {
                            // skip empty lines
                        } else {
                            atmlines.push_back( line );
                        }

                        std::getline( in, line );
                    }
                    // in.close();
                    // cout << "Found " << _headerlines.size() << " header
                    // lines" << endl; cout << "Found " << atmlines.size() << "
                    // data lines" << endl;

                    // parse them out
                    size_t nfields = _headerlines.size();
                    if ( nfields == 0 ) {
                        throw std::runtime_error(
                            "Atmosphere1D - No descriptive fields found." );
                    }

                    // hold contents
                    std::vector<std::string> keys, fields;
                    std::vector<size_t> column_numbers;
                    std::vector<double> values;
                    std::vector<units_ptr_t> units;

                    for ( auto it = _headerlines.cbegin();
                          it != _headerlines.cend(); ++it ) {
                        fields = NCPA::strings::split(
                            NCPA::strings::deblank( *it ), " ," );
                        if ( fields.size() != 3 && fields.size() != 4 ) {
                            oss << "Atmosphere1D - Error parsing descriptive "
                                   "line:"
                                << std::endl
                                << line << std::endl
                                << "Must be formatted as:" << std::endl
                                << "column,key,units[,value]" << std::endl
                                << "Use column=0 and specify value for scalar "
                                   "quantities."
                                << std::endl;
                            throw std::invalid_argument( oss.str() );
                        }

                        // process fields
                        // field line needs to be parseable as an integer
                        size_t col;
                        try {
                            col = (size_t)std::stoi( fields[ 0 ] );
                        } catch ( std::invalid_argument& e ) {
                            oss << "Atmosphere1D - Error parsing descriptive "
                                   "line:"
                                << std::endl
                                << line << std::endl
                                << "First field not parseable as an integer"
                                << std::endl;
                            throw std::invalid_argument( oss.str() );
                        }

                        double tempval = 0.0;
                        units_ptr_t tempunits
                            = NCPA::units::Units::from_string(
                                NCPA::strings::deblank( fields[ 2 ] ) );
                        if ( tempunits == nullptr ) {
                            throw std::invalid_argument(
                                "Unrecognized units string: " + fields[ 2 ] );
                        }
                        if ( fields.size() == 4 ) {
                            try {
                                tempval = std::stod(
                                    NCPA::strings::deblank( fields[ 3 ] ) );
                            } catch ( std::out_of_range& e ) {
                                oss << "Error converting data line:"
                                    << std::endl
                                    << *it << std::endl
                                    << "Value " << fields[ 3 ]
                                    << " is out of range for double "
                                       "precision.";
                                throw std::runtime_error( oss.str() );
                            }
                        }

                        // add to header vectors
                        column_numbers.push_back( col );
                        keys.push_back( NCPA::strings::deblank( fields[ 1 ] ) );
                        units.push_back( tempunits );
                        values.push_back( tempval );  // this will be ignored
                                                      // for vector quantities
                    }

                    // check for uniqueness in keys
                    std::vector<std::string> tempkeys( keys );
                    std::sort( tempkeys.begin(), tempkeys.end() );
                    std::vector<std::string>::iterator uit
                        = std::unique( tempkeys.begin(), tempkeys.end() );
                    if ( uit != tempkeys.end() ) {
                        throw std::invalid_argument(
                            "Atmosphere1D: key values are not unique" );
                    }


                    int nlines           = atmlines.size();
                    size_t ncols         = 0;
                    units_ptr_t depunits = nullptr;
                    for ( i = 0; i < column_numbers.size(); i++ ) {
                        ncols = column_numbers[ i ] > ncols
                                  ? column_numbers[ i ]
                                  : ncols;
                        if ( column_numbers[ i ] == 1 ) {
                            depunits = units[ i ];
                        }
                    }

                    std::vector<double *> columns( ncols );
                    for ( i = 0; i < ncols; i++ ) {
                        double *col  = new double[ nlines ];
                        columns[ i ] = col;
                    }


                    // step through the data lines
                    size_t row = 0;
                    std::vector<std::string>::const_iterator it;
                    for ( it = atmlines.cbegin(); it != atmlines.cend();
                          ++it ) {
                        fields.clear();
                        fields = NCPA::strings::split(
                            NCPA::strings::deblank( *it ), " \t," );
                        if ( fields.size() != ncols ) {
                            oss << "Error parsing data line:" << std::endl
                                << *it << std::endl
                                << ncols << " columns expected, "
                                << fields.size() << " columns found.";
                            for ( i = 0; i < ncols; i++ ) {
                                delete[] columns[ i ];
                            }
                            throw std::invalid_argument( oss.str() );
                        }
                        for ( i = 0; i < ncols; i++ ) {
                            try {
                                double *thiscol = columns[ i ];
                                thiscol[ row ]  = std::stod( fields[ i ] );
                            } catch ( std::out_of_range& e ) {
                                oss << "Error converting data line:"
                                    << std::endl
                                    << *it << std::endl
                                    << "Value " << fields[ i ]
                                    << " is out of range for double "
                                       "precision.";
                                for ( i = 0; i < ncols; i++ ) {
                                    delete[] columns[ i ];
                                }
                                throw std::runtime_error( oss.str() );
                            } catch ( std::invalid_argument& e ) {
                                oss << "Error parsing data line:" << std::endl
                                    << *it << std::endl
                                    << "Can't parse field " << fields[ i ]
                                    << " as a double";
                                for ( i = 0; i < ncols; i++ ) {
                                    delete[] columns[ i ];
                                }
                                throw std::runtime_error( oss.str() );
                            }
                        }
                        row++;
                    }

                    atm.set_altitude_vector(
                        vector_t( nlines, columns[ 0 ], depunits ) );
                    // z_ = new NCPA::VectorWithUnits( nlines, columns[ 0 ],
                    //                                 depunits );

                    for ( i = 0; i < keys.size(); i++ ) {
                        if ( column_numbers[ i ] == 0 ) {
                            atm.add_property(
                                keys[ i ],
                                scalar_t( values[ i ], *units[ i ] ) );
                        } else if ( column_numbers[ i ] > 1 ) {
                            atm.add_property(
                                keys[ i ],
                                vector_t( nlines,
                                          columns[ column_numbers[ i ] - 1 ],
                                          units[ i ] ) );
                        }
                    }

                    for ( i = 0; i < ncols; i++ ) {
                        delete[] columns[ i ];
                    }

                    return Atmosphere1D(
                        std::unique_ptr<details::abstract_atmosphere_1d>(
                            new details::stratified_atmosphere_1d( atm ) ) );
                }

            private:
                std::vector<std::string> _headerlines;
        };
    }  // namespace atmos
}  // namespace NCPA
