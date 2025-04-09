#pragma once

#ifdef _WIN32
const char file_separator = '\\';
#else
const char file_separator = '/';
#endif

#include "NCPA/atmosphere/Atmosphere1D.hpp"
#include "NCPA/atmosphere/Atmosphere2D.hpp"
#include "NCPA/atmosphere/Atmosphere3D.hpp"
// #include "NCPA/atmosphere/builders.hpp"
#include "NCPA/atmosphere/declarations.hpp"
#include "NCPA/files.hpp"

#include <fstream>
#include <iostream>
#include <memory>
#include <vector>

static void swap( NCPA::atmos::_abstract_atmosphere_reader_1d&,
                  NCPA::atmos::_abstract_atmosphere_reader_1d& ) noexcept;
static void swap( NCPA::atmos::ncpaprop_atmosphere_reader_1d&,
                  NCPA::atmos::ncpaprop_atmosphere_reader_1d& ) noexcept;
static void swap( NCPA::atmos::_abstract_atmosphere_reader_2d&,
                  NCPA::atmos::_abstract_atmosphere_reader_2d& ) noexcept;
static void swap( NCPA::atmos::ncpaprop_atmosphere_reader_2d&,
                  NCPA::atmos::ncpaprop_atmosphere_reader_2d& ) noexcept;
static void swap(
    NCPA::atmos::ncpaprop_atmosphere_reader_stratified_2d&,
    NCPA::atmos::ncpaprop_atmosphere_reader_stratified_2d& ) noexcept;
static void swap( NCPA::atmos::_abstract_atmosphere_reader_3d&,
                  NCPA::atmos::_abstract_atmosphere_reader_3d& ) noexcept;
static void swap( NCPA::atmos::ncpaprop_atmosphere_reader_3d&,
                  NCPA::atmos::ncpaprop_atmosphere_reader_3d& ) noexcept;
static void swap(
    NCPA::atmos::ncpaprop_atmosphere_reader_stratified_3d&,
    NCPA::atmos::ncpaprop_atmosphere_reader_stratified_3d& ) noexcept;

static void swap( NCPA::atmos::AtmosphereReader1D&,
                  NCPA::atmos::AtmosphereReader1D& ) noexcept;
static void swap( NCPA::atmos::AtmosphereReader2D&,
                  NCPA::atmos::AtmosphereReader2D& ) noexcept;
static void swap( NCPA::atmos::AtmosphereReader3D&,
                  NCPA::atmos::AtmosphereReader3D& ) noexcept;

namespace NCPA {
    namespace atmos {
        class _abstract_atmosphere_reader_1d {
            public:
                virtual ~_abstract_atmosphere_reader_1d() {}

                virtual std::unique_ptr<_abstract_atmosphere_reader_1d> clone()
                    const
                    = 0;

                virtual Atmosphere1D read( std::istream& in1 ) = 0;
                virtual Atmosphere1D read(
                    std::vector<std::istream *>& streams )
                    = 0;
                virtual Atmosphere1D read( const std::string& filename ) = 0;
                virtual Atmosphere1D read(
                    const std::vector<const std::string>& filenames )
                    = 0;
        };

        class ncpaprop_atmosphere_reader_1d
            : public _abstract_atmosphere_reader_1d {
            public:
                ncpaprop_atmosphere_reader_1d() :
                    _abstract_atmosphere_reader_1d() {}

                ncpaprop_atmosphere_reader_1d(
                    const ncpaprop_atmosphere_reader_1d& other ) :
                    ncpaprop_atmosphere_reader_1d() {}

                DECLARE_BOILERPLATE_METHODS( ncpaprop_atmosphere_reader_1d,
                                             _abstract_atmosphere_reader_1d )

                virtual Atmosphere1D read(
                    const std::string& filename ) override {
                    std::ifstream fs( filename, std::ios_base::in );
                    // std::cout << "Reading from single filename: " <<
                    // filename << std::endl;
                    return this->read( fs );
                }

                virtual Atmosphere1D read(
                    const std::vector<const std::string>& filenames )
                    override {
                    // std::cout << "Reading from vector of filenames" <<
                    // std::endl;
                    std::ifstream fs1, fs2;
                    std::vector<std::istream *> strs;
                    switch (filenames.size()) {
                        case 0:
                            throw std::logic_error(
                                "ncpaprop_atmosphere_reader_1d.read(): No "
                                "filenames provided!" );
                            break;
                        case 1:
                            fs1.open( filenames[ 0 ], std::ios_base::in );
                            return this->read( fs1 );
                            break;
                        default:
                            fs1.open( filenames[ 0 ], std::ios_base::in );
                            fs2.open( filenames[ 1 ], std::ios_base::in );
                            strs = { &fs1, &fs2 };
                            return this->read( strs );
                            break;
                    }
                }

                virtual Atmosphere1D read( std::istream& in1 ) override {
                    // std::cout << "Reading from single stream" << std::endl;
                    return _read_from_streams( &in1 );
                }

                virtual Atmosphere1D read(
                    std::vector<std::istream *>& streams ) override {
                    // std::cout << "Reading from vector of streams" <<
                    // std::endl;
                    switch (streams.size()) {
                        case 0:
                            throw std::logic_error(
                                "ncpaprop_atmosphere_reader_1d.read(): No "
                                "streams provided!" );
                            break;
                        case 1:
                            return _read_from_streams( streams[ 0 ] );
                            break;
                        default:
                            return _read_from_streams( streams[ 0 ],
                                                       streams[ 1 ] );
                            break;
                    }
                }

            protected:
                Atmosphere1D _read_from_streams( std::istream *in,
                                                 std::istream *header_in
                                                 = nullptr ) {
                    std::vector<std::string> headerlines;
                    if (header_in != nullptr) {
                        headerlines = _read_header_lines( *header_in );
                    }

                    if (!in->good()) {
                        throw std::runtime_error(
                            "ncpaprop_atmosphere_reader_1d._read_from_streams("
                            "): Input stream not in good state" );
                    }

                    std::string line;
                    std::vector<std::string> atmlines, scalarlines;
                    std::ostringstream oss;  // for exceptions
                    size_t i;                // repeated index variable
                    tuple_atmosphere_1d atm;

                    std::getline( *in, line );
                    while (in->good()) {
                        // lines will either be comments (# ), field
                        // descriptions (#% ), or field contents
                        line = NCPA::strings::deblank( line );
                        // std::cout << "Line is '" << line << "'" <<
                        // std::endl;
                        if (line[ 0 ] == '#') {
                            // check second character
                            if (header_in == nullptr && line.size() > 1
                                && line[ 1 ] == '%') {
                                headerlines.push_back( line.substr( 2 ) );
                                // std::cout << "Added '" << line << "' to
                                // headerlines" << std::endl;
                            }
                            // otherwise it's a regular comment and can be
                            // ignored

                        } else if (line.size() == 0) {
                            // skip empty lines
                        } else {
                            atmlines.push_back( line );
                        }

                        std::getline( *in, line );
                    }
                    // std::cout << "Found " << headerlines.size()
                    //           << " header lines " << std::endl;
                    // std::cout << " Found " << atmlines.size() << " data
                    // lines "
                    //           << std::endl;

                    // parse them out
                    size_t nfields = headerlines.size();
                    if (nfields == 0) {
                        throw std::runtime_error(
                            "ncpaprop_atmosphere_reader_1d._read_from_streams("
                            "): No descriptive fields found." );
                    }

                    // hold contents
                    std::vector<std::string> keys, fields;
                    std::vector<size_t> column_numbers;
                    std::vector<double> values;
                    std::vector<units_ptr_t> units;

                    for (auto it = headerlines.cbegin();
                         it != headerlines.cend(); ++it) {
                        fields = NCPA::strings::split(
                            NCPA::strings::deblank( *it ), " ," );
                        if (fields.size() != 3 && fields.size() != 4) {
                            oss << "ncpaprop_atmosphere_reader_1d._read_from_"
                                   "streams(): Error parsing descriptive line:"
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
                        } catch (std::invalid_argument& e) {
                            oss << "ncpaprop_atmosphere_reader_1d._read_from_"
                                   "streams(): Error parsing descriptive "
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
                        if (tempunits == nullptr) {
                            throw std::invalid_argument(
                                "ncpaprop_atmosphere_reader_1d._read_from_"
                                "streams(): Unrecognized units string: "
                                + fields[ 2 ] );
                        }
                        if (fields.size() == 4) {
                            try {
                                tempval = std::stod(
                                    NCPA::strings::deblank( fields[ 3 ] ) );
                            } catch (std::out_of_range& e) {
                                oss << "ncpaprop_atmosphere_reader_1d._read_"
                                       "from_"
                                       "streams(): Error converting data line:"
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
                        keys.push_back(
                            NCPA::strings::deblank( fields[ 1 ] ) );
                        units.push_back( tempunits );
                        values.push_back( tempval );  // this will be ignored
                                                      // for vector quantities
                    }

                    // check for uniqueness in keys
                    std::vector<std::string> tempkeys( keys );
                    std::sort( tempkeys.begin(), tempkeys.end() );
                    std::vector<std::string>::iterator uit
                        = std::unique( tempkeys.begin(), tempkeys.end() );
                    if (uit != tempkeys.end()) {
                        throw std::invalid_argument(
                            "ncpaprop_atmosphere_reader_1d._read_from_"
                            "streams(): key values are not unique" );
                    }


                    int nlines           = atmlines.size();
                    size_t ncols         = 0;
                    units_ptr_t depunits = nullptr;
                    for (i = 0; i < column_numbers.size(); i++) {
                        ncols = column_numbers[ i ] > ncols
                                  ? column_numbers[ i ]
                                  : ncols;
                        if (column_numbers[ i ] == 1) {
                            depunits = units[ i ];
                        }
                    }

                    std::vector<double *> columns( ncols );
                    for (i = 0; i < ncols; i++) {
                        double *col  = new double[ nlines ];
                        columns[ i ] = col;
                    }


                    // step through the data lines
                    size_t row = 0;
                    std::vector<std::string>::const_iterator it;
                    for (it = atmlines.cbegin(); it != atmlines.cend(); ++it) {
                        fields.clear();
                        fields = NCPA::strings::split(
                            NCPA::strings::deblank( *it ), " \t," );
                        if (fields.size() != ncols) {
                            oss << "ncpaprop_atmosphere_reader_1d._read_from_"
                                   "streams(): Error parsing data line:"
                                << std::endl
                                << *it << std::endl
                                << ncols << " columns expected, "
                                << fields.size() << " columns found.";
                            for (i = 0; i < ncols; i++) {
                                delete[] columns[ i ];
                            }
                            throw std::invalid_argument( oss.str() );
                        }
                        for (i = 0; i < ncols; i++) {
                            try {
                                double *thiscol = columns[ i ];
                                thiscol[ row ]  = std::stod( fields[ i ] );
                            } catch (std::out_of_range& e) {
                                oss << "ncpaprop_atmosphere_reader_1d._read_"
                                       "from_"
                                       "streams(): Error converting data line:"
                                    << std::endl
                                    << *it << std::endl
                                    << "Value " << fields[ i ]
                                    << " is out of range for double "
                                       "precision.";
                                for (i = 0; i < ncols; i++) {
                                    delete[] columns[ i ];
                                }
                                throw std::runtime_error( oss.str() );
                            } catch (std::invalid_argument& e) {
                                oss << "ncpaprop_atmosphere_reader_1d._read_"
                                       "from_"
                                       "streams(): Error parsing data line:"
                                    << std::endl
                                    << *it << std::endl
                                    << "Can't parse field " << fields[ i ]
                                    << " as a double";
                                for (i = 0; i < ncols; i++) {
                                    delete[] columns[ i ];
                                }
                                throw std::runtime_error( oss.str() );
                            }
                        }
                        row++;
                    }

                    atm.set_axis(
                        vector_u_t( nlines, columns[ 0 ], depunits ) );

                    for (i = 0; i < keys.size(); i++) {
                        if (column_numbers[ i ] == 0) {
                            atm.add_property(
                                keys[ i ],
                                scalar_u_t( values[ i ], *units[ i ] ) );
                        } else if (column_numbers[ i ] > 1) {
                            // AtmosphericProperty1D prop
                            //     = AtmosphereFactory::build(
                            //         atmospheric_property_1d_t::TUPLE );
                            tuple_atmospheric_property_1d prop;
                            prop.set(
                                atm.get_axis_vector(),
                                vector_u_t( nlines,
                                            columns[ column_numbers[ i ] - 1 ],
                                            units[ i ] ) );
                            atm.add_property( keys[ i ], prop );
                        }
                    }

                    for (i = 0; i < ncols; i++) {
                        delete[] columns[ i ];
                    }

                    // std::cout << "Read successful" << std::endl;

                    return Atmosphere1D(
                        std::unique_ptr<abstract_atmosphere_1d>(
                            new tuple_atmosphere_1d( atm ) ) );
                }

                std::vector<std::string> _read_header_lines(
                    std::istream& in ) {
                    std::string line;
                    std::vector<std::string> headerlines;
                    if (!in.good()) {
                        throw std::runtime_error(
                            "ncpaprop_atmosphere_reader_1d._read_header_lines("
                            "): Header input stream not in good "
                            "state" );
                    }

                    std::getline( in, line );
                    while (in.good()) {
                        // lines will either be comments (# ), field
                        // descriptions (#% ), or field contents
                        line = NCPA::strings::deblank( line );
                        if (line[ 0 ] == '#') {
                            // check second character
                            if (line.size() > 1 && line[ 1 ] == '%') {
                                headerlines.push_back( line.substr( 2 ) );
                            }  // otherwise it's a regular comment and can be
                               // ignored
                        }

                        std::getline( in, line );
                    }
                    return headerlines;
                }
        };

        class _abstract_atmosphere_reader_2d {
            public:
                virtual ~_abstract_atmosphere_reader_2d() {}

                virtual std::unique_ptr<_abstract_atmosphere_reader_2d> clone()
                    const
                    = 0;

                virtual Atmosphere2D read( const std::string& filename ) = 0;
                virtual Atmosphere2D read(
                    const std::vector<const std::string>& filenames )
                    = 0;
                virtual _abstract_atmosphere_reader_2d& set_axis_units(
                    units_ptr_t u )
                    = 0;
        };

        class ncpaprop_atmosphere_reader_2d
            : public _abstract_atmosphere_reader_2d {
            public:
                ncpaprop_atmosphere_reader_2d() :
                    _abstract_atmosphere_reader_2d(),
                    _axis_units { &NCPA::units::KILOMETERS } {}

                ncpaprop_atmosphere_reader_2d(
                    const ncpaprop_atmosphere_reader_2d& other ) :
                    ncpaprop_atmosphere_reader_2d() {
                    _axis_units = other._axis_units;
                }

                DECLARE_BOILERPLATE_METHODS( ncpaprop_atmosphere_reader_2d,
                                             _abstract_atmosphere_reader_2d )

                virtual bool stratified() const { return false; }

                virtual Atmosphere2D read(
                    const std::string& filename ) override {
                    // auto const pos = filename.find_last_of(
                    // NCPA::files::filesep() );
                    std::string basedir = NCPA::files::pathname( filename );
                    // if (pos == filename.npos) {
                    //     basedir = ".";
                    // } else {
                    //     basedir = filename.substr( 0, pos );
                    // }
                    std::ifstream summary( filename, std::ios_base::in );
                    vector_u_t ranges( 0, "km" );
                    std::vector<std::string> fns;
                    std::string line;
                    do {
                        std::getline( summary, line );
                        line = NCPA::strings::deblank( line );
                        if (line.length() > 0 && line[ 0 ] != '#') {
                            auto tokens = NCPA::strings::split( line );
                            if (tokens.size() != 2) {
                                std::ostringstream oss;
                                oss << "Atmosphere2D.read(): Malformed line:"
                                    << std::endl
                                    << "'" << line << "'" << std::endl
                                    << "Must be '<range> <filename>'";
                                throw std::invalid_argument( oss.str() );
                            }
                            ranges.push_back( std::stod( tokens[ 0 ] ) );
                            fns.push_back( NCPA::files::fullfile(
                                basedir, tokens[ 1 ] ) );
                        }
                    } while (!summary.eof());
                    summary.close();

                    std::vector<Atmosphere1D> atms;
                    ncpaprop_atmosphere_reader_1d reader;
                    for (auto it = fns.begin(); it != fns.end(); ++it) {
                        atms.push_back( reader.read( *it ) );
                    }
                    Atmosphere2D atm2d;
                    if (atms.size() == 1) {
                        atm2d = Atmosphere2D(
                            _atm_2d_ptr_t( new stratified_atmosphere_2d() ) );
                        atm2d.set( ranges, atms );
                    } else {
                        if (this->stratified()) {
                            atm2d = Atmosphere2D( _atm_2d_ptr_t(
                                new piecewise_stratified_atmosphere_2d() ) );
                        } else {
                            atm2d = Atmosphere2D(
                                _atm_2d_ptr_t( new grid_atmosphere_2d() ) );
                        }
                        atm2d.set( ranges, atms );
                    }
                    return atm2d;
                }

                virtual Atmosphere2D read(
                    const std::vector<const std::string>& filenames )
                    override {
                    if (filenames.size() != 2) {
                        throw std::range_error(
                            "Filename vector must have 2 elements: summary "
                            "file, header file" );
                    }
                    std::ifstream summary( filenames[ 0 ], std::ios_base::in );
                    vector_u_t ranges( 0, _axis_units );
                    std::vector<std::string> fns;
                    std::string basedir
                        = NCPA::files::pathname( filenames[ 1 ] );
                    std::string line;
                    do {
                        std::getline( summary, line );
                        line = NCPA::strings::deblank( line );
                        if (line.length() > 0 && line[ 0 ] != '#') {
                            auto tokens = NCPA::strings::split( line );
                            if (tokens.size() != 2) {
                                std::ostringstream oss;
                                oss << "ncpaprop_atmosphere_reader_2d.read(): "
                                       "Malformed line:"
                                    << std::endl
                                    << "'" << line << "'" << std::endl
                                    << "Must be '<range> <filename>'";
                                throw std::invalid_argument( oss.str() );
                            }
                            ranges.push_back( std::stod( tokens[ 0 ] ) );
                            fns.push_back( NCPA::files::fullfile(
                                basedir, tokens[ 1 ] ) );
                        }
                    } while (!summary.eof());
                    summary.close();

                    std::vector<Atmosphere1D> atms;
                    ncpaprop_atmosphere_reader_1d reader;
                    for (auto it = fns.begin(); it != fns.end(); ++it) {
                        std::vector<const std::string> fns { *it,
                                                             filenames[ 1 ] };
                        atms.push_back( reader.read( fns ) );
                    }
                    Atmosphere2D atm2d;
                    if (atms.size() == 1) {
                        atm2d = Atmosphere2D(
                            _atm_2d_ptr_t( new stratified_atmosphere_2d() ) );
                        atm2d.set( ranges, atms );
                    } else {
                        if (this->stratified()) {
                            atm2d = Atmosphere2D( _atm_2d_ptr_t(
                                new piecewise_stratified_atmosphere_2d() ) );
                        } else {
                            atm2d = Atmosphere2D(
                                _atm_2d_ptr_t( new grid_atmosphere_2d() ) );
                        }
                        atm2d.set( ranges, atms );
                    }
                    return atm2d;
                }

                virtual _abstract_atmosphere_reader_2d& set_axis_units(
                    units_ptr_t u ) override {
                    _axis_units = u;
                    return static_cast<_abstract_atmosphere_reader_2d&>(
                        *this );
                }

            protected:
                units_ptr_t _axis_units;
        };

        class ncpaprop_atmosphere_reader_stratified_2d
            : public ncpaprop_atmosphere_reader_2d {
            public:
                ncpaprop_atmosphere_reader_stratified_2d() :
                    ncpaprop_atmosphere_reader_2d() {}

                ncpaprop_atmosphere_reader_stratified_2d(
                    const ncpaprop_atmosphere_reader_stratified_2d& other ) :
                    ncpaprop_atmosphere_reader_2d() {}

                DECLARE_BOILERPLATE_METHODS(
                    ncpaprop_atmosphere_reader_stratified_2d,
                    _abstract_atmosphere_reader_2d )

                virtual bool stratified() const override { return true; }
        };

        class _abstract_atmosphere_reader_3d {
            public:
                virtual ~_abstract_atmosphere_reader_3d() {}

                virtual std::unique_ptr<_abstract_atmosphere_reader_3d> clone()
                    const
                    = 0;

                virtual Atmosphere3D read( const std::string& filename ) = 0;
                virtual Atmosphere3D read(
                    const std::vector<const std::string>& filenames )
                    = 0;
                virtual _abstract_atmosphere_reader_3d& set_axis_units(
                    size_t n, units_ptr_t u )
                    = 0;
        };

        class ncpaprop_atmosphere_reader_3d
            : public _abstract_atmosphere_reader_3d {
            public:
                ncpaprop_atmosphere_reader_3d() :
                    _abstract_atmosphere_reader_3d(),
                    _axis_units {
                        { &NCPA::units::KILOMETERS, &NCPA::units::KILOMETERS }
                } {}

                ncpaprop_atmosphere_reader_3d(
                    const ncpaprop_atmosphere_reader_3d& other ) :
                    ncpaprop_atmosphere_reader_3d() {
                    _axis_units = other._axis_units;
                }

                DECLARE_BOILERPLATE_METHODS( ncpaprop_atmosphere_reader_3d,
                                             _abstract_atmosphere_reader_3d )

                virtual bool stratified() const { return false; }

                virtual Atmosphere3D read(
                    const std::string& filename ) override {
                    vector_u_t x1( 0, _axis_units[ 0 ] ),
                        x2( 0, _axis_units[ 1 ] );
                    std::vector<size_t> x1inds, x2inds;
                    NCPA::arrays::ndvector<2, std::string> fns;
                    this->_parse_summary_file( filename, x1inds, x1, x2inds,
                                               x2, fns );
                    // std::cout << "Done parsing summary file" << std::endl;
                    NCPA::arrays::ndvector<2, Atmosphere1D> atms;
                    atms.reshape( fns.shape() );
                    // std::cout << "Created atmosphere 2d matrix of shape {"
                    //           << atms.shape()[ 0 ] << ", " << atms.shape()[
                    //           1 ]
                    //           << "}" << std::endl;
                    ncpaprop_atmosphere_reader_1d reader;
                    vector_u_t ux1 = x1, ux2 = x2;
                    ux1.resize( atms.shape()[ 0 ] );
                    ux2.resize( atms.shape()[ 1 ] );
                    for (size_t n = 0; n < x1.size(); ++n) {
                        Atmosphere1D atm1d
                            = reader.read( fns[ x1inds[ n ] ][ x2inds[ n ] ] );
                        atms[ x1inds[ n ] ][ x2inds[ n ] ] = atm1d;
                        ux1[ x1inds[ n ] ]                 = x1[ n ];
                        ux2[ x2inds[ n ] ]                 = x2[ n ];
                    }
                    Atmosphere3D atm3d;
                    if (x1.size() == 1 && x2.size() == 1) {
                        atm3d = Atmosphere3D(
                            _atm_3d_ptr_t( new stratified_atmosphere_3d() ) );
                        atm3d.set( atms[ 0 ][ 0 ] );
                    } else {
                        atm3d = Atmosphere3D(
                            _atm_3d_ptr_t( new grid_atmosphere_3d() ) );
                        atm3d.set( ux1, ux2, atms );
                    }
                    return atm3d;
                }

                virtual Atmosphere3D read(
                    const std::vector<const std::string>& filenames )
                    override {
                    if (filenames.size() != 2) {
                        throw std::range_error(
                            "Filename vector must have 2 elements: summary "
                            "file, header file" );
                    }
                    vector_u_t x1( 0, _axis_units[ 0 ] ),
                        x2( 0, _axis_units[ 1 ] );
                    std::vector<size_t> x1inds, x2inds;
                    NCPA::arrays::ndvector<2, std::string> fns;
                    this->_parse_summary_file( filenames[ 0 ], x1inds, x1,
                                               x2inds, x2, fns );
                    NCPA::arrays::ndvector<2, Atmosphere1D> atms;
                    atms.reshape( fns.shape() );
                    ncpaprop_atmosphere_reader_1d reader;
                    vector_u_t ux1 = x1, ux2 = x2;
                    ux1.resize( atms.shape()[ 0 ] );
                    ux2.resize( atms.shape()[ 1 ] );
                    for (size_t n = 0; n < fns.size(); ++n) {
                        atms[ x1inds[ n ] ][ x2inds[ n ] ]
                            = reader.read( { fns[ x1inds[ n ] ][ x2inds[ n ] ],
                                             filenames[ 1 ] } );
                        ux1[ x1inds[ n ] ] = x1[ n ];
                        ux2[ x2inds[ n ] ] = x2[ n ];
                    }
                    Atmosphere3D atm3d;
                    if (x1.size() && x2.size() == 1) {
                        atm3d = Atmosphere3D(
                            _atm_3d_ptr_t( new stratified_atmosphere_3d() ) );
                        atm3d.set( atms[ 0 ][ 0 ] );
                    } else {
                        atm3d = Atmosphere3D(
                            _atm_3d_ptr_t( new grid_atmosphere_3d() ) );
                        atm3d.set( ux1, ux2, atms );
                    }
                    return atm3d;
                }

                virtual _abstract_atmosphere_reader_3d& set_axis_units(
                    size_t n, units_ptr_t u ) override {
                    _axis_units[ n ] = u;
                    return static_cast<_abstract_atmosphere_reader_3d&>(
                        *this );
                }

            protected:
                std::array<units_ptr_t, 2> _axis_units;

                void _parse_summary_file(
                    const std::string& filename, std::vector<size_t>& x1inds,
                    vector_u_t& x1, std::vector<size_t>& x2inds,
                    vector_u_t& x2,
                    NCPA::arrays::ndvector<2, std::string>& fns ) {
                    std::string basedir = NCPA::files::pathname( filename );
                    std::string line;
                    std::ifstream summary( filename, std::ios_base::in );
                    std::vector<std::string> tmpfilenames;
                    do {
                        std::getline( summary, line );
                        line = NCPA::strings::deblank( line );
                        if (line.length() > 0 && line[ 0 ] != '#') {
                            auto tokens = NCPA::strings::split( line );
                            // std::cout << "Parsed line into { ";
                            // NCPA::arrays::write( std::cout, tokens );
                            // std::cout << " }" << std::endl;
                            if (tokens.size() != 5) {
                                std::ostringstream oss;
                                oss << "ncpaprop_atmosphere_reader_3d.read(): "
                                       "Malformed line:"
                                    << std::endl
                                    << line << std::endl
                                    << "Must be '<x1 index> <x1> <x2 index> "
                                       "<x2> <filename>'";
                                throw std::invalid_argument( oss.str() );
                            }
                            x1inds.push_back( std::stoi( tokens[ 0 ] ) );
                            x1.push_back( std::stod( tokens[ 1 ] ) );
                            x2inds.push_back( std::stoi( tokens[ 2 ] ) );
                            x2.push_back( std::stod( tokens[ 3 ] ) );
                            tmpfilenames.push_back( NCPA::files::fullfile(
                                basedir, tokens[ 4 ] ) );
                        }
                    } while (!summary.eof());
                    summary.close();

                    std::vector<size_t> ux1inds = x1inds, ux2inds = x2inds;
                    std::sort( ux1inds.begin(), ux1inds.end() );
                    std::sort( ux2inds.begin(), ux2inds.end() );
                    size_t nx1 = std::distance(
                        ux1inds.begin(),
                        std::unique( ux1inds.begin(), ux1inds.end() ) );
                    size_t nx2 = std::distance(
                        ux2inds.begin(),
                        std::unique( ux2inds.begin(), ux2inds.end() ) );
                    fns.reshape( { nx1, nx2 } );
                    for (size_t i = 0; i < tmpfilenames.size(); ++i) {
                        fns[ x1inds[ i ] ][ x2inds[ i ] ] = tmpfilenames[ i ];
                    }
                }
        };

        class ncpaprop_atmosphere_reader_stratified_3d
            : public ncpaprop_atmosphere_reader_3d {
            public:
                ncpaprop_atmosphere_reader_stratified_3d() :
                    ncpaprop_atmosphere_reader_3d() {}

                ncpaprop_atmosphere_reader_stratified_3d(
                    const ncpaprop_atmosphere_reader_stratified_3d& other ) :
                    ncpaprop_atmosphere_reader_stratified_3d() {}

                DECLARE_BOILERPLATE_METHODS(
                    ncpaprop_atmosphere_reader_stratified_3d,
                    _abstract_atmosphere_reader_3d )

                virtual bool stratified() const override { return true; }
        };

        // PUBLIC API
        class AtmosphereReader1D {
            public:
                AtmosphereReader1D() {}

                AtmosphereReader1D(
                    std::unique_ptr<_abstract_atmosphere_reader_1d> ptr ) :
                    AtmosphereReader1D() {
                    _ptr = std::move( ptr );
                }

                // copy constructor
                AtmosphereReader1D( const AtmosphereReader1D& other ) :
                    AtmosphereReader1D() {
                    _ptr = std::move( other._ptr->clone() );
                }

                DECLARE_WRAPPER_BOILERPLATE_METHODS( AtmosphereReader1D )

                virtual Atmosphere1D read( std::istream& in1 ) {
                    _check_pointer();
                    return _ptr->read( in1 );
                }

                virtual Atmosphere1D read(
                    std::vector<std::istream *>& streams ) {
                    _check_pointer();
                    return _ptr->read( streams );
                }

                virtual Atmosphere1D read( const std::string& filename ) {
                    _check_pointer();
                    return _ptr->read( filename );
                }

                virtual Atmosphere1D read(
                    const std::vector<const std::string>& filenames ) {
                    _check_pointer();
                    return _ptr->read( filenames );
                }

            protected:
                void _check_pointer() const {
                    if (!_ptr) {
                        throw std::logic_error(
                            "AtmosphereReader1D: No pointer set!" );
                    }
                }

                std::unique_ptr<_abstract_atmosphere_reader_1d> _ptr;
        };

        class AtmosphereReader2D {
            public:
                AtmosphereReader2D() {}

                AtmosphereReader2D(
                    std::unique_ptr<_abstract_atmosphere_reader_2d> ptr ) :
                    AtmosphereReader2D() {
                    _ptr = std::move( ptr );
                }

                // copy constructor
                AtmosphereReader2D( const AtmosphereReader2D& other ) :
                    AtmosphereReader2D() {
                    _ptr = std::move( other._ptr->clone() );
                }

                DECLARE_WRAPPER_BOILERPLATE_METHODS( AtmosphereReader2D )

                virtual Atmosphere2D read( const std::string& filename,
                                           bool stratified = false ) {
                    _check_pointer();
                    return _ptr->read( filename );
                }

                virtual Atmosphere2D read(
                    const std::vector<const std::string>& filenames,
                    bool stratified = false ) {
                    _check_pointer();
                    return _ptr->read( filenames );
                }

            protected:
                void _check_pointer() const {
                    if (!_ptr) {
                        throw std::logic_error(
                            "AtmosphereReader2D: No pointer set!" );
                    }
                }

                std::unique_ptr<_abstract_atmosphere_reader_2d> _ptr;
        };

        class AtmosphereReader3D {
            public:
                AtmosphereReader3D() {}

                AtmosphereReader3D(
                    std::unique_ptr<_abstract_atmosphere_reader_3d> ptr ) :
                    AtmosphereReader3D() {
                    _ptr = std::move( ptr );
                }

                // copy constructor
                AtmosphereReader3D( const AtmosphereReader3D& other ) :
                    AtmosphereReader3D() {
                    _ptr = std::move( other._ptr->clone() );
                }

                DECLARE_WRAPPER_BOILERPLATE_METHODS( AtmosphereReader3D )

                virtual Atmosphere3D read( const std::string& filename,
                                           bool stratified = false ) {
                    _check_pointer();
                    return _ptr->read( filename );
                }

                virtual Atmosphere3D read(
                    const std::vector<const std::string>& filenames,
                    bool stratified = false ) {
                    _check_pointer();
                    return _ptr->read( filenames );
                }

                virtual AtmosphereReader3D& set_axis_units( size_t n,
                                                            units_ptr_t u ) {
                    _check_pointer();
                    _ptr->set_axis_units( n, u );
                    return *this;
                }

            protected:
                void _check_pointer() const {
                    if (!_ptr) {
                        throw std::logic_error(
                            "AtmosphereReader3D: No pointer set!" );
                    }
                }

                std::unique_ptr<_abstract_atmosphere_reader_3d> _ptr;
        };


    }  // namespace atmos
}  // namespace NCPA

static void swap( NCPA::atmos::_abstract_atmosphere_reader_1d& a,
                  NCPA::atmos::_abstract_atmosphere_reader_1d& b ) noexcept {}

static void swap( NCPA::atmos::ncpaprop_atmosphere_reader_1d& a,
                  NCPA::atmos::ncpaprop_atmosphere_reader_1d& b ) noexcept {
    ::swap( dynamic_cast<NCPA::atmos::_abstract_atmosphere_reader_1d&>( a ),
            dynamic_cast<NCPA::atmos::_abstract_atmosphere_reader_1d&>( b ) );
}

static void swap( NCPA::atmos::_abstract_atmosphere_reader_2d& a,
                  NCPA::atmos::_abstract_atmosphere_reader_2d& b ) noexcept {}

static void swap( NCPA::atmos::ncpaprop_atmosphere_reader_2d& a,
                  NCPA::atmos::ncpaprop_atmosphere_reader_2d& b ) noexcept {
    ::swap( dynamic_cast<NCPA::atmos::_abstract_atmosphere_reader_2d&>( a ),
            dynamic_cast<NCPA::atmos::_abstract_atmosphere_reader_2d&>( b ) );
}

static void swap(
    NCPA::atmos::ncpaprop_atmosphere_reader_stratified_2d& a,
    NCPA::atmos::ncpaprop_atmosphere_reader_stratified_2d& b ) noexcept {
    ::swap( dynamic_cast<NCPA::atmos::ncpaprop_atmosphere_reader_2d&>( a ),
            dynamic_cast<NCPA::atmos::ncpaprop_atmosphere_reader_2d&>( b ) );
}

static void swap( NCPA::atmos::_abstract_atmosphere_reader_3d& a,
                  NCPA::atmos::_abstract_atmosphere_reader_3d& b ) noexcept {}

static void swap( NCPA::atmos::ncpaprop_atmosphere_reader_3d& a,
                  NCPA::atmos::ncpaprop_atmosphere_reader_3d& b ) noexcept {
    ::swap( dynamic_cast<NCPA::atmos::_abstract_atmosphere_reader_3d&>( a ),
            dynamic_cast<NCPA::atmos::_abstract_atmosphere_reader_3d&>( b ) );
}

static void swap(
    NCPA::atmos::ncpaprop_atmosphere_reader_stratified_3d& a,
    NCPA::atmos::ncpaprop_atmosphere_reader_stratified_3d& b ) noexcept {
    ::swap( dynamic_cast<NCPA::atmos::ncpaprop_atmosphere_reader_3d&>( a ),
            dynamic_cast<NCPA::atmos::ncpaprop_atmosphere_reader_3d&>( b ) );
}

static void swap( NCPA::atmos::AtmosphereReader1D& a,
                  NCPA::atmos::AtmosphereReader1D& b ) noexcept {
    using std::swap;
    std::swap( a._ptr, b._ptr );
}

static void swap( NCPA::atmos::AtmosphereReader2D& a,
                  NCPA::atmos::AtmosphereReader2D& b ) noexcept {
    using std::swap;
    std::swap( a._ptr, b._ptr );
}

static void swap( NCPA::atmos::AtmosphereReader3D& a,
                  NCPA::atmos::AtmosphereReader3D& b ) noexcept {
    using std::swap;
    std::swap( a._ptr, b._ptr );
}
