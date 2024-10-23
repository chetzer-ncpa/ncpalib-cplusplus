

#pragma once

#include <string>
#include <vector>
#include <ctime>
#include <cmath>
#include <cstring>
#include <sstream>
#include <algorithm>

namespace NCPA::strings {

	/**
	Removes leading and trailing whitespace from a string, using
	specified whitespace characters.
	@brief Removes leading and trailing whitespace from a string.
	@param orig The original string to deblank.
	@param whitespace A string containing all characters to treat as whitespace
	@returns The original string with leading and trailing whitespace removed.
	*/
	std::string deblank( const std::string& str, const std::string& whitespace = " \t\n\r" ) {
		const size_t strBegin = str.find_first_not_of( whitespace );
		if (strBegin == std::string::npos) {
			return "";
		}

		const size_t strEnd = str.find_last_not_of( whitespace );
		const size_t strRange = strEnd - strBegin + 1;

		return str.substr( strBegin, strRange );
	}

	/**
	Splits a string into tokens using supplied delimiter characters.
	Similar to the Perl split() function.
	@brief Breaks a string into tokens using delimiter characters.
	@param input The string to split.
	@param delimiters The delimiter characters.  Defaults to whitespace.
	@returns A vector of the string tokens.
	 */
	std::vector< std::string > split( std::string input, std::string delimiters = " \t\n\r" ) {
		std::vector< std::string > tokens;
		tokens.clear();
		size_t firstind, lastind;
		firstind = input.find_first_not_of( delimiters );
		lastind = input.find_first_of( delimiters, firstind );
		while( firstind != std::string::npos ) {
			tokens.push_back( input.substr( firstind, lastind - firstind ) );
			firstind = input.find_first_not_of( delimiters, lastind );
			lastind = input.find_first_of( delimiters, firstind );
		}
		return tokens;
	}

	/**
	Converts a numerical time to its equivalent time string in the
	format yyyy-mm-dd hh:mm:ss.mmm
	@brief Returns a human-readable time from an epoch time.
	@param d The time to print, as an epoch time_t value.
	@returns The string representation of the time.
	*/
	std::string timeAsString( double d ) {
	    time_t temptime = (time_t)d;
	    tm* uttime = std::gmtime( &temptime );
	    double ipart, fpart;
	    fpart = std::modf( d, &ipart);
	    char holder[50];

	    std::sprintf(holder,"%4d-%02d-%02d %02d:%02d:%06.3f GMT",
	            (unsigned char)(uttime->tm_year)+1900,
				(unsigned char)(uttime->tm_mon+1),
				(unsigned char)(uttime->tm_mday),
	            (unsigned char)(uttime->tm_hour),
				(unsigned char)(uttime->tm_min),
				(double)(uttime->tm_sec) + fpart );
	    std::string s = holder;
	    return s;
	}

	// Taken from https://en.cppreference.com/w/cpp/string/byte/tolower
	std::string to_lower( std::string s ) {
		std::transform( s.begin(), s.end(), s.begin(),
						[](unsigned char c){ return std::tolower(c); }
		);
		return s;
	}

	// Taken from https://en.cppreference.com/w/cpp/string/byte/tolower
	std::string to_upper( std::string s ) {
		std::transform( s.begin(), s.end(), s.begin(),
						[](unsigned char c){ return std::toupper(c); }
		);
		return s;
	}

//	/**
//	Converts the supplied string to lower case.
//	@brief Converts to lower case.
//	@param in The string to convert.
//	@returns The input string converted into all lower case.
//	*/
//	void toLowerCase( char *in ) {
//		unsigned int len = std::strlen(in);
//		for (unsigned int i = 0; i < len; i++) {
//			in[ i ] = (char)std::tolower( in[ i ] );
//		}
//	}
//
//	/**
//	Converts the supplied string to lower case.
//	@brief Converts to lower case.
//	@param in The string to convert.
//	@returns The input string converted into all lower case.
//	*/
//	std::string toLowerCase( const std::string in ) {
//		std::ostringstream oss("");
//		for (unsigned int i = 0; i < in.length(); i++) {
//			oss << (char)std::tolower( in[ i ] );
//		}
//		return oss.str();
//	}
//
//	/**
//	Converts the supplied string to upper case.
//	@brief Converts to upper case.
//	@param in The string to convert.
//	@returns The input string converted into all upper case.
//	*/
//	void toUpperCase( char *in ) {
//		unsigned int len = std::strlen(in);
//		for (unsigned int i = 0; i < len; i++) {
//			in[ i ] = (char)std::toupper( in[ i ] );
//		}
//	}
//
//	/**
//	Converts the supplied string to upper case.
//	@brief Converts to upper case.
//	@param in The string to convert.
//	@returns The input string converted into all upper case.
//	*/
//	std::string toUpperCase( const std::string in ) {
//		std::ostringstream oss("");
//		for (unsigned int i = 0; i < in.length(); i++) {
//			oss << (char)std::toupper( in[ i ] );
//		}
//		return oss.str();
//	}
}
