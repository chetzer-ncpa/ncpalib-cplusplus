#pragma once

#include <cmath>
#include <stdexcept>

namespace NCPA::geographic {

	// template<typename T = double>
	// class GeographicLimits {
	// 	public:
	// 		void check_latitude( T lat ) {

	// 		}
	// 	private:
	// 		T _minlat, _minlon, _maxlat, _maxlon;
	// };

	/**
	@brief Normalizes an azimuth into the range [0,360)
	@param in The azimuth value to normalize, in degrees.
	@returns The normalized azimuth, in degrees.
	*/
	template<typename T>
	T normalize_azimuth( T in, bool radians = false ) {
		double interval = radians ? 2.0 * M_PI : 360.0;
		double out = (double)in;
		while (out < 0.0) {
			out += interval;
		}
		while (out >= interval) {
			out -= interval;
		}
		return (T)out;
	}

}
