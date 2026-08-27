#pragma once


#include "NCPA/interpolation/abstract_spline_1d.hpp"
#include "NCPA/interpolation/abstract_spline_2d.hpp"
#include "NCPA/interpolation/abstract_spline_3d.hpp"
#include "NCPA/interpolation/builders.hpp"
#include "NCPA/interpolation/defines.hpp"
#include "NCPA/interpolation/Extrapolator1D.hpp"
#include "NCPA/interpolation/gsl.hpp"
#include "NCPA/interpolation/Interpolator1D.hpp"
#include "NCPA/interpolation/Interpolator2D.hpp"
#include "NCPA/interpolation/Interpolator3D.hpp"
#include "NCPA/interpolation/lanl.hpp"
#include "NCPA/interpolation/nearest_neighbor_spline_1d.hpp"
#include "NCPA/interpolation/nearest_neighbor_spline_2d.hpp"
#include "NCPA/interpolation/stratified_spline_2d.hpp"
#include "NCPA/interpolation/types.hpp"

#ifndef NCPA_INTERPOLATION_DEFAULT_1D
#  ifdef NCPA_INTERPOLATION_GSL_INTERPOLATION_AVAILABLE
#    ifdef NCPA_INTERPOLATION_GSL_STEFFEN_SPLINE_AVAILABLE
#      define NCPA_INTERPOLATION_DEFAULT_1D \
          NCPA::interpolation::interpolator_1d_type_t::GSL_STEFFEN
#    else
#      define NCPA_INTERPOLATION_DEFAULT_1D \
          NCPA::interpolation::interpolator_1d_type_t::GSL_CUBIC
#    endif
#  else
#    define NCPA_INTERPOLATION_DEFAULT_1D \
        NCPA::interpolation::interpolator_1d_type_t::LANL_CUBIC
#  endif
#endif

#ifndef NCPA_INTERPOLATION_DEFAULT_2D
#  define NCPA_INTERPOLATION_DEFAULT_2D \
      NCPA::interpolation::interpolator_2d_type_t::LANL_BICUBIC
#endif

#ifndef NCPA_INTERPOLATION_DEFAULT_3D
#  define NCPA_INTERPOLATION_DEFAULT_3D \
      NCPA::interpolation::interpolator_3d_type_t::LANL_HYBRID
#endif
