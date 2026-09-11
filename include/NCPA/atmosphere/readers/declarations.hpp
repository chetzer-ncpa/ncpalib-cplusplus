#pragma once

namespace NCPA {
    namespace atmos {

        enum class reader_1d_t { NCPAPROP };
        enum class reader_2d_t {
            NCPAPROP,
            NCPAPROP_STRATIFIED,
            NCPAPROP_PIECEWISE_STRATIFIED
        };
        enum class reader_3d_t { NCPAPROP, NCPAPROP_STRATIFIED };

        class AtmosphereReader1D;
        class AtmosphereReader2D;
        class AtmosphereReader3D;

        // readers
        class _abstract_atmosphere_reader_1d;
        class ncpaprop_atmosphere_reader_1d;
        class _abstract_atmosphere_reader_2d;
        class ncpaprop_atmosphere_reader_2d;
        class ncpaprop_atmosphere_reader_stratified_2d;
        class _abstract_atmosphere_reader_3d;
        class ncpaprop_atmosphere_reader_3d;
        class ncpaprop_atmosphere_reader_stratified_3d;
    }  // namespace atmos
}  // namespace NCPA
