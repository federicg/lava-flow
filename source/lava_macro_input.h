
//
#define SET_COEFFICIENTS 4

// 1: Depth-Integrated lava model,
// 2: Linear advection equation
#define FLUX_MODEL 1

// 11: WB test for FLUX_MODEL=1
// 12: Travelling vortex for FLUX_MODEL=1
// 13: Pouring of lava from a vent over a flat topography, with constant vent discharge
// 14: Pouring of lava from a vent over a flat topography, with time-dependent vent discharge
// 21: Traveling wave available just in case FLUX_MODEL=2
#define SET_TEST 12 
