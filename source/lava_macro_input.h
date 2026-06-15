// 1: 10.1016/j.jcp.2024.112798
// 2: 10.1016/j.camwa.2025.02.014
// 3: https://arxiv.org/abs/2509.09460 Table 2
// 4: https://arxiv.org/abs/2509.09460 Table 3
#define SET_COEFFICIENTS 4

// FLUX_MODEL: selects the physical model
//   1: Depth-integrated lava model
//   2: Linear advection equation
#define FLUX_MODEL 1

// SET_TEST: selects the test case (first digit must match FLUX_MODEL)
//   FLUX_MODEL=1:
//     11: Well-balanced test
//     12: Travelling vortex
//     13: Lava pouring from a vent over flat topography, constant discharge
//     14: Lava pouring from a vent over flat topography, time-dependent discharge
//   FLUX_MODEL=2:
//     21: Travelling wave
#define SET_TEST 12
