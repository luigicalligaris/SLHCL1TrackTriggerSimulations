#ifndef AMSimulation_HTRZAlgorithm_h_
#define AMSimulation_HTRZAlgorithm_h_

#include "SLHCL1TrackTriggerSimulations/AMSimulationDataFormats/interface/TTRoad.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/TTRoadReader.h"

// 2D HT matrix or 1D histogram?
enum HTRZAlgorithm_mode_t {HTRZ_1D_COTANTHETA, HTRZ_2D_COTANTHETA_Z0};
enum HTRZAlgorithm_stub_accept_policy_t {LOOSE_ALL_NEIGHBOURS, MEDIUM_NEAR_NEIGHBOUR, TIGHT_NO_NEIGHBOURS};

struct HTRZAlgorithmConfig
{
  HTRZAlgorithm_mode_t               mode              ;
  HTRZAlgorithm_stub_accept_policy_t stub_accept_policy;
  double                             max_z0            ;
  double                             min_z0            ;
  double                             max_cotantheta    ;
  double                             min_cotantheta    ;
  unsigned                           nbins_z0          ;
  unsigned                           nbins_cotantheta  ;
};


/***
** The class implementing the r-z Hough Transform Algorithm
**
**
** NOTES:
** Given two boundaries, one where the stub is acepted and the other 
** where it does not pass the selection criteria (whatever they are),
** with the stub line that crosses the first in a given vertical cell, 
** the second does it in another vertical cell 1 cell away. 
** The acceptance of the stub in the cells at the side of the 
** boundaries can follow these three strategies:
** 
** LOOSE_ALL_NEIGHBOURS
** | | | | | |                     | | | | | |
** | | | | | |                     | | | | | |
** |-|-|-|-|-|                     |-|-|-|-|-|
** N N | | | |                     N N*| | | |
** |-|-|-|-|-|  accepts stub = *   |-|-|-|-|-|
** | | Y Y | | ==================> | |*Y*Y*| |
** |-|-|-|-|-|                     |-|-|-|-|-|
** | | | | Y Y                     | | | |*Y*Y
** |-|-|-|-|-|                     |-|-|-|-|-|
** | | | | | |                     | | | | | |
** 
** MEDIUM_NEAR_NEIGHBOUR
** | | | | | |                     | | | | | |
** | | | | | |                     | | | | | |
** |-|-|-|-|-|                     |-|-|-|-|-|
** N N | | | |                     N N | | | |
** |-|-|-|-|-|  accepts stub = *   |-|-|-|-|-|
** | | Y Y | | ==================> | |*Y*Y*| |
** |-|-|-|-|-|                     |-|-|-|-|-|
** | | | | Y Y                     | | | |*Y*Y
** |-|-|-|-|-|                     |-|-|-|-|-|
** | | | | | |                     | | | | | |
** 
** TIGHT_NO_NEIGHBOURS
** | | | | | |                     | | | | | |
** | | | | | |                     | | | | | |
** |-|-|-|-|-|                     |-|-|-|-|-|
** N N | | | |                     N N | | | |
** |-|-|-|-|-|  accepts stub = *   |-|-|-|-|-|
** | | Y Y | | ==================> | | Y*Y*| |
** |-|-|-|-|-|                     |-|-|-|-|-|
** | | | | Y Y                     | | | |*Y*Y
** |-|-|-|-|-|                     |-|-|-|-|-|
** | | | | | |                     | | | | | |
**
**/
class HTRZAlgorithm
{
public:
  static constexpr const unsigned NLAYERS = 6;
  
  HTRZAlgorithm();
  HTRZAlgorithm(HTRZAlgorithmConfig const& config);
  HTRZAlgorithm(HTRZAlgorithm const& rhs);
  HTRZAlgorithm(HTRZAlgorithm&&      rhs);
  HTRZAlgorithm& operator=(HTRZAlgorithm const& rhs);
  HTRZAlgorithm& operator=(HTRZAlgorithm&&      rhs);
  ~HTRZAlgorithm();

  slhcl1tt::TTRoad Filter(slhcl1tt::TTRoad const& input_road, slhcl1tt::TTRoadReader const& reader);

private:
  inline double stub_and_cotantheta_to_z0(double const zstub, double const rstub, double const cotantheta) const
  {
    return zstub - rstub * cotantheta;
  }

  inline double stub_and_z0_to_cotantheta(double const zstub, double const rstub, double const z0) const
  {
    return (zstub - z0) / rstub;
  }

private:
  HTRZAlgorithm_mode_t               mode_              ;
  HTRZAlgorithm_stub_accept_policy_t stub_accept_policy_;

  double   max_z0_          ;
  double   min_z0_          ;
  double   max_cotantheta_  ;
  double   min_cotantheta_  ;
  unsigned nbins_z0_        ;
  unsigned nbins_cotantheta_;
};

#endif