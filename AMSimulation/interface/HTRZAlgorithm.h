#ifndef AMSimulation_HTRZAlgorithm_h_
#define AMSimulation_HTRZAlgorithm_h_

#include "SLHCL1TrackTriggerSimulations/AMSimulationDataFormats/interface/TTRoad.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/TTRoadReader.h"

class HTRZAlgorithm
{
public:
  static constexpr const unsigned NLAYERS = 6;
  
  HTRZAlgorithm();
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
  enum mode_t {HTRZ_1D_COTANTHETA, HTRZ_2D_COTANTHETA_Z0};
  mode_t mode_;

  double   max_z0_          ;
  double   min_z0_          ;
  double   max_cotantheta_  ;
  double   min_cotantheta_  ;
  unsigned nbins_z0_        ;
  unsigned nbins_cotantheta_;
};

#endif