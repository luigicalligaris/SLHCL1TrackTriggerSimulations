#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/HTRZAlgorithm.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationDataFormats/interface/TTRoad.h"

#include "boost/multi_array.hpp"

// #include <unordered_map>
#include <array>

using namespace slhcl1tt;


HTRZAlgorithm::HTRZAlgorithm():
  mode_            (HTRZ_1D_COTANTHETA),
  max_z0_          ( 15.0),
  min_z0_          (-15.0),
  max_cotantheta_  ( 13.5),
  min_cotantheta_  (-13.5),
  nbins_z0_        (    8),
  nbins_cotantheta_(    8)
{
  
}

HTRZAlgorithm::HTRZAlgorithm(HTRZAlgorithm const& rhs):
  mode_            (rhs.mode_            ),
  max_z0_          (rhs.max_z0_          ),
  min_z0_          (rhs.min_z0_          ),
  max_cotantheta_  (rhs.max_cotantheta_  ),
  min_cotantheta_  (rhs.min_cotantheta_  ),
  nbins_z0_        (rhs.nbins_z0_        ),
  nbins_cotantheta_(rhs.nbins_cotantheta_)
{
  
}

HTRZAlgorithm::HTRZAlgorithm(HTRZAlgorithm&&      rhs):
  mode_            (rhs.mode_            ),
  max_z0_          (rhs.max_z0_          ),
  min_z0_          (rhs.min_z0_          ),
  max_cotantheta_  (rhs.max_cotantheta_  ),
  min_cotantheta_  (rhs.min_cotantheta_  ),
  nbins_z0_        (rhs.nbins_z0_        ),
  nbins_cotantheta_(rhs.nbins_cotantheta_)
{
  
}

HTRZAlgorithm& HTRZAlgorithm::operator=(HTRZAlgorithm const& rhs)
{
  mode_             = rhs.mode_            ;
  max_z0_           = rhs.max_z0_          ;
  min_z0_           = rhs.min_z0_          ;
  max_cotantheta_   = rhs.max_cotantheta_  ;
  min_cotantheta_   = rhs.min_cotantheta_  ;
  nbins_z0_         = rhs.nbins_z0_        ;
  nbins_cotantheta_ = rhs.nbins_cotantheta_;

  return *this;
}

HTRZAlgorithm& HTRZAlgorithm::operator=(HTRZAlgorithm&&      rhs)
{
  mode_             = rhs.mode_            ;
  max_z0_           = rhs.max_z0_          ;
  min_z0_           = rhs.min_z0_          ;
  max_cotantheta_   = rhs.max_cotantheta_  ;
  min_cotantheta_   = rhs.min_cotantheta_  ;
  nbins_z0_         = rhs.nbins_z0_        ;
  nbins_cotantheta_ = rhs.nbins_cotantheta_;

  return *this;
}

HTRZAlgorithm::~HTRZAlgorithm()
{
  
}



struct HTRZCell
{
  std::array< std::vector< unsigned >, HTRZAlgorithm::NLAYERS > stubrefs_by_layer;
};


TTRoad HTRZAlgorithm::Filter(slhcl1tt::TTRoad const& input_road, TTRoadReader const& reader)
{
  boost::multi_array<HTRZCell, 2> ht_matrix(boost::extents[nbins_cotantheta_][nbins_z0_]);


  // Generate boundary vectors
  std::vector<double> ht_bounds_cotantheta;
  std::vector<double> ht_bounds_z0;

  for (unsigned iCot = 0; iCot < nbins_cotantheta_ + 1; ++iCot)
  {
    ht_bounds_cotantheta[iCot] = min_cotantheta_ * (double(nbins_cotantheta_ - iCot)/double(nbins_cotantheta_)) + max_cotantheta_ * (double(iCot)/double(nbins_cotantheta_));
  }

  for (unsigned iZ0 = 0; iZ0 < nbins_z0_ + 1; ++iZ0)
  {
    ht_bounds_cotantheta[iZ0] = min_cotantheta_ * (double(nbins_z0_ - iZ0)/double(nbins_z0_)) + max_cotantheta_ * (double(iZ0)/double(nbins_z0_));
  }




  // There is 1 superstrip per layer
  for (unsigned iLayer = 0; iLayer < NLAYERS; ++iLayer)
  {
    // Loop over stubs in superstrip
    for (unsigned stubRef : input_road.stubRefs[iLayer])
    {
      const float    stub_x          = reader.vb_x          ->at(stubRef);
      const float    stub_y          = reader.vb_y          ->at(stubRef);
      const float    stub_z          = reader.vb_z          ->at(stubRef);
      const float    stub_r          = reader.vb_r          ->at(stubRef);
      const float    stub_eta        = reader.vb_eta        ->at(stubRef);
      const float    stub_phi        = reader.vb_phi        ->at(stubRef);
      const float    stub_coordx     = reader.vb_coordx     ->at(stubRef);
      const float    stub_coordy     = reader.vb_coordy     ->at(stubRef);
      const float    stub_trigBend   = reader.vb_trigBend   ->at(stubRef);
      const float    stub_roughPt    = reader.vb_roughPt    ->at(stubRef);
      const float    stub_clusWidth0 = reader.vb_clusWidth0 ->at(stubRef);
      const float    stub_clusWidth1 = reader.vb_clusWidth1 ->at(stubRef);
      const unsigned stub_modId      = reader.vb_modId      ->at(stubRef);
      const int      stub_tpId       = reader.vb_tpId       ->at(stubRef);
      
      
      
      
    }
  }
  
  return input_road;
}

























