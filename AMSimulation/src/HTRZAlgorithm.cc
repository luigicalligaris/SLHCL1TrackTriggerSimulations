#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/HTRZAlgorithm.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationDataFormats/interface/TTRoad.h"

#include "boost/multi_array.hpp"

// #include <unordered_map>
#include <array>

using namespace slhcl1tt;


HTRZAlgorithm::HTRZAlgorithm():
  mode_              (  HTRZ_1D_COTANTHETA),
  stub_accept_policy_(LOOSE_ALL_NEIGHBOURS),
  max_z0_            ( 15.0),
  min_z0_            (-15.0),
  max_cotantheta_    ( 13.5),
  min_cotantheta_    (-13.5),
  nbins_z0_          (    8),
  nbins_cotantheta_  (    8)
{
  assert(max_z0_         >= min_z0_        );
  assert(max_cotantheta_ >= min_cotantheta_);
}

HTRZAlgorithm::HTRZAlgorithm(HTRZAlgorithmConfig const& config):
  mode_              (config.mode              ),
  stub_accept_policy_(config.stub_accept_policy),
  max_z0_            (config.max_z0            ),
  min_z0_            (config.min_z0            ),
  max_cotantheta_    (config.max_cotantheta    ),
  min_cotantheta_    (config.min_cotantheta    ),
  nbins_z0_          (config.nbins_z0          ),
  nbins_cotantheta_  (config.nbins_cotantheta  )
{
  assert(max_z0_         >= min_z0_        );
  assert(max_cotantheta_ >= min_cotantheta_);
}

HTRZAlgorithm::HTRZAlgorithm(HTRZAlgorithm const& rhs):
  mode_              (rhs.mode_              ),
  stub_accept_policy_(rhs.stub_accept_policy_),
  max_z0_            (rhs.max_z0_            ),
  min_z0_            (rhs.min_z0_            ),
  max_cotantheta_    (rhs.max_cotantheta_    ),
  min_cotantheta_    (rhs.min_cotantheta_    ),
  nbins_z0_          (rhs.nbins_z0_          ),
  nbins_cotantheta_  (rhs.nbins_cotantheta_  )
{
  assert(max_z0_         >= min_z0_        );
  assert(max_cotantheta_ >= min_cotantheta_);
}

HTRZAlgorithm::HTRZAlgorithm(HTRZAlgorithm&&      rhs):
  mode_              (rhs.mode_              ),
  stub_accept_policy_(rhs.stub_accept_policy_),
  max_z0_            (rhs.max_z0_            ),
  min_z0_            (rhs.min_z0_            ),
  max_cotantheta_    (rhs.max_cotantheta_    ),
  min_cotantheta_    (rhs.min_cotantheta_    ),
  nbins_z0_          (rhs.nbins_z0_          ),
  nbins_cotantheta_  (rhs.nbins_cotantheta_  )
{
  assert(max_z0_         >= min_z0_        );
  assert(max_cotantheta_ >= min_cotantheta_);
}

HTRZAlgorithm& HTRZAlgorithm::operator=(HTRZAlgorithm const& rhs)
{
  mode_                = rhs.mode_               ;
  stub_accept_policy_  = rhs.stub_accept_policy_ ;
  max_z0_              = rhs.max_z0_             ;
  min_z0_              = rhs.min_z0_             ;
  max_cotantheta_      = rhs.max_cotantheta_     ;
  min_cotantheta_      = rhs.min_cotantheta_     ;
  nbins_z0_            = rhs.nbins_z0_           ;
  nbins_cotantheta_    = rhs.nbins_cotantheta_   ;

  assert(max_z0_         >= min_z0_        );
  assert(max_cotantheta_ >= min_cotantheta_);

  return *this;
}

HTRZAlgorithm& HTRZAlgorithm::operator=(HTRZAlgorithm&&      rhs)
{
  mode_                = rhs.mode_               ;
  stub_accept_policy_  = rhs.stub_accept_policy_ ;
  max_z0_              = rhs.max_z0_             ;
  min_z0_              = rhs.min_z0_             ;
  max_cotantheta_      = rhs.max_cotantheta_     ;
  min_cotantheta_      = rhs.min_cotantheta_     ;
  nbins_z0_            = rhs.nbins_z0_           ;
  nbins_cotantheta_    = rhs.nbins_cotantheta_   ;

  assert(max_z0_         >= min_z0_        );
  assert(max_cotantheta_ >= min_cotantheta_);

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
  TTRoad output_road;
  
  switch (mode_)
  {
    case HTRZ_2D_COTANTHETA_Z0:
      {
        // Instantiate 2D HT matrix (x = cotan theta, y = z0)
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
          ht_bounds_z0[iZ0] = min_z0_ * (double(nbins_z0_ - iZ0)/double(nbins_z0_)) + max_z0_ * (double(iZ0)/double(nbins_z0_));
        }

        // The way we build the boundaries vectors makes them sorted in an increasing value order, no duplicates


        // There is 1 superstrip per layer
        for (unsigned iLayer = 0; iLayer < NLAYERS; ++iLayer)
        {
          // Loop over stubs in superstrip
          for (unsigned stubRef : input_road.stubRefs[iLayer])
          {
//             const float    stub_x          = reader.vb_x          ->at(stubRef);
//             const float    stub_y          = reader.vb_y          ->at(stubRef);
            const float    stub_z          = reader.vb_z          ->at(stubRef);
            const float    stub_r          = reader.vb_r          ->at(stubRef);
//             const float    stub_eta        = reader.vb_eta        ->at(stubRef);
//             const float    stub_phi        = reader.vb_phi        ->at(stubRef);
//             const float    stub_coordx     = reader.vb_coordx     ->at(stubRef);
//             const float    stub_coordy     = reader.vb_coordy     ->at(stubRef);
//             const float    stub_trigBend   = reader.vb_trigBend   ->at(stubRef);
//             const float    stub_roughPt    = reader.vb_roughPt    ->at(stubRef);
//             const float    stub_clusWidth0 = reader.vb_clusWidth0 ->at(stubRef);
//             const float    stub_clusWidth1 = reader.vb_clusWidth1 ->at(stubRef);
//             const unsigned stub_modId      = reader.vb_modId      ->at(stubRef);
//             const int      stub_tpId       = reader.vb_tpId       ->at(stubRef);
            
            // Instantiate & bootstrap the boundary vectors
            std::vector<        bool, nbins_cotantheta_ + 1> bound_accept = {};
            std::vector<unsigned int, nbins_cotantheta_ + 1> bound_iZ0    = {};
            
            // Loop on the (nbins_cotantheta_ + 1) bin edges
            for (unsigned iCot = 0; iCot < nbins_cotantheta_ + 1; ++iCot)
            {
              double const z0 = stub_and_cotantheta_to_z0(stub_z, stub_r, ht_bounds_cotantheta[iCot]);
              
              if      (z0 < ht_bounds_z0.front())
              {
                // The z0 at this column falls out of the bottom of the matrix
                bound_accept[iCot] = false;
                continue;
              }
              else if (z0 > ht_bounds_z0.back() )
              {
                // The z0 at this column falls out of the top of the matrix
                bound_accept[iCot] = false;
                continue;
              }
              else
              {
                if      (z0 == ht_bounds_z0.front())
                {
                  // z0 compatible with the bottom edge (i.e. we are at the very rare crossing point between 2 cells at the matrix bottom)
                  // reminder: we keep semi-open intervals [lower, upper)
                  bound_accept[iCot] = true;
                  bound_iZ0   [iCot] = 0;
                }
                else if (z0 == ht_bounds_z0.back() )
                {
                  // z0 compatible with the top edge (i.e. we are at the very rare crossing point between 2 cells at the matrix top)
                  // reminder: we keep semi-open intervals [lower, upper)
                  bound_accept[iCot] = false;
                  bound_iZ0   [iCot] = nbins_z0_ - 1;
                }
                else
                {
                  // We are strictly inside the matrix in the vertical direction
                  unsigned int const iZ0 = std::distance(ht_bounds_z0.begin(), std::lower_bound(ht_bounds_z0.begin(), ht_bounds_z0.end(), z0) );
                  
                  // This true statement can be changed to more complicated requirements.
                  // For now the matrix boundaries (-15.0,15.0) are enough to guarantee 
                  // the respect of the condition we care about: that the track candidate 
                  // is compatible with a real track from collisions.
                  bound_accept[iCot] = true;
                  bound_iZ0   [iCot] = iZ0;
                }
              }
            }
            
            // Accept/reject the stubs based on the values at the boundaries and the 
            // set acceptance policy.
            for (unsigned iCot = 0; iCot < nbins_cotantheta_; ++iCot)
            {
            
            }
            ht_matrix[iCot][            0].stubrefs_by_layer[iLayer].push_back(stubRef);
            ht_matrix[iCot][nbins_z0_ - 1].stubrefs_by_layer[iLayer].push_back(stubRef);
            ht_matrix[iCot][          iZ0].stubrefs_by_layer[iLayer].push_back(stubRef);
            
          }
        }
        
        return input_road;
      }
      break;
    case HTRZ_1D_COTANTHETA:
      {
        
      }
      break;
    default:
      break;
  }
  
  return output_road;
}

























