#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/HTRZAlgorithm.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationDataFormats/interface/TTRoad.h"

#include "boost/multi_array.hpp"

#include <set>
// #include <unordered_map>
#include <unordered_set>
#include <array>
#include <iostream>
#include <cstdlib>

using namespace slhcl1tt;


HTRZAlgorithm::HTRZAlgorithm():
  mode_              (HTRZ_2D_COTANTHETA_Z0),
  stub_accept_policy_( LOOSE_ALL_NEIGHBOURS),
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


bool majority_all_layers(HTRZCell const& cell, unsigned int const threshold)
{
//   std::cerr << "call to majority_all_layers" << std::endl;
  
  unsigned int count = 0;
  
  for (unsigned int i = 0; i < HTRZAlgorithm::NLAYERS; ++i)
  {
//     std::cerr << "Layer " << i << " stubrefs vector size = " << cell.stubrefs_by_layer[i].size()  << std::endl;
    
    if (!cell.stubrefs_by_layer[i].empty())
    {
//       std::cerr << "count = " << count << std::endl;
      ++count;
    }
  }
  
//   for (auto const& layer : cell.stubrefs_by_layer)
//     if (!layer.empty())
//     {
//       std::cerr << "count = " << count << std::endl;
//       ++count;
//     }
  
  return count >= threshold;
}


bool majority_ps_layers(HTRZCell const& cell, unsigned int const threshold)
{
  unsigned int count = 0;
  
  for (unsigned int i = 0; i < 3; ++i)
    if (!cell.stubrefs_by_layer[i].empty())
      ++count;
  
  return count >= threshold;
}



constexpr static const double lut_pixel_halfwidth_z_bylayer[] = {0.0722812068078, 0.0722812068078, 0.0722812068078, 2.5125007629395, 2.5125007629395, 2.5125007629395};



TTRoad HTRZAlgorithm::Filter(slhcl1tt::TTRoad const& input_road, TTRoadReader const& reader)
{
//   std::cerr 
//     << "In road nstubs: " 
//       << input_road.stubRefs.at(0).size() << " "
//       << input_road.stubRefs.at(1).size() << " "
//       << input_road.stubRefs.at(2).size() << " "
//       << input_road.stubRefs.at(3).size() << " "
//       << input_road.stubRefs.at(4).size() << " "
//       << input_road.stubRefs.at(5).size() << " ";
  
  TTRoad output_road;
  
  std::unordered_set<unsigned int> passing_stubrefs;
  
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
//           ht_bounds_cotantheta.push_back( min_cotantheta_ * (double(nbins_cotantheta_ - iCot)/double(nbins_cotantheta_)) + max_cotantheta_ * (double(iCot)/double(nbins_cotantheta_)) );
          ht_bounds_cotantheta.push_back( min_cotantheta_ + (max_cotantheta_ - min_cotantheta_) * double(iCot)/double(nbins_cotantheta_) );
        }

        for (unsigned iZ0 = 0; iZ0 < nbins_z0_ + 1; ++iZ0)
        {
          ht_bounds_z0.push_back( min_z0_ * (double(nbins_z0_ - iZ0)/double(nbins_z0_)) + max_z0_ * (double(iZ0)/double(nbins_z0_)) );
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
            std::vector<                       bool > bound_accept (nbins_cotantheta_ + 1,  false);
            std::vector< std::vector< unsigned int> > bound_iZ0    (nbins_z0_         + 1,  std::vector<unsigned>());
            
            
            
            // Loop on the (nbins_cotantheta_ + 1) bin edges
            for (unsigned iCot = 0; iCot < nbins_cotantheta_ + 1; ++iCot)
            {
              double const z0_lo = stub_and_cotantheta_to_z0(stub_z - lut_pixel_halfwidth_z_bylayer[iLayer], stub_r, ht_bounds_cotantheta[iCot]);
              double const z0_hi = stub_and_cotantheta_to_z0(stub_z + lut_pixel_halfwidth_z_bylayer[iLayer], stub_r, ht_bounds_cotantheta[iCot]);
              
              std::cerr << "z0_lo = " << z0_lo << "  z0_hi = " << z0_hi << std::endl;
              assert(z0_lo < z0_hi);
              
//               std::cerr << "stub_z " << stub_z << " stub_r " << stub_r << " cotan theta bound " << iCot << " with value " << ht_bounds_cotantheta[iCot] << " ---> " << z0 << std::endl;
              
              if      (z0_hi < ht_bounds_z0.front())
              {
                // The z0 at this column falls out of the bottom of the matrix
                bound_accept[iCot] = false;
              }
              else if (z0_lo > ht_bounds_z0.back() )
              {
                // The z0 at this column falls out of the top of the matrix
                bound_accept[iCot] = false;
              }
              else
              {
                // We are inside the matrix in the vertical direction
                unsigned int const iZ0_lo = z0_lo <= ht_bounds_z0.front() ? 0 : 
                (
                  z0_lo >= ht_bounds_z0.back() ? nbins_z0_ - 1 : 
                    std::distance(ht_bounds_z0.begin(), std::lower_bound(ht_bounds_z0.begin(), ht_bounds_z0.end(), z0_lo) )
                );
                
                unsigned int const iZ0_hi = z0_hi <= ht_bounds_z0.front() ? 0 : 
                (
                  z0_hi >= ht_bounds_z0.back() ? nbins_z0_ - 1 : 
                    std::distance(ht_bounds_z0.begin(), std::lower_bound(ht_bounds_z0.begin(), ht_bounds_z0.end(), z0_hi) )
                );
                
                std::cerr << "z0_lo = " << z0_lo << "  z0_hi = " << z0_hi << "  iZ0_lo = " << iZ0_lo << "  iZ0_hi = " << iZ0_hi << std::endl;
                assert(iZ0_lo <= iZ0_hi);
                
                // This true statement can be changed to more complicated requirements.
                // For now the matrix boundaries (-15.0,15.0) are enough to guarantee 
                // the respect of the condition we care about: that the track candidate 
                // is compatible with a real track from collisions.
                bound_accept[iCot] = true;
                
                for (unsigned int iZ0 = iZ0_lo; iZ0 <= iZ0_hi; ++iZ0)
                  bound_iZ0[iCot].push_back(iZ0);
                
                assert(!bound_iZ0[iCot].empty());
              }
            }
            
            // Accept/reject the stubs based on the values at the boundaries and the 
            // set acceptance policy.
            for (unsigned iCot = 0; iCot < nbins_cotantheta_; ++iCot)
            {
//               std::cerr << "layer " << iLayer << " stubRef " << stubRef << " cotan theta bin " << iCot;
              
              if (!bound_accept[iCot] && !bound_accept[iCot + 1])
              {
                // Both borders are invalid, skip the stub in this column
//                 std::cerr << " both borders invalid --> reject " << std::endl;
                continue;
              }
              else if (bound_accept[iCot] && bound_accept[iCot + 1])
              {
//                 std::cerr << " both borders valid --> accept " << std::endl;
                // Both borders are valid, save as many cells as are crossed by the line
//                 if (bound_iZ0[iCot] <= bound_iZ0[iCot + 1])
//                   for (unsigned iZ0 = bound_iZ0[iCot]; iZ0 < bound_iZ0[iCot + 1] + 1; ++iZ0)
//                     ht_matrix[iCot][iZ0].stubrefs_by_layer[iLayer].push_back(stubRef);
//                 else
//                 {
//                   for (unsigned iZ0 = bound_iZ0[iCot + 1]; iZ0 < bound_iZ0[iCot] + 1; ++iZ0)
//                     ht_matrix[iCot][iZ0].stubrefs_by_layer[iLayer].push_back(stubRef);
//                 }
                
                std::set<unsigned> merged_iZ0( bound_iZ0[iCot    ].begin(), bound_iZ0[iCot    ].end() );
                merged_iZ0            .insert( bound_iZ0[iCot + 1].begin(), bound_iZ0[iCot + 1].end() );
                
                for (unsigned iZ0 : merged_iZ0)
                  ht_matrix[iCot][iZ0].stubrefs_by_layer[iLayer].push_back(stubRef);
              }
              else
              {
//                 std::cerr << " one border valid, another invalid --> loookup policy ";
                // One border is valid, the other is not. Follow the set acceptance policy.
                switch (stub_accept_policy_)
                {
                  case LOOSE_ALL_NEIGHBOURS:
                    {
//                       std::cerr << " --> all neighbors added " << std::endl;
                      
//                       if (bound_iZ0[iCot] <= bound_iZ0[iCot + 1])
//                         for (unsigned iZ0 = bound_iZ0[iCot]; iZ0 < bound_iZ0[iCot + 1] + 1; ++iZ0)
//                           ht_matrix[iCot][iZ0].stubrefs_by_layer[iLayer].push_back(stubRef);
//                       else
//                       {
//                         for (unsigned iZ0 = bound_iZ0[iCot + 1]; iZ0 < bound_iZ0[iCot] + 1; ++iZ0)
//                           ht_matrix[iCot][iZ0].stubrefs_by_layer[iLayer].push_back(stubRef);
//                       }
                      
                      std::set<unsigned> merged_iZ0( bound_iZ0[iCot    ].begin(), bound_iZ0[iCot    ].end() );
                      merged_iZ0            .insert( bound_iZ0[iCot + 1].begin(), bound_iZ0[iCot + 1].end() );
                      
                      std::cerr << "A     merged_iZ0.size() = " << merged_iZ0.size() << std::endl;
                      unsigned iZ0_end = *(merged_iZ0.rbegin());
                      std::cerr << "B     iZ0_end = " << iZ0_end << std::endl;
                      for (unsigned iZ0 = *(merged_iZ0.begin()); iZ0 < iZ0_end; ++iZ0)
                      {
                        std::cerr << "C     iZ0 = " << iZ0 << std::endl;
                        ht_matrix[iCot][iZ0].stubrefs_by_layer[iLayer].push_back(stubRef);
                      }
                      
                    }
                    break;
                  case MEDIUM_NEAR_NEIGHBOUR:
                    {
//                       std::cerr << " --> near neighbor added " << std::endl;
                      
//                       if (bound_accept[iCot])
//                         ht_matrix[iCot][ bound_iZ0[iCot]     ].stubrefs_by_layer[iLayer].push_back(stubRef);
//                       else
//                         ht_matrix[iCot][ bound_iZ0[iCot + 1] ].stubrefs_by_layer[iLayer].push_back(stubRef);
                      
                      
                      if (bound_accept[iCot])
                        for (unsigned iZ0 : bound_iZ0[iCot    ])
                          ht_matrix[iCot][iZ0].stubrefs_by_layer[iLayer].push_back(stubRef);
                      else
                        for (unsigned iZ0 : bound_iZ0[iCot + 1])
                          ht_matrix[iCot][iZ0].stubrefs_by_layer[iLayer].push_back(stubRef);
                    }
                    break;
                  case TIGHT_NO_NEIGHBOURS:
                    {
//                       std::cerr << " --> none added " << std::endl;
                    }
                    break;
                  default:
                    break;
                }
              }
            }
          }
        }
        
        // Loop over the matrix cells and take note of stubrefs in cells activated my the majority logic
        for (unsigned int iCot = 0; iCot < nbins_cotantheta_; ++iCot)
          for (unsigned int iZ0 = 0; iZ0 < nbins_z0_; ++iZ0)
            if ( majority_all_layers(ht_matrix[iCot][iZ0], 5) )
            {
              for (unsigned iLayer = 0; iLayer < NLAYERS; ++iLayer)
                for(unsigned int const stubref : ht_matrix[iCot][iZ0].stubrefs_by_layer[iLayer])
                {
                  passing_stubrefs.insert(stubref);
                }
            }
        
      }
      break;
    case HTRZ_1D_COTANTHETA:
      {
        
      }
      break;
    default:
      break;
  }
  
  output_road.patternRef    = input_road.patternRef    ;
  output_road.tower         = input_road.tower         ;
  output_road.patternInvPt  = input_road.patternInvPt  ;
  output_road.superstripIds = input_road.superstripIds ;
  
  for (unsigned iLayer = 0; iLayer < NLAYERS; ++iLayer)
  {
    output_road.stubRefs.push_back( std::vector<unsigned>() );
    
    for (unsigned int const in_stubref : input_road.stubRefs[iLayer])
    {
      if (passing_stubrefs.find(in_stubref) != passing_stubrefs.end())
        output_road.stubRefs[iLayer].push_back(in_stubref);
    }
  }
  
  
  output_road.nstubs = passing_stubrefs.size();
  
//   std::cerr 
//     << "    Out road nstubs: " 
//       << output_road.stubRefs.at(0).size() << " "
//       << output_road.stubRefs.at(1).size() << " "
//       << output_road.stubRefs.at(2).size() << " "
//       << output_road.stubRefs.at(3).size() << " "
//       << output_road.stubRefs.at(4).size() << " "
//       << output_road.stubRefs.at(5).size() << " ";
//   
//   std::cerr 
//       << std::endl;
  
  return output_road;
}

























