#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/HTRZAlgorithm.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationDataFormats/interface/TTRoad.h"

#include "boost/multi_array.hpp"

#include <set>
// #include <unordered_map>
#include <unordered_set>
#include <array>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <iterator>


#define HTRZALGORITHM_DEFAULT_MODE                        HTRZ_2D_COTANTHETA_Z0
#define HTRZALGORITHM_DEFAULT_STUB_ACCEPT_POLICY          LOOSE_ALL_NEIGHBOURS
#define HTRZALGORITHM_DEFAULT_COLOR_OUTPUT                true
#define HTRZALGORITHM_DEFAULT_VERBOSE                     -1
#define HTRZALGORITHM_DEFAULT_NBINS_Z0                    8
#define HTRZALGORITHM_DEFAULT_NBINS_COTANTHETA            8
#define HTRZALGORITHM_DEFAULT_NBINS_THRESHOLD_ALL_LAYERS  4
#define HTRZALGORITHM_DEFAULT_NBINS_THRESHOLD_PS_LAYERS   2
#define HTRZALGORITHM_DEFAULT_MAX_Z0                       15.0
#define HTRZALGORITHM_DEFAULT_MIN_Z0                      -15.0
#define HTRZALGORITHM_DEFAULT_MAX_COTANTHETA               13.5
#define HTRZALGORITHM_DEFAULT_MIN_COTANTHETA              -13.5


using namespace slhcl1tt;

HTRZAlgorithm::HTRZAlgorithm():
  config_{
    HTRZALGORITHM_DEFAULT_MODE                      ,
    HTRZALGORITHM_DEFAULT_STUB_ACCEPT_POLICY        ,
    HTRZALGORITHM_DEFAULT_COLOR_OUTPUT              ,
    HTRZALGORITHM_DEFAULT_VERBOSE                   ,
    HTRZALGORITHM_DEFAULT_NBINS_Z0                  ,
    HTRZALGORITHM_DEFAULT_NBINS_COTANTHETA          ,
    HTRZALGORITHM_DEFAULT_NBINS_THRESHOLD_ALL_LAYERS,
    HTRZALGORITHM_DEFAULT_NBINS_THRESHOLD_PS_LAYERS ,
    HTRZALGORITHM_DEFAULT_MAX_Z0                    ,
    HTRZALGORITHM_DEFAULT_MIN_Z0                    ,
    HTRZALGORITHM_DEFAULT_MAX_COTANTHETA            ,
    HTRZALGORITHM_DEFAULT_MIN_COTANTHETA            }
{
  assert(config_.max_z0         >= config_.min_z0        );
  assert(config_.max_cotantheta >= config_.min_cotantheta);

  RegenerateBoundaries();
}


HTRZAlgorithm::HTRZAlgorithm(HTRZAlgorithm const& rhs):
  config_(rhs.config_)
{
  assert(config_.max_z0         >= config_.min_z0        );
  assert(config_.max_cotantheta >= config_.min_cotantheta);

  ht_bounds_cotantheta_ = rhs.ht_bounds_cotantheta_;
  ht_bounds_z0_         = rhs.ht_bounds_z0_        ;
}

HTRZAlgorithm::HTRZAlgorithm(HTRZAlgorithm&&      rhs):
  config_(rhs.config_)
{
  assert(config_.max_z0         >= config_.min_z0        );
  assert(config_.max_cotantheta >= config_.min_cotantheta);

  ht_bounds_cotantheta_ = std::move(rhs.ht_bounds_cotantheta_);
  ht_bounds_z0_         = std::move(rhs.ht_bounds_z0_        );
}

HTRZAlgorithm& HTRZAlgorithm::operator=(HTRZAlgorithm const& rhs)
{
  config_ = rhs.config_;

  assert(config_.max_z0         >= config_.min_z0        );
  assert(config_.max_cotantheta >= config_.min_cotantheta);

  ht_bounds_cotantheta_ = rhs.ht_bounds_cotantheta_;
  ht_bounds_z0_         = rhs.ht_bounds_z0_        ;

  return *this;
}

HTRZAlgorithm& HTRZAlgorithm::operator=(HTRZAlgorithm&&      rhs)
{
  config_ = rhs.config_;

  assert(config_.max_z0         >= config_.min_z0        );
  assert(config_.max_cotantheta >= config_.min_cotantheta);

  ht_bounds_cotantheta_ = std::move(rhs.ht_bounds_cotantheta_);
  ht_bounds_z0_         = std::move(rhs.ht_bounds_z0_        );

  return *this;
}

HTRZAlgorithm::~HTRZAlgorithm()
{
  
}


HTRZAlgorithm::HTRZAlgorithm(HTRZAlgorithmConfig const& config):
  config_(config)
{
  assert(config_.max_z0         >= config_.min_z0        );
  assert(config_.max_cotantheta >= config_.min_cotantheta);

  RegenerateBoundaries();
}

void HTRZAlgorithm::LoadConfig(HTRZAlgorithmConfig const& config)
{
  config_ = config;

  assert(config_.max_z0         >= config_.min_z0        );
  assert(config_.max_cotantheta >= config_.min_cotantheta);

  RegenerateBoundaries();
}



void HTRZAlgorithm::RegenerateBoundariesCotanTheta()
{
  // Generate boundary vectors

  for (unsigned iCot = 0; iCot < config_.nbins_cotantheta + 1; ++iCot)
  {
    ht_bounds_cotantheta_.push_back( config_.min_cotantheta + (config_.max_cotantheta - config_.min_cotantheta) * double(iCot)/double(config_.nbins_cotantheta) );
  }

  
  if (config_.verbose >= 4)
  {
    std::cout << "New bounds in cotan theta = [";
    for (auto val : ht_bounds_cotantheta_)
    {
      std::cout << val << ",";
    }
    std::cout << "]" << std::endl;
  }
}

void HTRZAlgorithm::RegenerateBoundariesZ0()
{
  for (unsigned iZ0 = 0; iZ0 < config_.nbins_z0 + 1; ++iZ0)
  {
    ht_bounds_z0_.push_back( config_.min_z0 * (double(config_.nbins_z0 - iZ0)/double(config_.nbins_z0)) + config_.max_z0 * (double(iZ0)/double(config_.nbins_z0)) );
  }
  
  if (config_.verbose >= 4)
  {
    std::cout << "New bounds in z0 = [";
    for (auto val : ht_bounds_z0_)
    {
      std::cout << val << ",";
    }
    std::cout << "]" << std::endl;
  }
}



struct HTRZCell
{
  std::array< std::vector< unsigned >, HTRZAlgorithm::NLAYERS > stubrefs_by_layer;
};




unsigned int count_stubs_all_layers(HTRZCell const& cell)
{
  unsigned int count = 0;
  
  for (unsigned int i = 0; i < HTRZAlgorithm::NLAYERS; ++i)
  {
    if (!cell.stubrefs_by_layer[i].empty())
    {
      ++count;
    }
  }
  
  return count;
}

unsigned int count_stubs_ps_layers(HTRZCell const& cell)
{
  unsigned int count = 0;
  
  for (unsigned int i = 0; i < 3; ++i)
  {
    if (!cell.stubrefs_by_layer[i].empty())
    {
      ++count;
    }
  }
  
  return count;
}


bool majority_all_layers(HTRZCell const& cell, unsigned int const threshold)
{
  return count_stubs_all_layers(cell) >= threshold;
}

bool majority_ps_layers(HTRZCell const& cell, unsigned int const threshold)
{
  return count_stubs_ps_layers(cell) >= threshold;
}



constexpr static const double lut_pixel_halfwidth_z_bylayer[] = {0.0722812068078, 0.0722812068078, 0.0722812068078, 2.5125007629395, 2.5125007629395, 2.5125007629395};


TTRoad HTRZAlgorithm::Filter(slhcl1tt::TTRoad const& input_road, TTRoadReader const& reader)
{
  TTRoad output_road;
  
  std::unordered_set<unsigned int> passing_stubrefs;
  
  switch (config_.mode)
  {
    case NULL_ALGO:
      output_road = input_road;
      break;
    
    case HTRZ_2D_COTANTHETA_Z0:
      {
        // Instantiate 2D HT matrix (x = cotan theta, y = z0)
        boost::multi_array<HTRZCell, 2> ht_matrix(boost::extents[config_.nbins_cotantheta][config_.nbins_z0]);

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
            std::vector<                       bool > bound_accept (config_.nbins_cotantheta + 1,  false);
            std::vector< std::vector< unsigned int> > bound_iZ0    (config_.nbins_cotantheta + 1,  std::vector<unsigned>());
            
            
            // Loop on the (config_.nbins_cotantheta + 1) bin edges
            for (unsigned iCot = 0; iCot < config_.nbins_cotantheta + 1; ++iCot)
            {
              double const z0_lo = stub_and_cotantheta_to_z0(stub_z - lut_pixel_halfwidth_z_bylayer[iLayer], stub_r, ht_bounds_cotantheta_[iCot]);
              double const z0_hi = stub_and_cotantheta_to_z0(stub_z + lut_pixel_halfwidth_z_bylayer[iLayer], stub_r, ht_bounds_cotantheta_[iCot]);
              
              if (config_.verbose >= 5)
              {
                std::cerr << "z0_lo = " << std::fixed << std::setprecision(4) << std::showpoint << std::showpos << std::setw(10) << z0_lo << "  z0_hi = " << std::fixed << std::setprecision(4) << std::showpoint << std::showpos << std::setw(10) << z0_hi << std::endl;
              }
              
              assert(z0_lo < z0_hi);
              
              if      (z0_hi <= ht_bounds_z0_.front())
              {
                // The z0 at this column falls out of the bottom of the matrix
                bound_accept[iCot] = false;
                
                if (config_.verbose >= 4)
                  std::cerr << "Boundary REJECTED z0_lo = " << z0_lo << "  z0_hi = " << z0_hi << std::endl;
              }
              else if (z0_lo >= ht_bounds_z0_.back() )
              {
                // The z0 at this column falls out of the top of the matrix
                bound_accept[iCot] = false;
                
                if (config_.verbose >= 4)
                  std::cerr << "Boundary REJECTED z0_lo = " << z0_lo << "  z0_hi = " << z0_hi << std::endl;
              }
              else
              {
                // We are inside the matrix in the vertical direction
                
                unsigned int const iZ0_lo = z0_lo <= ht_bounds_z0_.front() ? 0             : std::distance(ht_bounds_z0_.begin(), std::lower_bound(ht_bounds_z0_.begin(), ht_bounds_z0_.end(), z0_lo) ) - 1;
                
                unsigned int const iZ0_hi = z0_hi >= ht_bounds_z0_.back()  ? config_.nbins_z0 - 1 : std::distance(ht_bounds_z0_.begin(), std::lower_bound(ht_bounds_z0_.begin(), ht_bounds_z0_.end(), z0_hi) ) - 1;
                
                if (config_.verbose >= 4)
                  std::cerr << "Boundary ACCEPTED z0_lo = " << z0_lo << "  z0_hi = " << z0_hi << "  iZ0_lo = " << iZ0_lo << "  iZ0_hi = " << iZ0_hi << std::endl;
                
                assert(iZ0_lo <= iZ0_hi);
                
                
                // This true statement can be changed to more complicated requirements.
                // For now the matrix boundaries (-15.0,15.0) are enough to guarantee 
                // the respect of the condition we care about: that the track candidate 
                // is compatible with a real track from collisions.
                
                bound_accept[iCot] = true;
                
                for (unsigned int iZ0 = iZ0_lo; iZ0 <= iZ0_hi; ++iZ0)
                {
                  bound_iZ0[iCot].push_back(iZ0);
                }
                
                assert(!bound_iZ0[iCot].empty());
              }
            }
            
            
            // Accept/reject the stubs based on the values at the boundaries and the 
            // set acceptance policy.
            for (unsigned iCot = 0; iCot < config_.nbins_cotantheta; ++iCot)
            {
              if (config_.verbose >= 4)
                std::cerr << "layer " << iLayer << " stubRef " << stubRef << " cotan theta bin " << iCot;
              
              if (!bound_accept[iCot] && !bound_accept[iCot + 1])
              {
                // Both borders are invalid, skip the stub in this column
                if (config_.verbose >= 4)
                  std::cerr << " both borders invalid --> reject " << std::endl;
                
                continue;
              }
              else if (bound_accept[iCot] && bound_accept[iCot + 1])
              {
                if (config_.verbose >= 4)
                  std::cerr << " both borders valid --> accept " << std::endl;
                
                std::set<unsigned> merged_iZ0( bound_iZ0[iCot    ].begin(), bound_iZ0[iCot    ].end() );
                merged_iZ0            .insert( bound_iZ0[iCot + 1].begin(), bound_iZ0[iCot + 1].end() );
                
                for (unsigned iZ0 : merged_iZ0)
                  ht_matrix[iCot][iZ0].stubrefs_by_layer[iLayer].push_back(stubRef);
              }
              else
              {
                if (config_.verbose >= 4)
                  std::cerr << " one border valid, another invalid --> loookup policy ";
                // One border is valid, the other is not. Follow the set acceptance policy.
                switch (config_.stub_accept_policy)
                {
                  case LOOSE_ALL_NEIGHBOURS:
                    {
                      if (config_.verbose >= 4)
                        std::cerr << " --> all neighbors added " << std::endl;
                      
                      std::set<unsigned> merged_iZ0( bound_iZ0[iCot    ].begin(), bound_iZ0[iCot    ].end() );
                      merged_iZ0            .insert( bound_iZ0[iCot + 1].begin(), bound_iZ0[iCot + 1].end() );
                      
                      unsigned iZ0_end = *(merged_iZ0.rbegin());
                      
                      for (unsigned iZ0 = *(merged_iZ0.begin()); iZ0 < iZ0_end; ++iZ0)
                      {
//                         std::cerr << "C     iZ0 = " << iZ0 << std::endl;
                        ht_matrix[iCot][iZ0].stubrefs_by_layer[iLayer].push_back(stubRef);
                      }
                    }
                    break;
                  case MEDIUM_NEAR_NEIGHBOUR:
                    {
                      if (config_.verbose >= 4)
                        std::cerr << " --> near neighbor added " << std::endl;
                      
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
                      if (config_.verbose >= 4)
                        std::cerr << " --> none added " << std::endl;
                    }
                    break;
                  default:
                    break;
                }
              }
            }
          }
        }
        
        
        // Diagnostic print of the HT matrix
        
        if (config_.verbose >= 3)
        {
          std::cerr << std::endl << "HT Matrix of road from pattern "<< input_road.patternRef << ":" << std::endl;
          
          for (unsigned int iZ0 = 0; iZ0 < config_.nbins_z0; ++iZ0)
          {
            for (unsigned int iCot = 0; iCot < config_.nbins_cotantheta; ++iCot)
            {
              std::string color_str_begin = "";
              std::string color_str_end   = "";
              
              unsigned const count_all_lay = count_stubs_all_layers( ht_matrix[iCot][iZ0] );
              unsigned const count_ps_lay  = count_stubs_ps_layers ( ht_matrix[iCot][iZ0] );
              
              bool const maj_all_lay = count_all_lay >= config_.threshold_all_layers;
              bool const maj_ps_lay  = count_ps_lay  >= config_.threshold_ps_layers ;
              
              
              if (config_.color_output)
              {
                if (count_all_lay == 0)
                {
                  color_str_begin = "\033[2;37m"; // FAINT GREY
                  color_str_end   = "\033[0m";
                }
                else
                {
                  switch (maj_all_lay)
                  {
                    case true:
                    {
                      switch (maj_ps_lay)
                      {
                        case 1:
                          color_str_begin = "\033[1;31m"; // BOLD RED
                          color_str_end   = "\033[0m";
                          break;
                        case 0:
                          color_str_begin = "\033[1;32m"; // BOLD GREEN
                          color_str_end   = "\033[0m";
                          break;
                        default:
                          break;
                      }
                    }
                      break;
                    
                    case false:
                    {
                      switch (maj_ps_lay)
                      {
                        case 1:
                          color_str_begin = "\033[1;33m"; //BOLD YELLOW
                          color_str_end   = "\033[0m";
                          
                          break;
                        case 0:
                          color_str_begin = ""; //NO COLOR
                          color_str_end   = "";
                          break;
                        default:
                          break;
                      }
                    }
                      break;
                    
                    default:
                      break;
                  }
                }
              }
              
              std::cerr << color_str_begin << count_all_lay << color_str_end;
            }
            
            std::cerr << std::endl;
          }
        }
        
        // Loop over the matrix cells and take note of stubrefs in cells activated by the majority logic
        for (unsigned int iCot = 0; iCot < config_.nbins_cotantheta; ++iCot)
          for (unsigned int iZ0 = 0; iZ0 < config_.nbins_z0; ++iZ0)
            if ( majority_all_layers( ht_matrix[iCot][iZ0], config_.threshold_all_layers) && majority_ps_layers(ht_matrix[iCot][iZ0], config_.threshold_ps_layers) )
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
      {
        output_road.stubRefs[iLayer].push_back(in_stubref);
      }
    }
  }
  
  
  output_road.nstubs = passing_stubrefs.size();
  
  
  if (config_.verbose >= 4)
  {
    std::cout << "Number of stubs by layer for this road: ";
    for (auto& stubRefVec : output_road.stubRefs)
    {
      std::cout << stubRefVec.size() << " ";
    }
    std::cout << std::endl;
    
    if (config_.verbose >= 4)
    {
      std::cout << "Passing stubrefs in: ";
      for (unsigned iLayer = 0; iLayer < NLAYERS; ++iLayer)
      {
        std::cout << "layer " << iLayer << " ( ";
        for (unsigned int const in_stubref : input_road.stubRefs[iLayer])
        {
          if (passing_stubrefs.find(in_stubref) != passing_stubrefs.end())
          {
            std::cout << in_stubref << " ";
          }
        }
        std::cout << ")  ";
      }
      std::cout << std::endl;
      
      
      std::cout << "Failing stubrefs in: ";
      for (unsigned iLayer = 0; iLayer < NLAYERS; ++iLayer)
      {
        output_road.stubRefs.push_back( std::vector<unsigned>() );
        std::cout << "layer " << iLayer << " ( ";
        for (unsigned int const in_stubref : input_road.stubRefs[iLayer])
        {
          if (passing_stubrefs.find(in_stubref) == passing_stubrefs.end())
          {
            std::cout << in_stubref << " ";
          }
        }
        std::cout << ")  ";
      }
      std::cout << std::endl;
    }
  }
  
  return output_road;
}

























