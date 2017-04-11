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


#define HTRZALGORITHM_DEFAULT_MODE                        HTRZ_2D_COTANTHETA_ZT
#define HTRZALGORITHM_DEFAULT_STUB_ACCEPT_POLICY          LOOSE_ALL_NEIGHBOURS
#define HTRZALGORITHM_DEFAULT_COLOR_OUTPUT                true
#define HTRZALGORITHM_DEFAULT_VERBOSE                     -1
#define HTRZALGORITHM_DEFAULT_NLAYERS                     6
#define HTRZALGORITHM_DEFAULT_NBINS_ZT                    8
#define HTRZALGORITHM_DEFAULT_NBINS_COTANTHETA            8
#define HTRZALGORITHM_DEFAULT_NBINS_THRESHOLD_ALL_LAYERS  4
#define HTRZALGORITHM_DEFAULT_NBINS_THRESHOLD_PS_LAYERS   2
#define HTRZALGORITHM_DEFAULT_MAX_ZT                       15.0
#define HTRZALGORITHM_DEFAULT_MIN_ZT                      -15.0
#define HTRZALGORITHM_DEFAULT_MAX_COTANTHETA               13.5
#define HTRZALGORITHM_DEFAULT_MIN_COTANTHETA              -13.5
#define HTRZALGORITHM_DEFAULT_T_RADIUS                      0.0


using namespace slhcl1tt;

HTRZAlgorithm::HTRZAlgorithm():
  config_{
    HTRZALGORITHM_DEFAULT_MODE                      ,
    HTRZALGORITHM_DEFAULT_STUB_ACCEPT_POLICY        ,
    HTRZALGORITHM_DEFAULT_COLOR_OUTPUT              ,
    HTRZALGORITHM_DEFAULT_VERBOSE                   ,
    HTRZALGORITHM_DEFAULT_NLAYERS                   ,
    HTRZALGORITHM_DEFAULT_NBINS_ZT                  ,
    HTRZALGORITHM_DEFAULT_NBINS_COTANTHETA          ,
    HTRZALGORITHM_DEFAULT_NBINS_THRESHOLD_ALL_LAYERS,
    HTRZALGORITHM_DEFAULT_NBINS_THRESHOLD_PS_LAYERS ,
    HTRZALGORITHM_DEFAULT_MIN_ZT                    ,
    HTRZALGORITHM_DEFAULT_MAX_ZT                    ,
    HTRZALGORITHM_DEFAULT_MIN_COTANTHETA            ,
    HTRZALGORITHM_DEFAULT_MAX_COTANTHETA            ,
    HTRZALGORITHM_DEFAULT_T_RADIUS                  
  }
{
  assert(config_.max_zT         >= config_.min_zT        );
  assert(config_.max_cotantheta >= config_.min_cotantheta);

  RegenerateBoundaries();
  RegenerateMatrix2D();
}


HTRZAlgorithm::HTRZAlgorithm(HTRZAlgorithm const& rhs):
  config_(rhs.config_)
{
  assert(config_.max_zT         >= config_.min_zT        );
  assert(config_.max_cotantheta >= config_.min_cotantheta);

  ht_bounds_cotantheta_ = rhs.ht_bounds_cotantheta_;
  ht_bounds_zT_         = rhs.ht_bounds_zT_        ;
}

HTRZAlgorithm::HTRZAlgorithm(HTRZAlgorithm&&      rhs):
  config_(rhs.config_)
{
  assert(config_.max_zT         >= config_.min_zT        );
  assert(config_.max_cotantheta >= config_.min_cotantheta);

  ht_bounds_cotantheta_ = std::move(rhs.ht_bounds_cotantheta_);
  ht_bounds_zT_         = std::move(rhs.ht_bounds_zT_        );
}

HTRZAlgorithm& HTRZAlgorithm::operator=(HTRZAlgorithm const& rhs)
{
  config_ = rhs.config_;

  assert(config_.max_zT         >= config_.min_zT        );
  assert(config_.max_cotantheta >= config_.min_cotantheta);

  ht_bounds_cotantheta_ = rhs.ht_bounds_cotantheta_;
  ht_bounds_zT_         = rhs.ht_bounds_zT_        ;

  return *this;
}

HTRZAlgorithm& HTRZAlgorithm::operator=(HTRZAlgorithm&&      rhs)
{
  config_ = rhs.config_;

  assert(config_.max_zT         >= config_.min_zT        );
  assert(config_.max_cotantheta >= config_.min_cotantheta);

  ht_bounds_cotantheta_ = std::move(rhs.ht_bounds_cotantheta_);
  ht_bounds_zT_         = std::move(rhs.ht_bounds_zT_        );

  return *this;
}

HTRZAlgorithm::~HTRZAlgorithm()
{
  
}


HTRZAlgorithm::HTRZAlgorithm(HTRZAlgorithmConfig const& config):
  config_(config)
{
  assert(config_.max_zT         >= config_.min_zT        );
  assert(config_.max_cotantheta >= config_.min_cotantheta);

  RegenerateBoundaries();
  RegenerateMatrix2D();
}

void HTRZAlgorithm::LoadConfig(HTRZAlgorithmConfig const& config)
{
  config_ = config;

  assert(config_.max_zT         >= config_.min_zT        );
  assert(config_.max_cotantheta >= config_.min_cotantheta);

  RegenerateBoundaries();
  RegenerateMatrix2D();
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

void HTRZAlgorithm::RegenerateBoundariesZT()
{
  for (unsigned iZT = 0; iZT < config_.nbins_zT + 1; ++iZT)
  {
    ht_bounds_zT_.push_back( config_.min_zT * (double(config_.nbins_zT - iZT)/double(config_.nbins_zT)) + config_.max_zT * (double(iZT)/double(config_.nbins_zT)) );
  }
  
  if (config_.verbose >= 4)
  {
    std::cout << "New bounds in z0 = [";
    for (auto val : ht_bounds_zT_)
    {
      std::cout << val << ",";
    }
    std::cout << "]" << std::endl;
  }
}



void HTRZAlgorithm::ResetMatrix2D()
{
  for (unsigned iCot = 0; iCot < config_.nbins_cotantheta; ++iCot)
    for (unsigned iZT = 0; iZT < config_.nbins_zT; ++iZT)
      for (unsigned iLayer = 0; iLayer < config_.nlayers; ++iLayer)
      {
        ht_matrix_stubs_[iCot][iZT][iLayer].stubrefs.clear();
      }
  
}

void HTRZAlgorithm::RegenerateMatrix2D()
{
  ht_matrix_stubs_.resize(boost::extents[config_.nbins_cotantheta][config_.nbins_zT][config_.nlayers]);
  ResetMatrix2D();
}


unsigned int HTRZAlgorithm::count_stubs_all_layers(unsigned const iCot, unsigned const iZT)
{
  unsigned int count = 0;
  
  for (unsigned iLayer = 0; iLayer < config_.nlayers; ++iLayer)
  {
    if (!ht_matrix_stubs_[iCot][iZT][iLayer].stubrefs.empty())
    {
      ++count;
    }
  }
  
  return count;
}

unsigned int HTRZAlgorithm::count_stubs_ps_layers(unsigned const iCot, unsigned const iZT)
{
  unsigned int count = 0;
  
  for (unsigned iLayer = 0; iLayer < 3; ++iLayer)
  {
    if (!ht_matrix_stubs_[iCot][iZT][iLayer].stubrefs.empty())
    {
      ++count;
    }
  }
  
  return count;
}


bool HTRZAlgorithm::majority_all_layers(unsigned const iCot, unsigned const iZT)
{
  return count_stubs_all_layers(iCot, iZT) >= config_.threshold_all_layers;
}

bool HTRZAlgorithm::majority_ps_layers(unsigned const iCot, unsigned const iZT)
{
  return count_stubs_ps_layers(iCot, iZT) >= config_.threshold_ps_layers;
}



constexpr static const double lut_pixel_halfwidth_z_bylayer[] = {0.0722812068078, 0.0722812068078, 0.0722812068078, 2.5125007629395, 2.5125007629395, 2.5125007629395};


TTRoad HTRZAlgorithm::Filter(slhcl1tt::TTRoad const& input_road, TTRoadReader const& reader)
{
  ResetMatrix2D();
  
  TTRoad output_road;
  
  std::unordered_set<unsigned int> passing_stubrefs;
  
  switch (config_.mode)
  {
    case NULL_ALGO:
      output_road = input_road;
      break;
    
    
    case HTRZ_2D_COTANTHETA_ZT:
      {
        // There is 1 superstrip per layer
        for (unsigned iLayer = 0; iLayer < config_.nlayers; ++iLayer)
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
            std::vector< std::vector< unsigned int> > bound_iZT    (config_.nbins_cotantheta + 1,  std::vector<unsigned>());
            
            
            // Loop on the (config_.nbins_cotantheta + 1) bin edges
            for (unsigned iCot = 0; iCot < config_.nbins_cotantheta + 1; ++iCot)
            {
              double const zT_lo = stub_and_cotantheta_to_zT(config_.t_radius, stub_z - lut_pixel_halfwidth_z_bylayer[iLayer], stub_r, ht_bounds_cotantheta_[iCot]);
              double const zT_hi = stub_and_cotantheta_to_zT(config_.t_radius, stub_z + lut_pixel_halfwidth_z_bylayer[iLayer], stub_r, ht_bounds_cotantheta_[iCot]);
              
//               if (config_.verbose >= 5)
//               {
//                 std::cout << "zT_lo = " << std::fixed << std::setprecision(4) << std::showpoint << std::showpos << std::setw(10) << zT_lo << "  zT_hi = " << std::fixed << std::setprecision(4) << std::showpoint << std::showpos << std::setw(10) << zT_hi << std::endl;
//               }
              
              assert(zT_lo < zT_hi);
              
              if      (zT_hi <= ht_bounds_zT_.front())
              {
                // The z0 at this column falls out of the bottom of the matrix
                bound_accept[iCot] = false;
                
                if (config_.verbose >= 4)
                  std::cout << "layer " << iLayer << " stubRef " << stubRef << " cotan theta boundary " << iCot << " REJECTED (below HT bottom) zT_lo = " << std::fixed << std::setprecision(4) << std::showpoint << std::showpos << std::setw(10) << zT_lo 
                    <<                         "  zT_hi = " << std::fixed << std::setprecision(4) << std::showpoint << std::showpos << std::setw(10) << zT_hi << std::endl;
              }
              else if (zT_lo >= ht_bounds_zT_.back() )
              {
                // The z0 at this column falls out of the top of the matrix
                bound_accept[iCot] = false;
                
                if (config_.verbose >= 4)
                  std::cout << "layer " << iLayer << " stubRef " << stubRef << " cotan theta boundary " << iCot << " REJECTED (above HT top   ) zT_lo = " << std::fixed << std::setprecision(4) << std::showpoint << std::showpos << std::setw(10) << zT_lo 
                    <<                         "  zT_hi = " << std::fixed << std::setprecision(4) << std::showpoint << std::showpos << std::setw(10) << zT_hi << std::endl;
              }
              else
              {
                // We are inside the matrix in the vertical direction
                
                unsigned int const iZT_lo = zT_lo <= ht_bounds_zT_.front() ? 0             : std::distance(ht_bounds_zT_.begin(), std::lower_bound(ht_bounds_zT_.begin(), ht_bounds_zT_.end(), zT_lo) ) - 1;
                
                unsigned int const iZT_hi = zT_hi >= ht_bounds_zT_.back()  ? config_.nbins_zT - 1 : std::distance(ht_bounds_zT_.begin(), std::lower_bound(ht_bounds_zT_.begin(), ht_bounds_zT_.end(), zT_hi) ) - 1;
                
                if (config_.verbose >= 4)
                  std::cout << "layer " << iLayer << " stubRef " << stubRef << " cotan theta boundary " << iCot << " ACCEPTED zT_lo = " << zT_lo << "  zT_hi = " << zT_hi << "  iZT_lo = " << iZT_lo << "  iZT_hi = " << iZT_hi << std::endl;
                
                assert(iZT_lo <= iZT_hi);
                
                
                // This true statement can be changed to more complicated requirements.
                // For now the matrix boundaries (-15.0,15.0) are enough to guarantee 
                // the respect of the condition we care about: that the track candidate 
                // is compatible with a real track from collisions.
                
                bound_accept[iCot] = true;
                
                for (unsigned int iZT = iZT_lo; iZT <= iZT_hi; ++iZT)
                {
                  bound_iZT[iCot].push_back(iZT);
                }
                
                assert(!bound_iZT[iCot].empty());
              }
            }
            
            
            // Accept/reject the stubs based on the values at the boundaries and the 
            // set acceptance policy.
            for (unsigned iCot = 0; iCot < config_.nbins_cotantheta; ++iCot)
            {
              if (config_.verbose >= 4)
                std::cout << "layer " << iLayer << " stubRef " << stubRef << " cotan theta bin " << iCot;
              
              if (!bound_accept[iCot] && !bound_accept[iCot + 1])
              {
                // Both borders are invalid, skip the stub in this column
                if (config_.verbose >= 4)
                  std::cout << " both borders invalid --> reject " << std::endl;
                
                continue;
              }
              else if (bound_accept[iCot] && bound_accept[iCot + 1])
              {
                if (config_.verbose >= 4)
                  std::cout << " both borders valid --> accept Z0 bins ";
                
                std::set<unsigned> merged_iZT( bound_iZT[iCot    ].begin(), bound_iZT[iCot    ].end() );
                merged_iZT            .insert( bound_iZT[iCot + 1].begin(), bound_iZT[iCot + 1].end() );
                
                for (unsigned iZT : merged_iZT)
                {
                  if (config_.verbose >= 4)
                    std::cout << iZT << " ";
                  ht_matrix_stubs_[iCot][iZT][iLayer].stubrefs.push_back(stubRef);
                }
                
                if (config_.verbose >= 4)
                  std::cout << std::endl;
              }
              else
              {
                if (config_.verbose >= 4)
                  std::cout << " one border valid, another invalid --> loookup policy ";
                // One border is valid, the other is not. Follow the set acceptance policy.
                switch (config_.stub_accept_policy)
                {
                  case LOOSE_ALL_NEIGHBOURS:
                    {
                      if (config_.verbose >= 4)
                        std::cout << " --> add all neighbor Z0 bins ";
                      
                      std::set<unsigned> merged_iZT( bound_iZT[iCot    ].begin(), bound_iZT[iCot    ].end() );
                      merged_iZT            .insert( bound_iZT[iCot + 1].begin(), bound_iZT[iCot + 1].end() );
                      
                      unsigned iZT_end = *(merged_iZT.rbegin());
                      
                      for (unsigned iZT = *(merged_iZT.begin()); iZT < iZT_end; ++iZT)
                      {
                        if (config_.verbose >= 4)
                          std::cout << iZT << " ";
                        ht_matrix_stubs_[iCot][iZT][iLayer].stubrefs.push_back(stubRef);
                      }
                      
                      if (config_.verbose >= 4)
                        std::cout << std::endl;
                    }
                    break;
                  case MEDIUM_NEAR_NEIGHBOUR:
                    {
                      if (config_.verbose >= 4)
                        std::cout << " --> add near neighbor Z0 bins ";
                      
                      if (bound_accept[iCot])
                        for (unsigned iZT : bound_iZT[iCot    ])
                        {
                          if (config_.verbose >= 4)
                            std::cout << iZT << " ";
                          ht_matrix_stubs_[iCot][iZT][iLayer].stubrefs.push_back(stubRef);
                        }
                      else
                        for (unsigned iZT : bound_iZT[iCot + 1])
                        {
                          if (config_.verbose >= 4)
                            std::cout << iZT << " ";
                          ht_matrix_stubs_[iCot][iZT][iLayer].stubrefs.push_back(stubRef);
                        }
                      
                      if (config_.verbose >= 4)
                        std::cout << std::endl;
                    }
                    break;
                  case TIGHT_NO_NEIGHBOURS:
                    {
                      if (config_.verbose >= 4)
                        std::cout << " --> none added " << std::endl;
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
          std::cout << std::endl << "HT Matrix of road from pattern "<< input_road.patternRef << ":" << std::endl;
          
          for (unsigned int iZT = 0; iZT < config_.nbins_zT; ++iZT)
          {
            for (unsigned int iCot = 0; iCot < config_.nbins_cotantheta; ++iCot)
            {
              std::string color_str_begin = "";
              std::string color_str_end   = "";
              
              unsigned const count_all_lay = count_stubs_all_layers(iCot, iZT);
              unsigned const count_ps_lay  = count_stubs_ps_layers (iCot, iZT);
              
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
              
              std::cout << color_str_begin << count_all_lay << color_str_end;
            }
            
            std::cout << std::endl;
          }
        }
        
        // Loop over the matrix cells and take note of stubrefs in cells activated by the majority logic
        for (unsigned int iCot = 0; iCot < config_.nbins_cotantheta; ++iCot)
          for (unsigned int iZT = 0; iZT < config_.nbins_zT; ++iZT)
            if ( majority_all_layers(iCot, iZT))
            {
              for (unsigned iLayer = 0; iLayer < config_.nlayers; ++iLayer)
                for(unsigned int const stubref : ht_matrix_stubs_[iCot][iZT][iLayer].stubrefs)
                {
                  passing_stubrefs.insert(stubref);
                }
            }
        
      }
      break;
    
    
    default:
      break;
  }
  
  output_road.patternRef    = input_road.patternRef    ;
  output_road.tower         = input_road.tower         ;
  output_road.patternInvPt  = input_road.patternInvPt  ;
  output_road.superstripIds = input_road.superstripIds ;
  
  for (unsigned iLayer = 0; iLayer < config_.nlayers; ++iLayer)
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
      for (unsigned iLayer = 0; iLayer < config_.nlayers; ++iLayer)
      {
        std::cout << "Passing stubrefs in layer " << iLayer << ": ";
        for (unsigned int const in_stubref : input_road.stubRefs[iLayer])
        {
          if (passing_stubrefs.find(in_stubref) != passing_stubrefs.end())
          {
            std::cout << in_stubref << " ";
          }
        }
        std::cout << std::endl;
      }
      
      
      std::cout << "Failing stubrefs in: ";
      for (unsigned iLayer = 0; iLayer < config_.nlayers; ++iLayer)
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

























