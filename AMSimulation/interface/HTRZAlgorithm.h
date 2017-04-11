#ifndef AMSimulation_HTRZAlgorithm_h_
#define AMSimulation_HTRZAlgorithm_h_

#include "SLHCL1TrackTriggerSimulations/AMSimulationDataFormats/interface/TTRoad.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/TTRoadReader.h"

#include <vector>

#include "boost/multi_array.hpp"



namespace slhcl1tt
{
//   struct HTMatrix2DConfig
//   {
//     unsigned n_bins_x;
//     unsigned n_bins_y;
//     unsigned n_slots_per_bin;
//   };
//   
//   struct HTMatrix2DStubData
//   {
//     unsigned stubref;
//     unsigned layernum;
//   };
//   
//   class HTMatrix2D
//   {
//   private:
//     inline std::size_t cell_number(unsigned const bin_x, unsigned const bin_y)
//     {
//       // Use row-major addressing:
//       // 123
//       // 456
//       // 789
//       
//       assert(bin_x < n_bins_x_);
//       assert(bin_y < n_bins_y_);
//       
//       return bin_x + bin_y * n_bins_x_;
//     }
//     
//     inline bool cell_full(unsigned const bin_x, unsigned const bin_y)
//     {
//       return (*(matrix_occu_ + cell_number(bin_x, bin_y))) == n_slots_per_cell_;
//     }
//     
//     inline std::size_t free_slot_idx(unsigned const bin_x, unsigned const bin_y)
//     {
//       return cell_number(bin_x, bin_y) * n_slots_per_cell + occupancy(bin_x, bin_y);
//     }
//     
//     
//   public:
//     HTMatrix(HTMatrix2DConfig const& config):
//       n_bins_x_        (config.n_bins_x        ),
//       n_bins_y_        (config.n_bins_y        ),
//       n_slots_per_cell_(config.n_slots_per_cell),
//       matrix_occu_     (new unsigned[config.n_bins_x * config.n_bins_y])
//       matrix_data_     (new HTMatrix2DStubData[config.n_bins_x * config.n_bins_y * config.n_slots_per_bin])
//     {
//       
//     }
//     
//     ~HTMatrix();
//     {
//       delete[] matrix_data_;
//       matrix_data_ = nullptr;
//       
//       delete[] matrix_occu_;
//       matrix_occu_ = nullptr;
//     }
//     
//     unsigned occupancy(unsigned const bin_x, unsigned const bin_y)
//     {
//       return *(matrix_occu_ + cell_number(bin_x, bin_y));
//     }
//     
//     void reset_matrix()
//     {
//       for (size_t idx : matrix_reset_)
//       {
//         *(matrix_occu_ + idx) = 0u;
//       }
//     }
//     
//     bool store(unsigned const stub_ref)
//     {
//       if (occupancy()
//     }
//     
//   private:
//     HTMatrix(){}
//     
//     unsigned n_bins_x_;
//     unsigned n_bins_y_;
//     unsigned n_slots_per_cell_;
//     
//     unsigned* matrix_occu_ ;
//     HTMatrix2DStubData* matrix_data_ ;
//     std::vector<std::size_t> matrix_reset_;
//   };
  
  
  
  
  
  // 2D HT matrix or 1D histogram?
  enum HTRZAlgorithm_mode_t {NULL_ALGO, HTRZ_2D_COTANTHETA_ZT};
  enum HTRZAlgorithm_stub_accept_policy_t {LOOSE_ALL_NEIGHBOURS, MEDIUM_NEAR_NEIGHBOUR, TIGHT_NO_NEIGHBOURS};

  struct HTRZAlgorithmConfig
  {
    HTRZAlgorithm_mode_t               mode                 ;
    HTRZAlgorithm_stub_accept_policy_t stub_accept_policy   ;
    bool                               color_output         ;
    int                                verbose              ;
    unsigned                           nlayers              ;
    unsigned                           nbins_zT             ;
    unsigned                           nbins_cotantheta     ;
    unsigned                           threshold_all_layers ;
    unsigned                           threshold_ps_layers  ;
    double                             min_zT               ;
    double                             max_zT               ;
    double                             min_cotantheta       ;
    double                             max_cotantheta       ;
    double                             t_radius             ;
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
    friend class HTRZAlgorithmRecorder;
    
    struct HTRZCell
    {
      std::vector< unsigned > stubrefs;
    };
    
    HTRZAlgorithm();
    HTRZAlgorithm(HTRZAlgorithm const& rhs);
    HTRZAlgorithm(HTRZAlgorithm&&      rhs);
    HTRZAlgorithm& operator=(HTRZAlgorithm const& rhs);
    HTRZAlgorithm& operator=(HTRZAlgorithm&&      rhs);
    ~HTRZAlgorithm();

    HTRZAlgorithm(HTRZAlgorithmConfig const& config);
    void LoadConfig(HTRZAlgorithmConfig const& config);

    slhcl1tt::TTRoad Filter(slhcl1tt::TTRoad const& input_road, slhcl1tt::TTRoadReader const& reader);

  private:
    inline double stub_and_cotantheta_to_zT(double const t_radius, double const zstub, double const rstub, double const cotantheta) const noexcept
    {
      return zstub - (rstub - t_radius) * cotantheta;
    }

    inline double stub_and_zT_to_cotantheta(double const t_radius, double const zstub, double const rstub, double const zT) const noexcept
    {
      assert(rstub != t_radius);
      return (zstub - zT) / (rstub - t_radius);
    }

    void RegenerateBoundariesZT();
    void RegenerateBoundariesCotanTheta();
    inline void RegenerateBoundaries(){RegenerateBoundariesCotanTheta(); RegenerateBoundariesZT();}

    void RegenerateMatrix2D();
    void ResetMatrix2D();
    
    unsigned int count_stubs_all_layers(unsigned const iCot, unsigned const iZT);
    unsigned int count_stubs_ps_layers (unsigned const iCot, unsigned const iZT);
    bool majority_all_layers(unsigned const iCot, unsigned const iZT);
    bool majority_ps_layers (unsigned const iCot, unsigned const iZT);
    bool majority(unsigned const iCot, unsigned const iZT)
    {
      return majority_all_layers(iCot, iZT) && majority_ps_layers(iCot, iZT);
    }

  private:
    HTRZAlgorithmConfig config_;

    std::vector<double> ht_bounds_cotantheta_;
    std::vector<double> ht_bounds_zT_;
    
    
    boost::multi_array<HTRZCell, 3> ht_matrix_stubs_;
  };
}

#endif