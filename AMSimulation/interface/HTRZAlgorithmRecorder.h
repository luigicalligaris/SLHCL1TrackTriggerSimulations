#ifndef AMSimulation_HTRZAlgorithmRecorder_h_
#define AMSimulation_HTRZAlgorithmRecorder_h_

#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/HTRZAlgorithm.h"

#include "boost/multi_array.hpp"

#include "Rtypes.h"


class TDirectory;
class TFile;
class TTree;
class TH2C;
class TH2D;

namespace slhcl1tt
{
  class HTRZAlgorithm;

  class HTRZAlgorithmRecorder
  {
  public:
    HTRZAlgorithmRecorder (const char* const out_filename, HTRZAlgorithmConfig const& monitored_algo_config);
    ~HTRZAlgorithmRecorder();
    
    void Snapshot(HTRZAlgorithm& monitored_algo, unsigned const event_number, unsigned const tower_number, unsigned const road_number, unsigned const pattern_number);
    
  private:
    HTRZAlgorithmRecorder();
    
  private:
    UInt_t htmatrix_event_number_   ;
    UInt_t htmatrix_road_number_    ;
    UInt_t htmatrix_tower_number_   ;
    UInt_t htmatrix_pattern_number_ ;
    
    TFile*          tree_file_;
    TTree*          tree_;
    
    TH2C* htmatrix_layercount_allstubs_;
    TH2C* htmatrix_layercount_psstubs_ ;
    TH2C* htmatrix_layercount_2sstubs_ ;
    TH2C* htmatrix_activated_          ;
    
    TH2D* htmatrix_sums_layercount_allstubs_;
    TH2D* htmatrix_sums_layercount_psstubs_ ;
    TH2D* htmatrix_sums_layercount_2sstubs_ ;
    TH2D* htmatrix_sums_activated_          ;
  };
}

#endif