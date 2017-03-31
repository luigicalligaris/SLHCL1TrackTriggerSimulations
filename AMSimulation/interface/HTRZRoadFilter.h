#ifndef AMSimulation_HTRZRoadFilter_h_
#define AMSimulation_HTRZRoadFilter_h_

#include "TString.h"

#include "SLHCL1TrackTriggerSimulations/AMSimulationDataFormats/interface/Pattern.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/Helper.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/ProgramOption.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/TriggerTowerMap.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/SuperstripArbiter.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/AssociativeMemory.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/HitBuffer.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/ModuleOverlapMap.h"



class HTRZRoadFilter {
  public:
    HTRZRoadFilter(const slhcl1tt::ProgramOption& po);
    ~HTRZRoadFilter();

    // Main driver
    int run();


  private:
    // Prevent default constructors from being called (refer to rule-of-five)
    HTRZRoadFilter();
    HTRZRoadFilter(HTRZRoadFilter const& rhs);
    HTRZRoadFilter(HTRZRoadFilter&&      rhs);
    HTRZRoadFilter& operator=(HTRZRoadFilter const& rhs);
    HTRZRoadFilter& operator=(HTRZRoadFilter&&      rhs);

    // Do road filtering and write surviving roads
    int filterRoads(TString src, TString out);

  private:
    // Configurations
    TString  readRoadPrefix_ ;
    TString  readRoadSuffix_ ;
    TString writeRoadPrefix_ ;
    TString writeRoadSuffix_ ;

    // Program options
    slhcl1tt::ProgramOption po_;
};

#endif
