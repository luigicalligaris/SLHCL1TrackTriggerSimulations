#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/HTRZRoadFilter.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/HTRZAlgorithm.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/TTRoadReader.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/Helper.h"

using namespace slhcl1tt;



// These below are set to private in class declaration, you won't use them
HTRZRoadFilter::HTRZRoadFilter() :
   readRoadPrefix_ ("AMTTRoads_"        ),
   readRoadSuffix_ (""                  ),
  writeRoadPrefix_ ("HTRZFilteredRoads_"),
  writeRoadSuffix_ (""                  ),
               po_ (ProgramOption()     )
{
  
}

HTRZRoadFilter::HTRZRoadFilter(HTRZRoadFilter const& rhs) :
   readRoadPrefix_ (rhs. readRoadPrefix_),
   readRoadSuffix_ (rhs. readRoadSuffix_),
  writeRoadPrefix_ (rhs.writeRoadPrefix_),
  writeRoadSuffix_ (rhs.writeRoadSuffix_),
               po_ (rhs.             po_)
{
  
}

HTRZRoadFilter::HTRZRoadFilter(HTRZRoadFilter&&      rhs) :
   readRoadPrefix_ (rhs. readRoadPrefix_),
   readRoadSuffix_ (rhs. readRoadSuffix_),
  writeRoadPrefix_ (rhs.writeRoadPrefix_),
  writeRoadSuffix_ (rhs.writeRoadSuffix_),
               po_ (rhs.             po_)
{
  
}

HTRZRoadFilter& HTRZRoadFilter::operator=(HTRZRoadFilter const& rhs)
{
   readRoadPrefix_ = rhs. readRoadPrefix_;
   readRoadSuffix_ = rhs. readRoadSuffix_;
  writeRoadPrefix_ = rhs.writeRoadPrefix_;
  writeRoadSuffix_ = rhs.writeRoadSuffix_;
               po_ = rhs.             po_;
  return *this;
}

HTRZRoadFilter& HTRZRoadFilter::operator=(HTRZRoadFilter&&      rhs)
{
   readRoadPrefix_ = rhs. readRoadPrefix_;
   readRoadSuffix_ = rhs. readRoadSuffix_;
  writeRoadPrefix_ = rhs.writeRoadPrefix_;
  writeRoadSuffix_ = rhs.writeRoadSuffix_;
               po_ = rhs.             po_;
  return *this;
}



// The actual constructor we use
HTRZRoadFilter::HTRZRoadFilter(const ProgramOption& po) :
   readRoadPrefix_ ("AMTTRoads_"        ),
   readRoadSuffix_ (""                  ),
  writeRoadPrefix_ ("HTRZFilteredRoads_"),
  writeRoadSuffix_ (""                  ),
               po_ (po                  )
{
  
}

HTRZRoadFilter::~HTRZRoadFilter()
{
  
}


int HTRZRoadFilter::run()
{
  int exitcode = 0;
  Timing(1);

  exitcode = filterRoads(po_.input, po_.output);
  if (exitcode)
    return exitcode;

  Timing();

  return exitcode;
}

int HTRZRoadFilter::filterRoads(TString inputfilename, TString outputfilename)
{
  // For reading
  TTRoadReader reader(po_.verbose);
  if (reader.init(inputfilename, readRoadPrefix_, readRoadSuffix_))
  {
    std::cout << Error() << "Failed to initialize TTRoadReader." << std::endl;
    return 1;
  }

  // For writing
  TTRoadWriter writer(po_.verbose);
  if (writer.init(reader.getChain(), outputfilename, writeRoadPrefix_, writeRoadSuffix_))
  {
    std::cout << Error() << "Failed to initialize TTRoadWriter." << std::endl;
    return 1;
  }


  // Monitoring counters for read events and kept events
  long int nEvtRead = 0, nEvtKept = 0, nRoadRead = 0, nRoadKept = 0;

  // Output vector
  std::vector<TTRoad> out_roads;
  out_roads.reserve(300);



  // Set-up HT algo

  HTRZAlgorithmConfig algo_config;
  
  
  if      (po_.htRZMode == "HTRZ_2D_COTANTHETA_Z0")
    algo_config.mode = HTRZ_2D_COTANTHETA_Z0 ;
  else if (po_.htRZMode == "HTRZ_1D_COTANTHETA")
    algo_config.mode = HTRZ_1D_COTANTHETA;
  else // (po_.htRZMode == "NULL_ALGO")
    algo_config.mode = NULL_ALGO;
  
  if      (po_.htRZStubAcceptPolicy == "MEDIUM_NEAR_NEIGHBOUR")
    algo_config.stub_accept_policy = MEDIUM_NEAR_NEIGHBOUR;
  else if (po_.htRZStubAcceptPolicy == "TIGHT_NO_NEIGHBOURS")
    algo_config.stub_accept_policy = TIGHT_NO_NEIGHBOURS;
  else // (po_.htRZStubAcceptPolicy == "LOOSE_ALL_NEIGHBOURS")
    algo_config.stub_accept_policy = LOOSE_ALL_NEIGHBOURS ;
  
  
//   algo_config.verbose               =                     1 ;
//   algo_config.max_z0                =                 +15.0 ;
//   algo_config.min_z0                =                 -15.0 ;
//           algo_config.max_cotantheta        =                 +13.5 ;
//           algo_config.min_cotantheta        =                 -13.5 ;
//   algo_config.max_cotantheta        =                  +1.5 ;
//   algo_config.min_cotantheta        =                  -1.0 ;
//   algo_config.nbins_z0              =                     4 ;
//   algo_config.nbins_cotantheta      =                    16 ;
//   algo_config.threshold_all_layers  =                     4 ;
//   algo_config.threshold_ps_layers   =                     1 ;
  
  algo_config.color_output          = IsColorEnabled();
  algo_config.verbose               = po_.verbose               ;
  algo_config.nbins_z0              = po_.htRZZ0Bins            ;
  algo_config.nbins_cotantheta      = po_.htRZCotanThetaBins    ;
  algo_config.threshold_all_layers  = po_.htRZThresholdLayerAll ;
  algo_config.threshold_ps_layers   = po_.htRZThresholdLayerPS  ;
  algo_config.min_cotantheta        = po_.htRZCotanThetaMin     ;
  algo_config.max_cotantheta        = po_.htRZCotanThetaMax     ;
  algo_config.min_z0                = po_.htRZZ0Min             ;
  algo_config.max_z0                = po_.htRZZ0Max             ;
  
  
  HTRZAlgorithm algo(algo_config);


  // Run on the events in tree
  for (long long ievt = 0; ievt < po_.maxEvents; ++ievt)
  {
    if (reader.loadTree(ievt) < 0)
      break;

    ++nEvtRead;

    reader.getEntry(ievt);

    const unsigned nroads = reader.vr_patternRef->size();

    if (po_.verbose >= 2 && ievt % 100 == 0)
      std::cout << Debug() << Form("... Processing event: %7lld", ievt) << std::endl;

    if (po_.verbose >= 3)
      std::cout << Debug() << "... evt: " << ievt << " # roads: " << nroads << std::endl;

    // Optimization: short circuit if no roads
    if (!nroads)
    {
      writer.fill(std::vector<TTRoad>());
      ++nEvtKept;
      continue;
    }
  
    out_roads.clear();
  
    // Process the input roads
    for (unsigned iroad = 0; iroad < nroads; ++iroad)
    {
      if (iroad >= (unsigned) po_.maxRoads)
        break;
      
      ++nRoadRead;
      
      TTRoad in_road;
      
      in_road.patternRef    = reader.vr_patternRef    ->at(iroad);
      in_road.tower         = reader.vr_tower         ->at(iroad);
      in_road.nstubs        = reader.vr_nstubs        ->at(iroad);
      in_road.patternInvPt  = reader.vr_patternInvPt  ->at(iroad);
      in_road.superstripIds = reader.vr_superstripIds ->at(iroad);
      in_road.stubRefs      = reader.vr_stubRefs      ->at(iroad);
      
      TTRoad out_road = algo.Filter(in_road, reader);
      
      if (out_road.nstubs > 0)
      {
        ++nRoadKept;
        out_roads.push_back(out_road);
      }
    }
    
    if (!out_roads.empty())
    {
      ++nEvtKept;
      writer.fill(out_roads);
    }
  }
  
  
  
  if (nEvtRead == 0)
  {
    std::cout << Error() << "Failed to read any event." << std::endl;
    return 1;
  }
  
  if (po_.verbose >= 0)
  {
    std::cout << Info() << Form("Events Read: %7ld, triggered: %7ld   Roads Read: %7ld, passing: %7ld", nEvtRead, nEvtKept, nRoadRead, nRoadKept) << "   Road eff: ";
    if (nRoadRead > 0)
      std::cout << Form("%7f", double(nRoadKept)/double(nRoadRead)) << std::endl;
    else
      std::cout << "N/A" << std::endl;
  }
  
  long long const nentries = writer.writeTree();
//   assert(nentries == nEvtRead);
  
  return 0;
}