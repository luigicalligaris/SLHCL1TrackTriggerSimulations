#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/HTRZRoadFilter.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/TTRoadReader.h"

using namespace slhcl1tt;

// These below are set to private in class declaration
HTRZRoadFilter::HTRZRoadFilter() :
             mode_ (NULL_ALGO           ),
          verbose_ (0                   ),
        maxEvents_ (0                   ),
               po_ (ProgramOption()     ),
   readRoadPrefix_ ("AMTTRoads_"        ),
   readRoadSuffix_ (""                  ),
  writeRoadPrefix_ ("HTRZFilteredRoads_"),
  writeRoadSuffix_ (""                  ) {
  
}

HTRZRoadFilter::HTRZRoadFilter(HTRZRoadFilter const& rhs) :
             mode_ (rhs.           mode_),
          verbose_ (rhs.        verbose_),
        maxEvents_ (rhs.      maxEvents_),
               po_ (rhs.             po_),
   readRoadPrefix_ (rhs. readRoadPrefix_),
   readRoadSuffix_ (rhs. readRoadSuffix_),
  writeRoadPrefix_ (rhs.writeRoadPrefix_),
  writeRoadSuffix_ (rhs.writeRoadSuffix_) {
  
}

HTRZRoadFilter::HTRZRoadFilter(HTRZRoadFilter&&      rhs) :
             mode_ (rhs.           mode_),
          verbose_ (rhs.        verbose_),
        maxEvents_ (rhs.      maxEvents_),
               po_ (rhs.             po_),
   readRoadPrefix_ (rhs. readRoadPrefix_),
   readRoadSuffix_ (rhs. readRoadSuffix_),
  writeRoadPrefix_ (rhs.writeRoadPrefix_),
  writeRoadSuffix_ (rhs.writeRoadSuffix_) {
  
}

HTRZRoadFilter& HTRZRoadFilter::operator=(HTRZRoadFilter const& rhs) {
             mode_ = rhs.           mode_;
          verbose_ = rhs.        verbose_;
        maxEvents_ = rhs.      maxEvents_;
               po_ = rhs.             po_;
   readRoadPrefix_ = rhs. readRoadPrefix_;
   readRoadSuffix_ = rhs. readRoadSuffix_;
  writeRoadPrefix_ = rhs.writeRoadPrefix_;
  writeRoadSuffix_ = rhs.writeRoadSuffix_;
  return *this;
}

HTRZRoadFilter& HTRZRoadFilter::operator=(HTRZRoadFilter&&      rhs) {
             mode_ = rhs.           mode_;
          verbose_ = rhs.        verbose_;
        maxEvents_ = rhs.      maxEvents_;
               po_ = rhs.             po_;
   readRoadPrefix_ = rhs. readRoadPrefix_;
   readRoadSuffix_ = rhs. readRoadSuffix_;
  writeRoadPrefix_ = rhs.writeRoadPrefix_;
  writeRoadSuffix_ = rhs.writeRoadSuffix_;
  return *this;
}



// The actual constructor we use
HTRZRoadFilter::HTRZRoadFilter(const ProgramOption& po) :
             mode_ (NULL_ALGO           ),
          verbose_ (po.verbose          ),
        maxEvents_ (po.maxEvents        ),
               po_ (po                  ),
   readRoadPrefix_ ("AMTTRoads_"        ),
   readRoadSuffix_ (""                  ),
  writeRoadPrefix_ ("HTRZFilteredRoads_"),
  writeRoadSuffix_ (""                  ) {
  
}

HTRZRoadFilter::~HTRZRoadFilter() {
  
}


int HTRZRoadFilter::run() {
  int exitcode = 0;
  Timing(1);

  exitcode = filterRoads(po_.input, po_.output);
  if (exitcode)
    return exitcode;

  Timing();

  return exitcode;
}

int HTRZRoadFilter::filterRoads(TString inputfilename, TString outputfilename) {

  // For reading
  TTRoadReader reader(verbose_);
  if (reader.init(inputfilename, readRoadPrefix_, readRoadSuffix_)) {
    std::cout << Error() << "Failed to initialize TTRoadReader." << std::endl;
    return 1;
  }

  // For writing
  TTRoadWriter writer(verbose_);
  if (writer.init(reader.getChain(), outputfilename, writeRoadPrefix_, writeRoadSuffix_)) {
    std::cout << Error() << "Failed to initialize TTRoadWriter." << std::endl;
    return 1;
  }


  // Monitoring counters for read events and kept events
  long int nRead = 0, nKept = 0;

  // Output vector
  std::vector<TTRoad> out_roads;
  out_roads.reserve(300);

  for (long long ievt = 0; ievt < maxEvents_; ++ievt) {
    if (reader.loadTree(ievt) < 0)
      break;

    reader.getEntry(ievt);

    const unsigned nroads = reader.vr_patternRef->size();

    if (verbose_ > 1 && ievt % 100 == 0)
      std::cout << Debug() << Form("... Processing event: %7lld, fitting: %7ld", ievt, nKept) << std::endl;

    if (verbose_ > 2)
      std::cout << Debug() << "... evt: " << ievt << " # roads: " << nroads << std::endl;

    // Optimization: short circuit if no roads
    if (!nroads) {  // skip if no road
      writer.fill(std::vector<TTRoad>());
      ++nRead;
      continue;
    }

    // Not needed
    //out_roads.clear();

    // Process the input roads
    for (unsigned iroad = 0; iroad < nroads; ++iroad) {

      if (iroad >= (unsigned) po_.maxRoads)
        break;

      // Instantiate output TTRoad
      TTRoad out_road;

      switch (mode_)
      {
        case FILTER_ROADS:
        
        case NULL_ALGO:
        default:
        {
          out_road.patternRef    = reader.vr_patternRef    ->at(iroad);
          out_road.tower         = reader.vr_tower         ->at(iroad);
          out_road.nstubs        = reader.vr_nstubs        ->at(iroad);
          out_road.patternInvPt  = reader.vr_patternInvPt  ->at(iroad);
          out_road.superstripIds = reader.vr_superstripIds ->at(iroad);
          out_road.stubRefs      = reader.vr_stubRefs      ->at(iroad);
        }
        break;
        
      }
      
      out_roads.push_back(out_road);
    }

    writer.fill(out_roads);
    ++nRead;
  }


  if (nRead == 0) {
    std::cout << Error() << "Failed to read any event." << std::endl;
    return 1;
  }

  if (verbose_)
    std::cout << Info() << Form("Read: %7ld, triggered: %7ld", nRead, nKept) << std::endl;

  long long const nentries = writer.writeTree();
  assert(nentries == nRead);

  return 0;
}