#ifndef AMSimulation_ParameterDuplicateRemoval_h_
#define AMSimulation_ParameterDuplicateRemoval_h_

#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/TTTrackReader.h"
#include <vector>
#include <TVector2.h>

namespace slhcl1tt{
  class ParameterDuplicateRemoval{
    public:
      ParameterDuplicateRemoval() {}
      ~ParameterDuplicateRemoval() {}

      void ReduceTracks(std::vector<TTTrack2>& Tracks);
  };

  //structure containing necessary original Track information
  struct TrackCloud{
    std::vector<unsigned> TrackIterators;
    unsigned  BestTrack;
    bool BestRoadCategory; //false = 5/x, true = 6/6
    float BestChi2;
    float BestEta;
    float BestPhi;
    float BestPt; //test

    bool Add(unsigned TrackIterator_, const std::vector<unsigned> &combination_, float Chi2_, float Eta_, float Phi_, float Pt_){ //test
      bool ContainsNoDummy=true;
      for(unsigned i=0; i<combination_.size(); ++i){
        if(combination_[i]==999999999) ContainsNoDummy=false;
      }
      if(TrackIterators.size()==0){
        TrackIterators.push_back(TrackIterator_);
        BestChi2=Chi2_;
        BestEta=Eta_;
        BestPhi=Phi_;
        BestPt=Pt_; //test
        BestRoadCategory=ContainsNoDummy;
        BestTrack=TrackIterator_;
        return true;
      }
      else if(fabs(Eta_-BestEta)<0.0601 && fabs(TVector2::Phi_mpi_pi(Phi_-BestPhi))<0.00576 && ((BestPt>100. && Pt_>100.) || fabs(BestPt-Pt_)<10.)){ //test
        if(Chi2_<BestChi2){    
          BestChi2=Chi2_;
          BestEta=Eta_;
          BestPhi=Phi_;
	  BestPt=Pt_; //test
	  BestRoadCategory=ContainsNoDummy;
	  BestTrack=TrackIterator_;
        }
        TrackIterators.push_back(TrackIterator_);
        return true;
      }
      else return false;
    }
  };
}

#endif
