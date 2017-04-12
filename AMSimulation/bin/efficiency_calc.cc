#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/TTRoadReader.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/TTStubReader.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/TTStubPlusTPReader.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/TTTrackReader.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/Helper.h"

#include <cmath>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TDirectoryFile.h"


std::vector<double> generateLogarithmicallySpacedBinBoundaries(unsigned const nbins, double const xlow, double const xhi)
{
  using namespace std;
  
  assert(nbins > 0);
  assert(xhi  > 0.0);
  assert(xlow > 0.0);
  assert(xlow < xhi);
  
  vector<double> ret;
  
  double const log_xlow = log(xlow);
  double const log_xhi  = log(xhi );
  
  for (int iBin = 0; iBin < int(nbins + 1u); ++iBin)
  {
    ret.push_back( exp(log_xlow + iBin * (log_xhi - log_xlow) / double(nbins)) );
  }
  
  return ret;
}


int main(int argc, char* argv[])
{
  using namespace std;
  using namespace slhcl1tt;
  
  const string inputfilename(argv[1]);
  const string outputfilename(argv[2]);
  
  const string denomRoadPrefix ("AMTTRoads_");
  const string denomRoadSuffix ("");
  const string numerRoadPrefix ("HTRZFilteredRoads_");
  const string numerRoadSuffix ("");
  
  int const verbosity = 1;
  long long const neventsmax = 100000;
  
  
  
  // For reading
  TTRoadReader       denomReader(verbosity);
  TTRoadReader       numerReader(verbosity);
  TTStubPlusTPReader tpartReader(verbosity);
  
  
  if (denomReader.init(inputfilename, denomRoadPrefix, denomRoadSuffix))
  {
    std::cout << Error() << "Failed to initialize TTRoadReader for denominator road collection." << std::endl;
    return 1;
  }
  
  if (numerReader.init(inputfilename, numerRoadPrefix, numerRoadSuffix))
  {
    std::cout << Error() << "Failed to initialize TTRoadReader for numerator road collection." << std::endl;
    return 1;
  }
  
  if (tpartReader.init(inputfilename))
  {
    std::cout << Error() << "Failed to initialize TTStubPlusTPReader." << std::endl;
    return 1;
  }
  
  
  vector<double> pt_log_binning = generateLogarithmicallySpacedBinBoundaries(20, 2.0, 100.0);
  
//   TEfficiency tpeff_pt     ("tpeff_pt"     ,"HT r-z filter efficiency for TPs;p_{T};#epsilon",200, 0.0 ,100.0);
  TEfficiency tpeff_pt     ("tpeff_pt"     ,"HT r-z filter efficiency for TPs;p_{T};#epsilon",pt_log_binning.size() - 1, pt_log_binning.data());
  TEfficiency tpeff_eta    ("tpeff_eta"    ,"HT r-z filter efficiency for TPs;#eta;#epsilon" ,100,-2.5 ,  2.5);
  
//   TEfficiency roadeff_pt     ("roadeff_pt"     ,"HT r-z filter efficiency for single roads associated to good post-AM TPs;p_{T};#epsilon",200, 0.0 ,100.0);
  TEfficiency roadeff_pt     ("roadeff_pt"     ,"HT r-z filter efficiency for single roads associated to good post-AM TPs;p_{T};#epsilon",pt_log_binning.size() - 1, pt_log_binning.data());
  TEfficiency roadeff_eta    ("roadeff_eta"    ,"HT r-z filter efficiency for single roads associated to good post-AM TPs;#eta;#epsilon" ,100,-2.5 ,  2.5);
  
  TH1F     nstubs_denom       ("nstubs_denom"       ,"Number of stubs per road (denominator);N_{stubs};N_{roads}",15, 0.0, 15.0);
  TH1F     nstubs_numer       ("nstubs_numer"       ,"Number of stubs per road (numerator);N_{stubs};N_{roads}"  ,15, 0.0, 15.0);
  TH1F     fstubs_denom       ("fstubs_denom"       ,"Fraction of stubs per road (denominator);N_{stubs};f_{roads}",15, 0.0, 15.0);
  TH1F     fstubs_numer       ("fstubs_numer"       ,"fraction of stubs per road (numerator);N_{stubs};f_{roads}"  ,15, 0.0, 15.0);
  
  TH1F     nPSstubs_denom     ("nPSstubs_denom"     ,"Number of PS stubs per road (denominator);N_{stubs};N_{roads}",15, 0.0, 15.0);
  TH1F     nPSstubs_numer     ("nPSstubs_numer"     ,"Number of PS stubs per road (numerator);N_{stubs};N_{roads}"  ,15, 0.0, 15.0);
  TH1F     fPSstubs_denom     ("fPSstubs_denom"     ,"Fraction of PS stubs per road (denominator);N_{stubs};f_{roads}",15, 0.0, 15.0);
  TH1F     fPSstubs_numer     ("fPSstubs_numer"     ,"fraction of PS stubs per road (numerator);N_{stubs};f_{roads}"  ,15, 0.0, 15.0);
  
  TProfile roads_per_tp       ("roads_per_tp"       , "Average number of roads per TP",2, 0.0, 2.0);
  
  TProfile n_roads_total      ("n_roads_total"      , "Total number of roads"         ,2, 0.0, 2.0);
  TProfile n_roads_good       ("n_roads_good"       , "Total number of good roads"    ,2, 0.0, 2.0);
  TProfile n_roads_nongood    ("n_roads_nongood"    , "Total number of non-good roads",2, 0.0, 2.0);
  
  TProfile n_roads_5lay ("n_roads_5lay" , "Total number of roads with stubs in 5 layers or more",2, 0.0, 2.0);
  
  
  
  // Run on the events in tree
  for (long long ievt = 0; ievt < neventsmax; ++ievt)
  {
    if (denomReader.loadTree(ievt) < 0) break;
    if (numerReader.loadTree(ievt) < 0) break;
    if (tpartReader.loadTree(ievt) < 0) break;

    denomReader.getEntry(ievt);
    numerReader.getEntry(ievt);
    tpartReader.getEntry(ievt);

    const size_t denomNRoads = denomReader.vr_patternRef->size();
    const size_t numerNRoads = numerReader.vr_patternRef->size();
    const size_t nTP         = tpartReader.vp2_pdgId    ->size();
    const size_t nStubs      = tpartReader.vb_tpId      ->size();

    if (verbosity >= 1 && ievt % 100 == 0)
      cout << Debug() << Form("... Processing event: %7lld", ievt) << endl;

    if (verbosity >= 3)
      cout << Debug() << "... evt: " << ievt << " # denom roads: " << denomNRoads << " # numer roads: " << numerNRoads << endl;
    
    
    
    
//     unordered_map<size_t, vector<size_t> > tpIdx_2_stubIdxes;
    unordered_map<size_t, unordered_set<size_t> > tpIdx_2_stubIdxes;
    
    for (size_t iStub = 0; iStub < nStubs; ++iStub)
    {
      int iTP = tpartReader.vb_tpId->at(iStub);
      
      if (iTP < 0)
        continue;
      
//       tpIdx_2_stubIdxes[ size_t(iTP) ].push_back( iStub );
      tpIdx_2_stubIdxes[ size_t(iTP) ].insert( iStub );
    }
    
    
    
    
//     unordered_map<size_t, vector<size_t> > stubIdx_2_denomRoadIdxes;
    unordered_map<size_t, unordered_set<size_t> > stubIdx_2_denomRoadIdxes;
    
    for (size_t iRoad = 0; iRoad < denomNRoads; ++iRoad)
    {
      unsigned const nLayers = (*(denomReader.vr_stubRefs))[iRoad].size();
      
      for (unsigned iLayer = 0; iLayer < nLayers; ++iLayer)
      {
        for (unsigned iStub  : (*(denomReader.vr_stubRefs))[iRoad][iLayer] )
        {
//           stubIdx_2_denomRoadIdxes[ iStub ].push_back( iRoad );
          stubIdx_2_denomRoadIdxes[ iStub ].insert( iRoad );
        }
      }
    }
    
    
//     unordered_map<size_t, vector<size_t> > stubIdx_2_numerRoadIdxes;
    unordered_map<size_t, unordered_set<size_t> > stubIdx_2_numerRoadIdxes;
    
    for (size_t iRoad = 0; iRoad < numerNRoads; ++iRoad)
    {
      unsigned const nLayers = (*(numerReader.vr_stubRefs))[iRoad].size();
      
      for (unsigned iLayer = 0; iLayer < nLayers; ++iLayer)
      {
        for (unsigned iStub  : (*(numerReader.vr_stubRefs))[iRoad][iLayer] )
        {
//           stubIdx_2_numerRoadIdxes[ iStub ].push_back( iRoad );
          stubIdx_2_numerRoadIdxes[ iStub ].insert( iRoad );
        }
      }
    }
    
    
    
    
//     unordered_map<size_t, vector<size_t> > tpIdx_2_denomRoadIdxes;
    unordered_map<size_t, unordered_set<size_t> > tpIdx_2_denomRoadIdxes;
    
    for (size_t iTP = 0; iTP < nTP; ++iTP)
    {
      if ( !tpartReader.vp2_intime->at(iTP) || !tpartReader.vp2_primary->at(iTP) )
        continue;
      
      auto& stubsIdxesForTP = tpIdx_2_stubIdxes[ iTP ];
      
      unordered_set<unsigned> roadIdxesForTP;
      
      for (unsigned iStub  : stubsIdxesForTP )
      {
        auto& roadIdxesForStub = stubIdx_2_denomRoadIdxes[ iStub ];
        
        for (size_t iRoad : roadIdxesForStub)
        {
          roadIdxesForTP.insert( iRoad );
        }
      }
      
      for (auto elem : roadIdxesForTP)
      {
//         tpIdx_2_denomRoadIdxes[ iTP ].push_back( elem );
        tpIdx_2_denomRoadIdxes[ iTP ].insert( elem );
      }
    }
    
    
//     unordered_map<size_t, vector<size_t> > tpIdx_2_numerRoadIdxes;
    unordered_map<size_t, unordered_set<size_t> > tpIdx_2_numerRoadIdxes;
    
    for (size_t iTP = 0; iTP < nTP; ++iTP)
    {
      if ( !tpartReader.vp2_intime->at(iTP) || !tpartReader.vp2_primary->at(iTP) )
        continue;
      
      auto& stubsIdxesForTP = tpIdx_2_stubIdxes[ iTP ];
      
      unordered_set<unsigned> roadIdxesForTP;
      
      for (unsigned iStub  : stubsIdxesForTP )
      {
        auto& roadIdxesForStub = stubIdx_2_numerRoadIdxes[ iStub ];
        
        for (size_t iRoad : roadIdxesForStub)
        {
          roadIdxesForTP.insert( iRoad );
        }
      }
      
      for (auto elem : roadIdxesForTP)
      {
//         tpIdx_2_numerRoadIdxes[ iTP ].push_back( elem );
        tpIdx_2_numerRoadIdxes[ iTP ].insert( elem );
      }
    }
    
    
    
    
    
    // Generate set of "good" TPs that will be used for the efficiency measurement
    
    unsigned const nMinimumLayersWithStubs = 5u;
    
    
    unordered_set<size_t> goodDenomTPIdxes;
    
    for (size_t iTP = 0; iTP < nTP; ++iTP)
    {
      // We are interested in primary particles and we ignore oot pile-up
      if ( !tpartReader.vp2_intime->at(iTP) && !tpartReader.vp2_primary->at(iTP) )
        continue;
      
      // Consider only pT > 3.0
      if ( tpartReader.vp2_pt->at(iTP) < 3.0)
        continue;
      
      // We require that at least one road has stubs in at least 5 layers
      bool satisfiesLayerRequirement = false;
      for (auto iRoad : tpIdx_2_denomRoadIdxes[iTP])
      {
        unsigned nLayersWithStubs = 0u;
        
        unsigned const nLayers = (*(denomReader.vr_stubRefs))[iRoad].size();
        
        for (unsigned iLayer = 0; iLayer < nLayers; ++iLayer)
        {
          bool foundStubFromTP = false;
          
          for (unsigned iStub  : (*(denomReader.vr_stubRefs))[iRoad][iLayer] )
          {
            if (tpIdx_2_stubIdxes[ iTP ].find(iStub) != tpIdx_2_stubIdxes[ iTP ].end())
            {
//               cout << "DENOM foundStubFromTP at iTP " << iTP << " iRoad " << iRoad << " iLayer " << iLayer << " iStub " << iStub << endl;
              foundStubFromTP = true;
              break;
            }
          }
          
          if (foundStubFromTP)
            ++nLayersWithStubs;
        }
        
        if (nLayersWithStubs >= nMinimumLayersWithStubs)
          satisfiesLayerRequirement = true;
      }
      
      if (!satisfiesLayerRequirement)
        continue;
      
      goodDenomTPIdxes.insert(iTP);
    }
    
    
    unordered_set<size_t> goodNumerTPIdxes;
    
    for (size_t iTP = 0; iTP < nTP; ++iTP)
    {
      // We are interested in primary particles and we ignore oot pile-up
      if ( !tpartReader.vp2_intime->at(iTP) && !tpartReader.vp2_primary->at(iTP) )
        continue;
      
      // Consider only pT > 3.0
      if ( tpartReader.vp2_pt->at(iTP) < 3.0)
        continue;
      
      // We require that at least one road has stubs in at least 5 layers
      bool satisfiesLayerRequirement = false;
      for (auto iRoad : tpIdx_2_numerRoadIdxes[iTP])
      {
        unsigned nLayersWithStubs = 0u;
        
        unsigned const nLayers = (*(numerReader.vr_stubRefs))[iRoad].size();
        
        for (unsigned iLayer = 0; iLayer < nLayers; ++iLayer)
        {
          bool foundStubFromTP = false;
          
          for (unsigned iStub  : (*(numerReader.vr_stubRefs))[iRoad][iLayer] )
          {
            if (tpIdx_2_stubIdxes[ iTP ].find(iStub) != tpIdx_2_stubIdxes[ iTP ].end())
            {
//               cout << "NUMER foundStubFromTP at iTP " << iTP << " iRoad " << iRoad << " iLayer " << iLayer << " iStub " << iStub << endl;
              foundStubFromTP = true;
              break;
            }
          }
          
          if (foundStubFromTP)
            ++nLayersWithStubs;
        }
        
        if (nLayersWithStubs >= nMinimumLayersWithStubs)
          satisfiesLayerRequirement = true;
      }
      
      if (!satisfiesLayerRequirement)
        continue;
      
      goodNumerTPIdxes.insert(iTP);
    }
    
//     cout << "Number of good roads in denominator sample: " << goodNumerTPIdxes.size() << endl;
//     cout << "Number of good roads in numerator sample  : " << goodDenomTPIdxes.size() << endl;
    
    
    
    
    
    
    // Efficiency for TPs
    for (auto iTP : goodDenomTPIdxes)
    {
      bool const pass = goodNumerTPIdxes.find(iTP) != goodNumerTPIdxes.end();
      
      tpeff_pt   .Fill(pass, tpartReader.vp2_pt ->at(iTP));
      tpeff_eta  .Fill(pass, tpartReader.vp2_eta->at(iTP));
    }
    
    
    // Efficiency for roads
    for (auto iTP : goodDenomTPIdxes)
//     for (size_t iTP = 0; iTP < nTP; ++iTP)
    {
      for (auto iRoad : tpIdx_2_denomRoadIdxes[iTP])
      {
        bool const pass = tpIdx_2_numerRoadIdxes[iTP].find(iRoad) != tpIdx_2_numerRoadIdxes[iTP].end();
        
        roadeff_pt .Fill(pass, tpartReader.vp2_pt ->at(iTP));
        roadeff_eta.Fill(pass, tpartReader.vp2_eta->at(iTP));
      }
    }
    
    
    // Number of roads per good TP
    for (auto iTP : goodDenomTPIdxes)
    {
      roads_per_tp.Fill("pre-filter", tpIdx_2_denomRoadIdxes[iTP].size());
    }
    
    for (auto iTP : goodNumerTPIdxes)
    {
      roads_per_tp.Fill("post-filter", tpIdx_2_numerRoadIdxes[iTP].size());
    }
    
    
    // Total number of roads
    unordered_set<size_t> goodDenomRoadsIdxes;
    for (auto iTP : goodDenomTPIdxes)
    {
      for (auto iRoad : tpIdx_2_denomRoadIdxes[iTP])
      {
        goodDenomRoadsIdxes.insert(iRoad);
      }
    }
    
    unordered_set<size_t> goodNumerRoadsIdxes;
    for (auto iTP : goodNumerTPIdxes)
    {
      for (auto iRoad : tpIdx_2_numerRoadIdxes[iTP])
      {
        goodNumerRoadsIdxes.insert(iRoad);
      }
    }
    
    
    {
      size_t const nRoadsTotDenom  = denomReader.vp2_charge->size();
      size_t const nRoadsTotNumer  = numerReader.vp2_charge->size();
      
      size_t const nGoodRoadsDenom = goodDenomRoadsIdxes.size();
      size_t const nGoodRoadsNumer = goodNumerRoadsIdxes.size();
      
      size_t const nBadRoadsDenom  = nRoadsTotDenom - nGoodRoadsDenom;
      size_t const nBadRoadsNumer  = nRoadsTotNumer - nGoodRoadsNumer;
      
      n_roads_total  .Fill("pre-filter" , nRoadsTotDenom );
      n_roads_total  .Fill("post-filter", nRoadsTotNumer );
      
      n_roads_good   .Fill("pre-filter" , nGoodRoadsDenom);
      n_roads_good   .Fill("post-filter", nGoodRoadsNumer);
      
      n_roads_nongood.Fill("pre-filter" , nBadRoadsDenom );
      n_roads_nongood.Fill("post-filter", nBadRoadsNumer );
    }
    
    {
      long nRoads5LayDenom  = 0.0;
      for (size_t iRoad = 0; iRoad < denomNRoads; ++iRoad)
      {
        unsigned nLayersWithStubs = 0u;
        
        unsigned const nLayers = (*(denomReader.vr_stubRefs))[iRoad].size();
        
        for (unsigned iLayer = 0; iLayer < nLayers; ++iLayer)
        {
          if (! (*(denomReader.vr_stubRefs))[iRoad][iLayer].empty() )
          {
            ++nLayersWithStubs;
          }
        }
        
        if (nLayersWithStubs >= 5)
          ++nRoads5LayDenom;
      }
      
      n_roads_5lay.Fill("pre-filter" , nRoads5LayDenom);
      
      
      long nRoads5LayNumer  = 0.0;
      for (size_t iRoad = 0; iRoad < numerNRoads; ++iRoad)
      {
        unsigned nLayersWithStubs = 0u;
        
        unsigned const nLayers = (*(numerReader.vr_stubRefs))[iRoad].size();
        
        for (unsigned iLayer = 0; iLayer < nLayers; ++iLayer)
        {
          if (! (*(numerReader.vr_stubRefs))[iRoad][iLayer].empty() )
          {
            ++nLayersWithStubs;
          }
        }
        
        if (nLayersWithStubs >= 5)
          ++nRoads5LayNumer;
      }
      
      n_roads_5lay.Fill("post-filter", nRoads5LayNumer);
    }
    
    
    
    
    
    // Number of all stubs
    for (size_t iRoad = 0; iRoad < denomNRoads; ++iRoad)
    {
      nstubs_denom.Fill( denomReader.vr_nstubs->at(iRoad) );
      fstubs_denom.Fill( denomReader.vr_nstubs->at(iRoad) );
    }
    
    for (size_t iRoad = 0; iRoad < numerNRoads; ++iRoad)
    {
      nstubs_numer.Fill( numerReader.vr_nstubs->at(iRoad) );
      fstubs_numer.Fill( numerReader.vr_nstubs->at(iRoad) );
    }
    
    
    // Number of PS stubs
    for (size_t iRoad = 0; iRoad < denomNRoads; ++iRoad)
    {
      unsigned const nLayers = (*(denomReader.vr_stubRefs))[iRoad].size();
      
      unsigned nPSStubs = 0u;
      for (unsigned iLayer = 0; iLayer < nLayers && iLayer < 3; ++iLayer)
      {
        nPSStubs += (*(denomReader.vr_stubRefs))[iRoad][iLayer].size();
      }
      
      nPSstubs_denom.Fill( nPSStubs );
      fPSstubs_denom.Fill( nPSStubs );
    }
    
    for (size_t iRoad = 0; iRoad < numerNRoads; ++iRoad)
    {
      unsigned const nLayers = (*(numerReader.vr_stubRefs))[iRoad].size();
      
      unsigned nPSStubs = 0u;
      for (unsigned iLayer = 0; iLayer < nLayers && iLayer < 3; ++iLayer)
      {
        nPSStubs += (*(numerReader.vr_stubRefs))[iRoad][iLayer].size();
      }
      
      nPSstubs_numer.Fill( nPSStubs );
      fPSstubs_numer.Fill( nPSStubs );
    }
    
  }
  
  
  fstubs_denom.Scale(1.0/fstubs_denom.Integral());
  fstubs_numer.Scale(1.0/fstubs_numer.Integral());
  
  
  fPSstubs_denom.Scale(1.0/fPSstubs_denom.Integral());
  fPSstubs_numer.Scale(1.0/fPSstubs_numer.Integral());
  
  
  
  
  TFile fout(outputfilename.c_str(),"RECREATE");
  fout.cd();
  
  tpeff_pt .Write();
  tpeff_eta.Write();
  
  roadeff_pt .Write();
  roadeff_eta.Write();
  
  nstubs_denom.Write();
  nstubs_numer.Write();
  fstubs_denom.Write();
  fstubs_numer.Write();
  
  nPSstubs_denom.Write();
  nPSstubs_numer.Write();
  fPSstubs_denom.Write();
  fPSstubs_numer.Write();
  
  roads_per_tp.Write();
  
  n_roads_total  .Write();
  n_roads_good   .Write();
  n_roads_nongood.Write();
  
  n_roads_5lay.Write();
  
  fout.Close();
  
  return 0;
}