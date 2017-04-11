#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/HTRZAlgorithmRecorder.h"
#include "SLHCL1TrackTriggerSimulations/AMSimulation/interface/HTRZAlgorithm.h"

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TDirectoryFile.h"
#include "TH1.h"
#include "TH2.h"

#include <iostream>

using namespace slhcl1tt;

HTRZAlgorithmRecorder::HTRZAlgorithmRecorder():
  htmatrix_event_number_         (0),
  htmatrix_road_number_          (0),
  htmatrix_tower_number_         (0),
  htmatrix_pattern_number_       (0),
  tree_file_                     (nullptr),
  tree_                          (nullptr),
  htmatrix_layercount_allstubs_  (nullptr),
  htmatrix_layercount_psstubs_   (nullptr),
  htmatrix_layercount_2sstubs_   (nullptr),
  htmatrix_activated_            (nullptr)
{
  
}


HTRZAlgorithmRecorder::HTRZAlgorithmRecorder(const char* const out_filename, HTRZAlgorithmConfig const& monitored_algo_config)
{
  htmatrix_layercount_allstubs_ = new TH2C("htmatrix_layercount_allstubs" ,"Number of all layers with stubs in HT matrix cells for a single road",
    monitored_algo_config.nbins_cotantheta,
    monitored_algo_config.min_cotantheta  ,
    monitored_algo_config.max_cotantheta  ,
    monitored_algo_config.nbins_zT        ,
    monitored_algo_config.min_zT          ,
    monitored_algo_config.max_zT          
  );
  
  htmatrix_layercount_psstubs_  = new TH2C("htmatrix_layercount_psstubs"  ,"Number of PS layers with stubs in HT matrix cells for a single road",
    monitored_algo_config.nbins_cotantheta,
    monitored_algo_config.min_cotantheta  ,
    monitored_algo_config.max_cotantheta  ,
    monitored_algo_config.nbins_zT        ,
    monitored_algo_config.min_zT          ,
    monitored_algo_config.max_zT          
  );
  
  htmatrix_layercount_2sstubs_  = new TH2C("htmatrix_layercount_2sstubs"  ,"Number of 2S layers with stubs in HT matrix cells for a single road",
    monitored_algo_config.nbins_cotantheta,
    monitored_algo_config.min_cotantheta  ,
    monitored_algo_config.max_cotantheta  ,
    monitored_algo_config.nbins_zT        ,
    monitored_algo_config.min_zT          ,
    monitored_algo_config.max_zT          
  );
  
  
  htmatrix_activated_           = new TH2C("htmatrix_layercount_activated","Activated cells in HT matrix cells for a single road",
    monitored_algo_config.nbins_cotantheta,
    monitored_algo_config.min_cotantheta  ,
    monitored_algo_config.max_cotantheta  ,
    monitored_algo_config.nbins_zT        ,
    monitored_algo_config.min_zT          ,
    monitored_algo_config.max_zT          
  );
  
  
  
  
  htmatrix_layercount_allstubs_  ->SetDirectory(nullptr);
  htmatrix_layercount_psstubs_   ->SetDirectory(nullptr);
  htmatrix_layercount_2sstubs_   ->SetDirectory(nullptr);
  htmatrix_activated_            ->SetDirectory(nullptr);
  
  tree_file_ = new TFile(out_filename,"RECREATE");
  tree_file_->cd();
  
  
  htmatrix_sums_layercount_allstubs_ = new TH2D("htmatrix_sums_layercount_allstubs" ,"Number of all layers with stubs in HT matrix cells for all processed roads",
    monitored_algo_config.nbins_cotantheta      ,
                                            -0.5,
    monitored_algo_config.nbins_cotantheta - 0.5,
    monitored_algo_config.nbins_zT              ,
                                            -0.5,
    monitored_algo_config.nbins_zT         - 0.5
  );
  
  htmatrix_sums_layercount_psstubs_  = new TH2D("htmatrix_sums_layercount_psstubs"  ,"Number of PS layers with stubs in HT matrix cells for all processed roads",
    monitored_algo_config.nbins_cotantheta      ,
                                            -0.5,
    monitored_algo_config.nbins_cotantheta - 0.5,
    monitored_algo_config.nbins_zT              ,
                                            -0.5,
    monitored_algo_config.nbins_zT         - 0.5
  );
  
  htmatrix_sums_layercount_2sstubs_  = new TH2D("htmatrix_sums_layercount_2sstubs"  ,"Number of 2S layers with stubs in HT matrix cells for all processed roads",
    monitored_algo_config.nbins_cotantheta      ,
                                            -0.5,
    monitored_algo_config.nbins_cotantheta - 0.5,
    monitored_algo_config.nbins_zT              ,
                                            -0.5,
    monitored_algo_config.nbins_zT         - 0.5
  );
  
  
  htmatrix_sums_activated_           = new TH2D("htmatrix_sums_layercount_activated","Activated cells in HT matrix cells for all processed roads",
    monitored_algo_config.nbins_cotantheta      ,
                                            -0.5,
    monitored_algo_config.nbins_cotantheta - 0.5,
    monitored_algo_config.nbins_zT              ,
                                            -0.5,
    monitored_algo_config.nbins_zT         - 0.5
  );
  
  
  tree_ = new TTree("playback_tree","tree used to store the playback of HT matrix");
  
  tree_->Branch("htmatrix_event_number/i"   , &htmatrix_event_number_  );
  tree_->Branch("htmatrix_road_number/i"    , &htmatrix_road_number_   );
  tree_->Branch("htmatrix_tower_number/i"   , &htmatrix_tower_number_  );
  tree_->Branch("htmatrix_pattern_number/i" , &htmatrix_pattern_number_);
  
  tree_->Branch("htmatrix_layercount_allstubs" ,"TH2C", htmatrix_layercount_allstubs_ );
  tree_->Branch("htmatrix_layercount_psstubs"  ,"TH2C", htmatrix_layercount_psstubs_  );
  tree_->Branch("htmatrix_layercount_2sstubs"  ,"TH2C", htmatrix_layercount_2sstubs_  );
  tree_->Branch("htmatrix_layercount_activated","TH2C", htmatrix_activated_           );
}

HTRZAlgorithmRecorder::~HTRZAlgorithmRecorder()
{
  tree_file_->cd();
  tree_     ->Write(tree_->GetName(), TObject::kOverwrite);
  
  htmatrix_sums_layercount_allstubs_->Write(htmatrix_sums_layercount_allstubs_->GetName(), TObject::kOverwrite);
  htmatrix_sums_layercount_psstubs_ ->Write(htmatrix_sums_layercount_psstubs_ ->GetName(), TObject::kOverwrite);
  htmatrix_sums_layercount_2sstubs_ ->Write(htmatrix_sums_layercount_2sstubs_ ->GetName(), TObject::kOverwrite);
  htmatrix_sums_activated_          ->Write(htmatrix_sums_activated_          ->GetName(), TObject::kOverwrite);
  
  tree_file_->Flush();
  tree_file_->Close();
  
  // Apparently, the TTree is destroyed upon file closure
//   if (tree_                         ) {delete tree_                          ; tree_                          = nullptr;}
  if (tree_file_                    ) {delete tree_file_                     ; tree_file_                     = nullptr;}
  if (htmatrix_layercount_allstubs_ ) {delete htmatrix_layercount_allstubs_  ; htmatrix_layercount_allstubs_  = nullptr;}
  if (htmatrix_layercount_psstubs_  ) {delete htmatrix_layercount_psstubs_   ; htmatrix_layercount_psstubs_   = nullptr;}
  if (htmatrix_layercount_2sstubs_  ) {delete htmatrix_layercount_2sstubs_   ; htmatrix_layercount_2sstubs_   = nullptr;}
  if (htmatrix_activated_           ) {delete htmatrix_activated_            ; htmatrix_activated_            = nullptr;}
}


void HTRZAlgorithmRecorder::Snapshot(HTRZAlgorithm& monitored_algo, unsigned const event_number, unsigned const tower_number, unsigned const road_number, unsigned const pattern_number)
{
  assert(htmatrix_layercount_allstubs_->GetXaxis()->GetNbins() == int(monitored_algo.config_.nbins_cotantheta));
  assert(htmatrix_layercount_psstubs_ ->GetXaxis()->GetNbins() == int(monitored_algo.config_.nbins_cotantheta));
  assert(htmatrix_layercount_2sstubs_ ->GetXaxis()->GetNbins() == int(monitored_algo.config_.nbins_cotantheta));
  assert(htmatrix_activated_          ->GetXaxis()->GetNbins() == int(monitored_algo.config_.nbins_cotantheta));
  
  assert(htmatrix_layercount_allstubs_->GetYaxis()->GetNbins() == int(monitored_algo.config_.nbins_zT));
  assert(htmatrix_layercount_psstubs_ ->GetYaxis()->GetNbins() == int(monitored_algo.config_.nbins_zT));
  assert(htmatrix_layercount_2sstubs_ ->GetYaxis()->GetNbins() == int(monitored_algo.config_.nbins_zT));
  assert(htmatrix_activated_          ->GetYaxis()->GetNbins() == int(monitored_algo.config_.nbins_zT));
  
  
  htmatrix_event_number_   = event_number  ;
  htmatrix_road_number_    = road_number   ;
  htmatrix_tower_number_   = tower_number  ;
  htmatrix_pattern_number_ = pattern_number;
  
  htmatrix_layercount_allstubs_ ->Reset("ICES");
  htmatrix_layercount_psstubs_  ->Reset("ICES");
  htmatrix_layercount_2sstubs_  ->Reset("ICES");
  htmatrix_activated_           ->Reset("ICES");
  htmatrix_layercount_allstubs_ ->SetEntries(1);
  htmatrix_layercount_psstubs_  ->SetEntries(1);
  htmatrix_layercount_2sstubs_  ->SetEntries(1);
  htmatrix_activated_           ->SetEntries(1);
  
  
  for (unsigned iCot = 0; iCot < monitored_algo.config_.nbins_cotantheta; ++iCot)
    for (unsigned iZT = 0; iZT < monitored_algo.config_.nbins_zT; ++iZT)
    {
      unsigned const n_stubs_all = monitored_algo.count_stubs_all_layers(iCot,iZT);
      unsigned const n_stubs_ps  = monitored_algo.count_stubs_ps_layers (iCot,iZT);
      
      
      htmatrix_layercount_allstubs_->SetBinContent(iCot, iZT, n_stubs_all             );
      htmatrix_layercount_psstubs_ ->SetBinContent(iCot, iZT, n_stubs_ps              );
      htmatrix_layercount_2sstubs_ ->SetBinContent(iCot, iZT, n_stubs_all - n_stubs_ps);
      htmatrix_activated_          ->SetBinContent(iCot, iZT, monitored_algo.majority(iCot,iZT));
      
      htmatrix_sums_layercount_allstubs_->Fill(iCot, iZT, n_stubs_all             );
      htmatrix_sums_layercount_psstubs_ ->Fill(iCot, iZT, n_stubs_ps              );
      htmatrix_sums_layercount_2sstubs_ ->Fill(iCot, iZT, n_stubs_all - n_stubs_ps);
      htmatrix_sums_activated_          ->Fill(iCot, iZT, monitored_algo.majority(iCot,iZT));
    }
  
//   tree_->Fill();
}

