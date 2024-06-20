
// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Common/DataModel/EventSelection.h"
//#include "Common/CCDB/EventSelectionParams.h"
//#include "Common/CCDB/TriggerAliases.h"
#include "CCDB/BasicCCDBManager.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"

#include <cmath>
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"

#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TObjArray.h"
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
//using namespace o2::analysis;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct ZDCAnalysis {
  
  //Filters
  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtMin) && (aod::track::pt < cfgCutPtMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls);
  
  
  //define output
  HistogramRegistry registry;
  
  //define my.....
  using myTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>>;
  
  using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
  
  //define global parameters
  double collisionCounter;
  
  int counter=0;
  
  void init(InitContext const&)
  {
    

    //Tower mean energies vs. centrality
    registry.add("hZNA_mean_t1_cent","", {HistType::kTProfile, {{axisCent}}});
    registry.add("hZNA_mean_t2_cent","", {HistType::kTProfile, {{axisCent}}});
    registry.add("hZNA_mean_t3_cent","", {HistType::kTProfile, {{axisCent}}});
    registry.add("hZNA_mean_t4_cent","", {HistType::kTProfile, {{axisCent}}});
    registry.add("hZNA_mean_common_cent","", {HistType::kTProfile, {{axisCent}}});
    
    registry.add("hZNC_mean_t1_cent","", {HistType::kTProfile, {{axisCent}}});
    registry.add("hZNC_mean_t2_cent","", {HistType::kTProfile, {{axisCent}}});
    registry.add("hZNC_mean_t3_cent","", {HistType::kTProfile, {{axisCent}}});
    registry.add("hZNC_mean_t4_cent","", {HistType::kTProfile, {{axisCent}}});
    registry.add("hZNC_mean_common_cent","", {HistType::kTProfile, {{axisCent}}});
    
  }
  
  void processRun2(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>>::iterator const& collision, aod::Zdcs const&, myTracks const& tracks) //

  {
      int Ntot = tracks.size();
      if (Ntot < 1)
        return;
      
      if (!collision.sel7())
        return;
      auto cent = collision.centRun2V0M();
      if(cent<0 || cent>100)
        return;
      
      if (collision.foundZDCId() >= 0) { //
        
        registry.fill(HIST("hZNA_mean_common_cent"),int(cent),collision.foundZDC().energyCommonZNA());
        registry.fill(HIST("hZNA_mean_t1_cent"),int(cent),collision.foundZDC().energySectorZNA()[0]);
        registry.fill(HIST("hZNA_mean_t2_cent"),int(cent),collision.foundZDC().energySectorZNA()[1]);
        registry.fill(HIST("hZNA_mean_t3_cent"),int(cent),collision.foundZDC().energySectorZNA()[2]);
        registry.fill(HIST("hZNA_mean_t4_cent"),int(cent),collision.foundZDC().energySectorZNA()[3]);
        
        registry.fill(HIST("hZNC_mean_common_cent"),int(cent),collision.foundZDC().energyCommonZNC());
        registry.fill(HIST("hZNC_mean_t1_cent"),int(cent),collision.foundZDC().energySectorZNC()[0]);
        registry.fill(HIST("hZNC_mean_t2_cent"),int(cent),collision.foundZDC().energySectorZNC()[1]);
        registry.fill(HIST("hZNC_mean_t3_cent"),int(cent),collision.foundZDC().energySectorZNC()[2]);
        registry.fill(HIST("hZNC_mean_t4_cent"),int(cent),collision.foundZDC().energySectorZNC()[3]);
 
      }
  }
  PROCESS_SWITCH(ZDCAnalysis, processRun2, "Processing ZDC mean energy for run 2", false);
  
  void processRun3(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>>::iterator const& collision, aod::BCsWithTimestamps const&, BCsRun3 const& bcs, aod::Zdcs const&, myTracks const& tracks) //

  {
      int Ntot = tracks.size();
      if (Ntot < 1)
        return;
      
      if (!collision.sel8())
        return;
      auto cent = collision.centFT0Cs();
      if(cent<0 || cent>90)
        return;
      
    const auto& foundBC = collision.foundBC_as<BCsRun3>();
    if(foundBC.has_zdc()){

      zdcCol = foundBC.zdc();
      
        registry.fill(HIST("hZNA_mean_common_cent"),int(cent),zdcCol.energyCommonZNA());
        registry.fill(HIST("hZNA_mean_t1_cent"),int(cent),zdcCol.energySectorZNA()[0]);
        registry.fill(HIST("hZNA_mean_t2_cent"),int(cent),zdcCol.energySectorZNA()[1]);
        registry.fill(HIST("hZNA_mean_t3_cent"),int(cent),zdcCol.energySectorZNA()[2]);
        registry.fill(HIST("hZNA_mean_t4_cent"),int(cent),zdcCol.energySectorZNA()[3]);
        
        registry.fill(HIST("hZNC_mean_common_cent"),int(cent),zdcCol.energyCommonZNC());
        registry.fill(HIST("hZNC_mean_t1_cent"),int(cent),zdcCol.energySectorZNC()[0]);
        registry.fill(HIST("hZNC_mean_t2_cent"),int(cent),zdcCol.energySectorZNC()[1]);
        registry.fill(HIST("hZNC_mean_t3_cent"),int(cent),zdcCol.energySectorZNC()[2]);
        registry.fill(HIST("hZNC_mean_t4_cent"),int(cent),zdcCol.energySectorZNC()[3]);
 
      } // end of has_zdc()
  }
  PROCESS_SWITCH(ZDCAnalysis, processRun3, "Processing ZDC mean energy for run 3", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ZDCAnalysis>(cfgc)};

  
}
