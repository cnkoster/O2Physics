
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

namespace o2::analysis::test{
std::shared_ptr<TProfile2D> histograms[10] = {{nullptr}};
}
using namespace o2::analysis::test;
#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct ZDCAnalysis {
  
  ConfigurableAxis axisCent{"axisCent", {90,0,90},"Centrality axis in 1% bins"};
  
  O2_DEFINE_CONFIGURABLE(cfgCutVertex,        float,      10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin,         float,      0.2f, "Minimal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax,         float,      3.0f, "Maximal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta,           float,      0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls,  float,      2.5, "Chi2 per TPC clusters")
  
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
    //Tower mean energies vs. centrality ( common=t0 )
    for(int tower = 0; tower<5; tower++){
      histograms[tower] = registry.add<TProfile2D>(Form("hZNA_mean_t%i_cent", tower),"", kTProfile2D, {{1,0,1},axisCent});
      histograms[tower+5] = registry.add<TProfile2D>(Form("hZNC_mean_t%i_cent", tower),"", kTProfile2D, {{1,0,1},axisCent});
    }
    
  }
  
  void processRun2(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>>::iterator const& collision, 
                   aod::Zdcs const&,
                   myTracks const& tracks)
  {
    int Ntot = tracks.size();
    if (Ntot < 1)
      return;
    
    if (!collision.sel7())
      return;
    auto cent = collision.centRun2V0M();
    if(cent<0 || cent>90)
      return;
    
    int runnumber = collision.bc().runNumber();
    
    if (collision.foundZDCId() >= 0) { //
    
      for(int tower = 0; tower<5; tower++){
        if(tower==0){
          histograms[tower]->Fill(Form("%d", runnumber), int(cent),collision.foundZDC().energyCommonZNA(),1);
          histograms[tower+5]->Fill(Form("%d", runnumber), int(cent),collision.foundZDC().energyCommonZNA(),1);
        } else {
          histograms[tower]->Fill(Form("%d", runnumber), int(cent),collision.foundZDC().energySectorZNA()[tower-1],1);
          histograms[tower+5]->Fill(Form("%d", runnumber), int(cent),collision.foundZDC().energySectorZNC()[tower-1],1);
        }
      }
    }
  }
  PROCESS_SWITCH(ZDCAnalysis, processRun2, "Processing ZDC mean energy for run 2", false);
  
  void processRun3(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>>::iterator const& collision,
                   BCsRun3 const& bcs,
                   aod::Zdcs const&,
                   myTracks const& tracks)
  {
    int Ntot = tracks.size();
    if (Ntot < 1)
      return;
    
    if (!collision.sel8())
      return;
    auto cent = collision.centFT0C();
    if(cent<0 || cent>90)
      return;
    
    const auto& foundBC = collision.foundBC_as<BCsRun3>();
   
    if(foundBC.has_zdc()){
      
      int runnumber = foundBC.runNumber();
      auto zdcCol = foundBC.zdc();
    
      for(int tower = 0; tower<5; tower++){
        if(tower==0){
          if(zdcCol.energyCommonZNA()>0)histograms[tower]->Fill(Form("%d", runnumber), int(cent), zdcCol.energyCommonZNA(),1);
          if(zdcCol.energyCommonZNA()>0)histograms[tower+5]->Fill(Form("%d", runnumber), int(cent), zdcCol.energyCommonZNA(),1);
        } else {
          if(zdcCol.energySectorZNA()[tower-1]>0)histograms[tower]->Fill(Form("%d", runnumber), int(cent), zdcCol.energySectorZNA()[tower-1],1);
          if(zdcCol.energySectorZNC()[tower-1]>0)histograms[tower+5]->Fill(Form("%d", runnumber), int(cent), zdcCol.energySectorZNC()[tower-1],1);
        }
      }
      
    } // end of has_zdc()
  }
  PROCESS_SWITCH(ZDCAnalysis, processRun3, "Processing ZDC mean energy for run 3", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ZDCAnalysis>(cfgc)};
  
  
}
