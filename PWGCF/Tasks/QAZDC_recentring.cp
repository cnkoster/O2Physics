
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
#include "TProfile.h"
#include "TObjArray.h"
#include <vector>
#include <stdlib.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
//using namespace o2::analysis;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};



struct ZDCAnalysis {
  
  O2_DEFINE_CONFIGURABLE(cfgCutVertex,        float,      10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMin,      float,      0.2f, "Minimal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMax,      float,      10.0f, "Maximal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin,         float,      0.2f, "Minimal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax,         float,      3.0f, "Maximal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta,           float,      0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls,  float,      2.5, "Chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgUseNch,           bool,       false, "Use Nch for flow observables")
  O2_DEFINE_CONFIGURABLE(cfgNbootstrap,       int,        10, "Number of subsamples")
  O2_DEFINE_CONFIGURABLE(nBins1,              int,        100, "NBins1")
  O2_DEFINE_CONFIGURABLE(nBins2,              int,        100, "NBins2")
  O2_DEFINE_CONFIGURABLE(nBins3,              int,        100, "NBins3")
  O2_DEFINE_CONFIGURABLE(MaxZN,               float,      150000, "MaxZN")
  O2_DEFINE_CONFIGURABLE(MaxMultNTracks,      float,      1000 ,"MaxMultNTracks")
  O2_DEFINE_CONFIGURABLE(centmin,             double,     0, "Minimal centrality")
  O2_DEFINE_CONFIGURABLE(centmax,             double,     100, "Maximal centrality")
  O2_DEFINE_CONFIGURABLE(fillHistograms,      bool,       true, "Fill Histogram registry")
  O2_DEFINE_CONFIGURABLE(dirMean,             std::string, "/Users/noorkoster/alice/Analysis/ZDC/AnalysisResultsCalib.root", "path to files with mean histos (RBR)")
//  O2_DEFINE_CONFIGURABLE(fInputFileName, TString, "", "Input file for calibration histos")
  
  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisPhi{"axisPhi", {60, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis axisEta{"axisEta", {40, -1., 1.}, "eta axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.2, 0.25, 0.30, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00}, "pt axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90}, "centrality axis for histograms"};
  ConfigurableAxis axisZDCSector{"axisZDCSector",{8,0,8}, "ZDC sectors"};
  ConfigurableAxis axisEnergy{"axisEnergy",{10000,0,1e5},"Energy deposititon in ZDC sector"};
  ConfigurableAxis axisEnergyCom{"axisEnergyCom",{40000,0,400000},"Energy deposititon for common ZDC sector"};
  ConfigurableAxis axisTmin{"axisTmin", {100,-30,30}, "Zdc timing A-C"};
  ConfigurableAxis axisTplus{"axisTplus", {1000,0,1000},"Zdc timing A+C"};
  ConfigurableAxis axisQ{"axisQ", {100,-2,2},"Q vector (xy) in ZDC"};
  ConfigurableAxis axisPosCollisionY{"axisPosCollisionY", {2000,-0.3,0.45},"for Pos X Y of collision"};
  ConfigurableAxis axisPosCollisionX{"axisPosCollisionX", {2000,-0.06,0.1},"for Pos X Y of collision"};
  ConfigurableAxis axisPosCollisionZ{"axisPosCollisionZ", {2000,10,10},"for vz of collision"};
  ConfigurableAxis axisCent{"axisCent", {100,0,100},"Centrality axis in 1% bins"};
  
  //Filters
  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtMin) && (aod::track::pt < cfgCutPtMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls);
  
  
  //define output
  HistogramRegistry registry{"Registry"};
  //  HistogramRegistry registryMean{"RegistryMean"};
  
  //define my.....
  using myCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>>;
  // snap niet precies wat Run2MatchedSparse doet.. maar laten we het testen (ze gebruiken het in de tut)
  using myTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>>;
  
  //define global parameters
  double collisionCounter;
  float sumQxA; float sumQxC; float sumQyA; float sumQyC;
  
  //recentered Q-vectors
  float recQxA; float recQyA; float recQxC; float recQyC;
  
  int counter=0;
  //number of towers on a and c side. (4 per side)
  int nTowers=4;
  
  TFile *f = TFile::Open(dirMean->c_str());
  
  int bin=0;
  int lastRunNumber = -1;
  int numFoundint = -1;
  
  void init(InitContext const&)
  {
    if(doprocessCalib){
      //Tower mean energies vs. centrality
      registry.add("hZNA_mean_t1_cent","", HistType::kTProfile2D, {{1,0.,1.} ,axisCent });
      registry.add("hZNA_mean_t2_cent","", HistType::kTProfile2D, {{1,0.,1.} ,axisCent });
      registry.add("hZNA_mean_t3_cent","", HistType::kTProfile2D, {{1,0.,1.} ,axisCent });
      registry.add("hZNA_mean_t4_cent","", HistType::kTProfile2D, {{1,0.,1.} ,axisCent });
      registry.add("hZNA_mean_common_cent","", HistType::kTProfile2D, {{1,0.,1.} ,axisCent });
      
      registry.add("hZNC_mean_t1_cent","", HistType::kTProfile2D, {{1,0.,1.} ,axisCent });
      registry.add("hZNC_mean_t2_cent","", HistType::kTProfile2D, {{1,0.,1.} ,axisCent });
      registry.add("hZNC_mean_t3_cent","", HistType::kTProfile2D, {{1,0.,1.} ,axisCent });
      registry.add("hZNC_mean_t4_cent","", HistType::kTProfile2D, {{1,0.,1.} ,axisCent });
      registry.add("hZNC_mean_common_cent","", HistType::kTProfile2D, {{1,0.,1.} ,axisCent });
    }
    
    if(doprocessData){
           
      //ZDC
      registry.add("ZNAcoll", "ZNAcoll", {HistType::kTH1F, {{nBins1, -10., MaxZN}}});
      registry.add("ZNCcoll", "ZNCcoll", {HistType::kTH1F, {{nBins1, -10., MaxZN}}});
      registry.add("ZNAvsZNCcoll", "ZNAvsZNCcoll", {HistType::kTH2F, {{{nBins1, -10., MaxZN}, {nBins1, -10., MaxZN}}}});
      registry.add("ZNCadcvstdccoll", "ZNCadcvstdccoll", {HistType::kTH2F, {{{400, -50., 50.}, {nBins1, -10., MaxZN}}}});
      registry.add("ZNAadcvstdccoll", "ZNAadcvstdccoll", {HistType::kTH2F, {{{400, -50., 50.}, {nBins1, -10., MaxZN}}}});
      
      registry.add("hAllTowersEnergy_processMean","", {HistType::kTH2F, {axisZDCSector,axisEnergy}});
      registry.add("hAllTowersEnergy","", {HistType::kTH2F, {axisZDCSector,axisEnergy}});
      registry.add("hAllTowersEnergy_noEq","", {HistType::kTH2F, {axisZDCSector,axisEnergy}});
      
      registry.add("hZNA_s1","", {HistType::kTH1F, {axisEnergy}});
      registry.add("hZNA_s2","", {HistType::kTH1F, {axisEnergy}});
      registry.add("hZNA_s3","", {HistType::kTH1F, {axisEnergy}});
      registry.add("hZNA_s4","", {HistType::kTH1F, {axisEnergy}});
      registry.add("hZNAC_common","", {HistType::kTH2F, {{2,0,2}, axisEnergy}});
      registry.add("hZNA_common","", {HistType::kTH2F, {axisEnergy, axisEnergy}});
      
      registry.add("hZNC_s1","", {HistType::kTH1F, {axisEnergy}});
      registry.add("hZNC_s2","", {HistType::kTH1F, {axisEnergy}});
      registry.add("hZNC_s3","", {HistType::kTH1F, {axisEnergy}});
      registry.add("hZNC_s4","", {HistType::kTH1F, {axisEnergy}});
      registry.add("hZNC_common","", {HistType::kTH2F, {axisEnergy, axisEnergy}});
      
      registry.add("hZNA_Qx","hZNA_Qx", {HistType::kTH1F, {axisQ}});
      registry.add("hZNA_Qy","hZNA_Qy", {HistType::kTH1F, {axisQ}});
      
      registry.add("hZNC_Qx","hZNC_Qx", {HistType::kTH1F, {axisQ}});
      registry.add("hZNC_Qy","hZNC_Qy", {HistType::kTH1F, {axisQ}});
      
      registry.add("hZNCvsZNA_Qx","", {HistType::kTH2F, {axisQ, axisQ}});
      registry.add("hZNCvsZNA_Qy","", {HistType::kTH2F, {axisQ, axisQ}});
      
      //Tower energy vs. cent
      registry.add("hZNA_t1_cent","", {HistType::kTH2F, {{axisCent,axisEnergy}}});
      registry.add("hZNA_t2_cent","", {HistType::kTH2F, {{axisCent,axisEnergy}}});
      registry.add("hZNA_t3_cent","", {HistType::kTH2F, {{axisCent,axisEnergy}}});
      registry.add("hZNA_t4_cent","", {HistType::kTH2F, {{axisCent,axisEnergy}}});
      registry.add("hZNA_common_cent","", {HistType::kTH2F, {{axisCent,axisEnergyCom}}});
      
      registry.add("hZNC_t1_cent","", {HistType::kTH2F, {{axisCent,axisEnergy}}});
      registry.add("hZNC_t2_cent","", {HistType::kTH2F, {{axisCent,axisEnergy}}});
      registry.add("hZNC_t3_cent","", {HistType::kTH2F, {{axisCent,axisEnergy}}});
      registry.add("hZNC_t4_cent","", {HistType::kTH2F, {{axisCent,axisEnergy}}});
      registry.add("hZNC_common_cent","", {HistType::kTH2F, {{axisCent,axisEnergyCom}}});
      
      //Calibrated energy vs centrality
      registry.add("hZNA_t1_cent_cal","", {HistType::kTH2F, {{axisCent,axisEnergy}}});
      registry.add("hZNA_t2_cent_cal","", {HistType::kTH2F, {{axisCent,axisEnergy}}});
      registry.add("hZNA_t3_cent_cal","", {HistType::kTH2F, {{axisCent,axisEnergy}}});
      registry.add("hZNA_t4_cent_cal","", {HistType::kTH2F, {{axisCent,axisEnergy}}});
      
      registry.add("hZNC_t1_cent_cal","", {HistType::kTH2F, {{axisCent,axisEnergy}}});
      registry.add("hZNC_t2_cent_cal","", {HistType::kTH2F, {{axisCent,axisEnergy}}});
      registry.add("hZNC_t3_cent_cal","", {HistType::kTH2F, {{axisCent,axisEnergy}}});
      registry.add("hZNC_t4_cent_cal","", {HistType::kTH2F, {{axisCent,axisEnergy}}});
      
      
      //Q-vector dependence
      registry.add("hZNA_Qx_vs_Centrality","hZNA_Qx_vs_Centrality", {HistType::kTH2F, {axisQ, axisMultiplicity}});
      registry.add("hZNA_Qy_vs_Centrality","hZNA_Qy_vs_Centrality", {HistType::kTH2F, {axisQ, axisMultiplicity}});
      
      registry.add("hZNC_Qx_vs_Centrality","hZNC_Qx_vs_Centrality", {HistType::kTH2F, {axisQ, axisMultiplicity}});
      registry.add("hZNC_Qy_vs_Centrality","hZNC_Qy_vs_Centrality", {HistType::kTH2F, {axisQ, axisMultiplicity}});
      
      registry.add("hZNA_Qx_vs_vx","hZNA_Qx_vs_vx", {HistType::kTH2F, {axisQ, axisPosCollisionX}});
      registry.add("hZNA_Qy_vs_vx","hZNA_Qy_vs_vx", {HistType::kTH2F, {axisQ, axisPosCollisionX}});
      registry.add("hZNA_Qx_vs_vy","hZNA_Qx_vs_vy", {HistType::kTH2F, {axisQ, axisPosCollisionY}});
      registry.add("hZNA_Qy_vs_vy","hZNA_Qy_vs_vy", {HistType::kTH2F, {axisQ, axisPosCollisionY}});
      registry.add("hZNA_Qx_vs_vz","hZNA_Qx_vs_vz", {HistType::kTH2F, {axisQ, axisVertex}});
      registry.add("hZNA_Qy_vs_vz","hZNA_Qy_vs_vz", {HistType::kTH2F, {axisQ, axisVertex}});
      
      registry.add("hZNC_Qx_vs_vx","hZNC_Qx_vs_vx", {HistType::kTH2F, {axisQ, axisPosCollisionX}});
      registry.add("hZNC_Qy_vs_vx","hZNC_Qy_vs_vx", {HistType::kTH2F, {axisQ, axisPosCollisionX}});
      registry.add("hZNC_Qx_vs_vy","hZNC_Qx_vs_vy", {HistType::kTH2F, {axisQ, axisPosCollisionY}});
      registry.add("hZNC_Qy_vs_vy","hZNC_Qy_vs_vy", {HistType::kTH2F, {axisQ, axisPosCollisionY}});
      registry.add("hZNC_Qx_vs_vz","hZNC_Qx_vs_vz", {HistType::kTH2F, {axisQ, axisPosCollisionZ}});
      registry.add("hZNC_Qy_vs_vz","hZNC_Qy_vs_vz", {HistType::kTH2F, {axisQ, axisPosCollisionZ}});
    }
    
    //tracks
    registry.add("hPhi", "", {HistType::kTH1F, {axisPhi}});
    registry.add("hEta", "", {HistType::kTH1F, {axisEta}});
    registry.add("hPt", "", {HistType::kTH1F, {axisPt}});
    
    //collisions
    registry.add("hVtxZ", "", {HistType::kTH1F, {axisVertex}});
    registry.add("hVtxY", "", {HistType::kTH1F, {axisPosCollisionY}});
    registry.add("hVtxX", "", {HistType::kTH1F, {axisPosCollisionX}});
    registry.add("hMult", "", {HistType::kTH1F, {{3000, 0.5, 3000.5}}});
    registry.add("hCent", "", {HistType::kTH1F, {{90, 0, 90}}});
    


    
  }
  
  void processCalib(myCollisions::iterator const& collision, aod::Zdcs const& zdcs, aod::BCs const& bcs, myTracks const& tracks)
  //Belangrijk dus dat je wel subscribed hier naar aod::Zdcs om toegang te krijgen tot zdc informatie
  {
    
    int Ntot = tracks.size();
    if (Ntot < 1)
      return;
    if (!collision.sel7())
      return;
    auto cent = collision.centRun2V0M();
    if(cent<0 || cent>100)
      return;
    
    if (collision.foundZDCId() >= 0) {
      
      // check if we have the right gaineq files
      
      int runNumber = collision.bc().runNumber();
      //      LOGF(info, "RUNNUMBER = %i", runNumber);
      

      
      // if calib returns false, there is no calibration data for the given run
      float EcomZNA = collision.foundZDC().energyCommonZNA();
      float EcomZNC = collision.foundZDC().energyCommonZNC();
      float E1ZNA = collision.foundZDC().energySectorZNA()[0];
      float E2ZNA = collision.foundZDC().energySectorZNA()[1];
      float E3ZNA = collision.foundZDC().energySectorZNA()[2];
      float E4ZNA = collision.foundZDC().energySectorZNA()[3];
      float E1ZNC = collision.foundZDC().energySectorZNC()[0];
      float E2ZNC = collision.foundZDC().energySectorZNC()[1];
      float E3ZNC = collision.foundZDC().energySectorZNC()[2];
      float E4ZNC = collision.foundZDC().energySectorZNC()[3];
      
      registry.get<TProfile2D>(HIST("hZNA_mean_common_cent"))->Fill(Form("%d",runNumber),int(cent),EcomZNA,1);
      registry.get<TProfile2D>(HIST("hZNA_mean_t1_cent"))->Fill(Form("%d", runNumber),int(cent),E1ZNA,1);
      registry.get<TProfile2D>(HIST("hZNA_mean_t2_cent"))->Fill(Form("%d", runNumber),int(cent),E2ZNA,1);
      registry.get<TProfile2D>(HIST("hZNA_mean_t3_cent"))->Fill(Form("%d", runNumber),int(cent),E3ZNA,1);
      registry.get<TProfile2D>(HIST("hZNA_mean_t4_cent"))->Fill(Form("%d", runNumber),int(cent),E4ZNA,1);
      
      registry.get<TProfile2D>(HIST("hZNC_mean_common_cent"))->Fill(Form("%d", runNumber),int(cent),EcomZNC,1);
      registry.get<TProfile2D>(HIST("hZNC_mean_t1_cent"))->Fill(Form("%d", runNumber),int(cent),E1ZNC,1);
      registry.get<TProfile2D>(HIST("hZNC_mean_t2_cent"))->Fill(Form("%d", runNumber),int(cent),E2ZNC,1);
      registry.get<TProfile2D>(HIST("hZNC_mean_t3_cent"))->Fill(Form("%d", runNumber),int(cent),E3ZNC,1);
      registry.get<TProfile2D>(HIST("hZNC_mean_t4_cent"))->Fill(Form("%d", runNumber),int(cent),E4ZNC,1);
    }
    
  }
  
  PROCESS_SWITCH(ZDCAnalysis, processCalib, "Process Calibration", true);
  
  void processData(myCollisions::iterator const& collision, aod::Zdcs const& zdcs, aod::BCs const& bcs, myTracks const& tracks)
  //Belangrijk dus dat je wel subscribed hier naar aod::Zdcs om toegang te krijgen tot zdc informatie
  {
    //https://alice-notes.web.cern.ch/system/files/notes/analysis/620/2017-May-31-analysis_note-ALICE_analysis_note_v2.pdf
    float ZDC_px[4] = {-1.75, 1.75,-1.75, 1.75};
    float ZDC_py[4] = {-1.75, -1.75,1.75, 1.75};
    float alpha = 0.395;
    //tot hier van die analysis note

    float sumZNA = 0;
    float sumZNC = 0;

    float xEnZNA = 0;
    float xEnZNC = 0;
    float yEnZNA = 0;
    float yEnZNC = 0;

    float QxA = 0;
    float QyA = 0;
    float QxC = 0;
    float QyC = 0;

    int Ntot = tracks.size();
    if (Ntot < 1)
      return;
    if (!collision.sel7())
      return;
    auto cent = collision.centRun2V0M();
    if(cent<0 || cent>100)
      return;

    registry.fill(HIST("hVtxX"), collision.posX());
    registry.fill(HIST("hVtxY"), collision.posY());
    registry.fill(HIST("hVtxZ"), collision.posZ());
    registry.fill(HIST("hMult"), Ntot);
    registry.fill(HIST("hCent"), cent);

    if (collision.foundZDCId() >= 0) {

      // check if we have the right gaineq files

      int runNumber = collision.bc().runNumber();
      //      LOGF(info, "RUNNUMBER = %i", runNumber);


      registry.get<TH1>(HIST("ZNAcoll"))->Fill(collision.foundZDC().amplitudeZNA());
      registry.get<TH1>(HIST("ZNCcoll"))->Fill(collision.foundZDC().amplitudeZNC());
      registry.get<TH2>(HIST("ZNAvsZNCcoll"))->Fill(collision.foundZDC().amplitudeZNC(), collision.foundZDC().amplitudeZNA());
      registry.get<TH2>(HIST("ZNCadcvstdccoll"))->Fill(collision.foundZDC().timeZNC(), collision.foundZDC().amplitudeZNC());
      registry.get<TH2>(HIST("ZNAadcvstdccoll"))->Fill(collision.foundZDC().timeZNA(), collision.foundZDC().amplitudeZNA());

      float EcomZNA = collision.foundZDC().energyCommonZNA();
      float EcomZNC = collision.foundZDC().energyCommonZNC();
      float E1ZNA = collision.foundZDC().energySectorZNA()[0];
      float E2ZNA = collision.foundZDC().energySectorZNA()[1];
      float E3ZNA = collision.foundZDC().energySectorZNA()[2];
      float E4ZNA = collision.foundZDC().energySectorZNA()[3];
      float E1ZNC = collision.foundZDC().energySectorZNC()[0];
      float E2ZNC = collision.foundZDC().energySectorZNC()[1];
      float E3ZNC = collision.foundZDC().energySectorZNC()[2];
      float E4ZNC = collision.foundZDC().energySectorZNC()[3];

      registry.fill(HIST("hZNA_common_cent"), int(cent),EcomZNA);
      registry.fill(HIST("hZNA_t1_cent"), int(cent),E1ZNA);
      registry.fill(HIST("hZNA_t2_cent"), int(cent),E2ZNA);
      registry.fill(HIST("hZNA_t3_cent"), int(cent),E3ZNA);
      registry.fill(HIST("hZNA_t4_cent"), int(cent),E4ZNA);

      registry.fill(HIST("hZNC_common_cent"), int(cent),EcomZNC);
      registry.fill(HIST("hZNC_t1_cent"), int(cent),E1ZNC);
      registry.fill(HIST("hZNC_t2_cent"), int(cent),E2ZNC);
      registry.fill(HIST("hZNC_t3_cent"), int(cent),E3ZNC);
      registry.fill(HIST("hZNC_t4_cent"), int(cent),E4ZNC);

      if(!f || f->IsZombie()|| !f->IsOpen()) {
        LOGF(fatal, "Mean file does not exist | run in processCalib mode!!!!");
        return;
      }

//      TList* list = (TList*)f->Get("z-d-c-analysis");
      TProfile2D* hZNA_centCom = (TProfile2D*)f->Get("z-d-c-analysis/hZNA_mean_common_cent");
      
      if(!hZNA_centCom){
        LOGF(error, "Histo not found! Abort mission");
        return; }
//      else{LOGF(info, "Histo found!");}
      
      if (runNumber != lastRunNumber) {
        lastRunNumber = runNumber;
        for(Int_t i=0; i<hZNA_centCom->GetXaxis()->GetNbins(); i++){
          const char* label = hZNA_centCom->GetXaxis()->GetBinLabel(i+1);
//          LOGF(info, "RUN NUMBER %i | looking in bin %s", runNumber, label);
          int labelInt = atoi(label);
          if(labelInt == runNumber){bin=i+1;}
//          LOGF(info, "looking in bin %i", bin);
        }
        numFoundint = atoi(hZNA_centCom->GetXaxis()->GetBinLabel(bin));
        LOGF(info, "RUN NUMBER %i | looking in bin %i : (%i)", runNumber, bin, numFoundint);
      }
        if(numFoundint!=runNumber)LOGF(error, "No match found for RUN NUMBER %i", runNumber);
      
      float ZNA_centCom = hZNA_centCom->GetBinContent(int(bin),int(cent)+1);
      TProfile2D* hZNC_centCom= (TProfile2D*)f->Get(Form("z-d-c-analysis/hZNC_mean_common_cent"));
      float ZNC_centCom = hZNC_centCom->GetBinContent(int(bin),int(cent)+1);

      float e1A=0;
      float e2A=0;
      float e3A=0;
      float e4A=0;

      float e1C=0;
      float e2C=0;
      float e3C=0;
      float e4C=0;
      
      

      TProfile2D* hmeanE1A= (TProfile2D*)f->Get("z-d-c-analysis/hZNA_mean_t1_cent");
      float meanE1A = hmeanE1A->GetBinContent(int(bin), int(cent)+1);
      TProfile2D* hmeanE2A= (TProfile2D*)f->Get("z-d-c-analysis/hZNA_mean_t2_cent");
      float meanE2A = hmeanE2A->GetBinContent(int(bin), int(cent)+1);
      TProfile2D* hmeanE3A= (TProfile2D*)f->Get("z-d-c-analysis/hZNA_mean_t3_cent");
      float meanE3A = hmeanE3A->GetBinContent(int(bin), int(cent)+1);
      TProfile2D* hmeanE4A= (TProfile2D*)f->Get("z-d-c-analysis/hZNA_mean_t4_cent");
      float meanE4A = hmeanE4A->GetBinContent(int(bin), int(cent)+1);

      TProfile2D* hmeanE1C= (TProfile2D*)f->Get("z-d-c-analysis/hZNC_mean_t1_cent");
      float meanE1C = hmeanE1C->GetBinContent(int(bin), int(cent)+1);
      TProfile2D* hmeanE2C= (TProfile2D*)f->Get("z-d-c-analysis/hZNC_mean_t2_cent");
       float meanE2C = hmeanE2C->GetBinContent(int(bin), int(cent)+1);
      TProfile2D* hmeanE3C= (TProfile2D*)f->Get("z-d-c-analysis/hZNC_mean_t3_cent");
      float meanE3C = hmeanE3C->GetBinContent(int(bin), int(cent)+1);
      TProfile2D* hmeanE4C= (TProfile2D*)f->Get("z-d-c-analysis/hZNC_mean_t4_cent");
      float meanE4C = hmeanE4C->GetBinContent(int(bin), int(cent)+1);

//      f->Close();
      
      if(meanE1A > 0)  e1A = E1ZNA * (0.25*ZNA_centCom) / meanE1A;
      if(meanE2A > 0)  e2A = E2ZNA * (0.25*ZNA_centCom) / meanE2A;
      if(meanE3A > 0)  e3A = E3ZNA * (0.25*ZNA_centCom) / meanE3A;
      if(meanE4A > 0)  e4A = E4ZNA * (0.25*ZNA_centCom) / meanE4A;

      if(meanE1C > 0)  e1C = E1ZNC * (0.25*ZNC_centCom) / meanE1C;
      if(meanE2C > 0)  e2C = E2ZNC * (0.25*ZNC_centCom) / meanE2C;
      if(meanE3C > 0)  e3C = E3ZNC * (0.25*ZNC_centCom) / meanE3C;
      if(meanE4C > 0)  e4C = E4ZNC * (0.25*ZNC_centCom) / meanE4C;

      //      if(meanE1A > 0)LOGF(info, "meanE1A = %f cent is %f",meanE1A,cent);
      //      if(meanE2A > 0)LOGF(info, "meanE2A = %f cent is %f",meanE2A,cent);
      //      if(meanE3A > 0)LOGF(info, "meanE3A = %f cent is %f",meanE3A,cent);
      //      if(meanE4A > 0)LOGF(info, "meanE4A = %f cent is %f",meanE4A,cent);
      //      if(meanE1C > 0)LOGF(info, "meanE1C = %f cent is %f",meanE1C,cent);
      //      if(meanE2C > 0)LOGF(info, "meanE2C = %f cent is %f",meanE2C,cent);
      //      if(meanE3C > 0)LOGF(info, "meanE3C = %f cent is %f",meanE3C,cent);
      //      if(meanE4C > 0)LOGF(info, "meanE4C = %f cent is %f",meanE4C,cent);

      //      LOGF(info, "Energy was = %f | Energy is = %f", collision.foundZDC().energySectorZNA()[0], e1A);
      //      LOGF(info, "Ecom/4/<Ej> =  %f",(0.25*ZNA_centCom) / registry.get<TH1>(HIST("hZNA_mean_t1_cent"))->GetBinContent(int(cent)+1));
      std::vector<float> calE = {e1A, e2A, e3A , e4A, e1C, e2C, e3C , e4C};

      // pas op: gebruik niet deze info om door te rekenen maar e1 etc.

      registry.fill(HIST("hZNA_t1_cent_cal"), int(cent), calE[0]);
      registry.fill(HIST("hZNA_t2_cent_cal"), int(cent), calE[1]);
      registry.fill(HIST("hZNA_t3_cent_cal"), int(cent), calE[2]);
      registry.fill(HIST("hZNA_t4_cent_cal"), int(cent), calE[3]);

      registry.fill(HIST("hZNC_t1_cent_cal"), int(cent), calE[4]);
      registry.fill(HIST("hZNC_t2_cent_cal"), int(cent), calE[5]);
      registry.fill(HIST("hZNC_t3_cent_cal"), int(cent), calE[6]);
      registry.fill(HIST("hZNC_t4_cent_cal"), int(cent), calE[7]);

      for(int i=0; i<4; i++){
        // for full centrality
        registry.fill(HIST("hAllTowersEnergy_noEq"),i,collision.foundZDC().energySectorZNA()[i]);
        registry.fill(HIST("hAllTowersEnergy_noEq"),i+4,collision.foundZDC().energySectorZNC()[i]);
        registry.fill(HIST("hAllTowersEnergy"),i,calE[i]); //fills from ZNA towers
        registry.fill(HIST("hAllTowersEnergy"),i+4,calE[i+4]); //fills from ZNC towers

        //for special centrality bin use output of hZN*_*_cent and make projection or somthing.
        //No need to implement that over here!

        //With calibrated energies!!
        sumZNA += TMath::Power(calE[i],alpha);
        sumZNC += TMath::Power(calE[i+4],alpha);

        xEnZNA += ZDC_px[i] * TMath::Power(calE[i],alpha);
        yEnZNA += ZDC_py[i] * TMath::Power(calE[i],alpha);

        xEnZNC += ZDC_px[i] * TMath::Power(calE[i+4],alpha);
        yEnZNC += ZDC_py[i] * TMath::Power(calE[i+4],alpha);

      } // end of for(int i=0; i<4; i++)


      QxA = xEnZNA / sumZNA ;
      QyA = yEnZNA / sumZNA ;

      QxC = xEnZNC / sumZNC ;
      QyC = yEnZNC / sumZNC ;

      registry.fill(HIST("hZNA_Qx"), QxA);
      registry.fill(HIST("hZNA_Qy"), QyA);

      registry.fill(HIST("hZNC_Qx"), QxC);
      registry.fill(HIST("hZNC_Qy"), QyC);

      registry.fill(HIST("hZNCvsZNA_Qx"), QxC, QxA);
      registry.fill(HIST("hZNCvsZNA_Qy"), QyC, QyA);

      registry.fill(HIST("hZNAC_common"),0,collision.foundZDC().energyCommonZNA());
      registry.fill(HIST("hZNAC_common"),1,collision.foundZDC().energyCommonZNC());

      registry.fill(HIST("hZNA_s1"),collision.foundZDC().energySectorZNA()[0]);
      registry.fill(HIST("hZNC_s1"),collision.foundZDC().energySectorZNC()[0]);

      registry.fill(HIST("hZNA_s2"),collision.foundZDC().energySectorZNA()[1]);
      registry.fill(HIST("hZNC_s2"),collision.foundZDC().energySectorZNC()[1]);

      registry.fill(HIST("hZNA_s3"),collision.foundZDC().energySectorZNA()[2]);
      registry.fill(HIST("hZNC_s3"),collision.foundZDC().energySectorZNC()[2]);

      registry.fill(HIST("hZNA_s4"),collision.foundZDC().energySectorZNA()[3]);
      registry.fill(HIST("hZNC_s4"),collision.foundZDC().energySectorZNC()[3]);

      registry.fill(HIST("hZNA_common"),collision.foundZDC().energyCommonZNA(),collision.foundZDC().energySectorZNA()[0]+collision.foundZDC().energySectorZNA()[1]+collision.foundZDC().energySectorZNA()[2]+collision.foundZDC().energySectorZNA()[3]);
      registry.fill(HIST("hZNC_common"),collision.foundZDC().energyCommonZNC(),collision.foundZDC().energySectorZNC()[0]+collision.foundZDC().energySectorZNC()[1]+collision.foundZDC().energySectorZNC()[2]+collision.foundZDC().energySectorZNC()[3]);

      registry.fill(HIST("hZNA_Qx_vs_vx"),QxA,collision.posX());
      registry.fill(HIST("hZNA_Qy_vs_vx"),QyA,collision.posX());
      registry.fill(HIST("hZNA_Qx_vs_vy"),QxA,collision.posY());
      registry.fill(HIST("hZNA_Qy_vs_vy"),QyA,collision.posY());
      registry.fill(HIST("hZNA_Qx_vs_vz"),QxA,collision.posZ());
      registry.fill(HIST("hZNA_Qy_vs_vz"),QyA,collision.posZ());

      registry.fill(HIST("hZNC_Qx_vs_vx"),QxC,collision.posX());
      registry.fill(HIST("hZNC_Qy_vs_vx"),QyC,collision.posX());
      registry.fill(HIST("hZNC_Qx_vs_vy"),QxC,collision.posY());
      registry.fill(HIST("hZNC_Qy_vs_vy"),QyC,collision.posY());
      registry.fill(HIST("hZNC_Qx_vs_vz"),QxC,collision.posZ());
      registry.fill(HIST("hZNC_Qy_vs_vz"),QyC,collision.posZ());

      registry.fill(HIST("hZNA_Qx_vs_Centrality"),QxA,cent);
      registry.fill(HIST("hZNA_Qy_vs_Centrality"),QyA,cent);
      registry.fill(HIST("hZNC_Qx_vs_Centrality"),QxC,cent);
      registry.fill(HIST("hZNC_Qy_vs_Centrality"),QyC,cent);

    }
    //    }
    //
    for(auto track : tracks){
      registry.fill(HIST("hPhi"), track.phi());
      registry.fill(HIST("hEta"), track.eta());
      registry.fill(HIST("hPt"), track.pt());
    }
//    f->Close();
    //    if(int(cent)==1) {
    //      LOGF(info, "*** E4ZNC[1] = %f |  E4ZNCcounter[1] = %i", E4ZNC_mean[1], E4ZNCcounter[1]);
    //    }

  }

  PROCESS_SWITCH(ZDCAnalysis, processData, "Process data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ZDCAnalysis>(cfgc)};
  
}
