
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
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "CCDB/BasicCCDBManager.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "ReconstructionDataFormats/Track.h"

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
#include "TF1.h"
#include "TFitResult.h"
#include "TCanvas.h"
#include <iostream>
#include "TSystem.h"
#include <TROOT.h>


#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;
using namespace o2::aod::evsel;

//define my.....
using myCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>>;
//run2
//soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>>;
using myTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>>;
using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;

namespace o2::analysis::myanalysistask{

bool step0_open=false;
bool step1_open=false;
bool step2_open=false;
bool step3_open=false;
bool step4_open=false;
bool step5_open=false;

int binRunNumber=-1;
int numFoundint = -1;
int lastRunNumber = -1;

int bin1=0;
int lastRunNumber_step1 = -1;
int counter =0;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


//here for each step open histos (if possible) and set bool for each step to true. Then in process check bool and proceed or start collecting q-vectors

//step0
std::vector<TProfile2D*> hZN_mean(10, nullptr);
//step1
std::vector<TProfile2D*> mean_1perCent_Run(4, nullptr);// hQXA, hQYA, hQXC, hQYC
//step2
std::vector<THnSparseD*> mean_10perCent_v(4, nullptr); // hQXA, hQYA, hQXC, hQYC
std::vector<TProfile*> v_mean_Run(3, nullptr); //hvx, hvy, hvz
//step3 [fit]
std::vector<TH2D*> hFit_Run_Cent_vx(2, nullptr); // hQYA, hQYC
std::vector<TH2D*> hFit_Run_Cent_vy(2, nullptr); // hQXA, hQXC
std::vector<TH2D*> hFit_Run_Cent_vz(4, nullptr); // hQXA, hQYA, hQXC, hQYC
std::vector<THnSparseD*> hQXA_Run_10perCent_vy(1, nullptr); //for check bin step3
//step4
std::vector<TProfile2D*> hmeanN_1perCent_Run(1, nullptr);
std::vector<THnSparseD*> mean_Magnet_10perCent_v(4, nullptr);
//step5
std::vector<THnSparseD*> mean_run_cent10_Mult(4, nullptr);

}

using namespace o2::analysis::myanalysistask;

struct ZDCqvectors{
  ConfigurableAxis axisCent{"axisCent", {90,0,90},"Centrality axis in 1% bins"};
  ConfigurableAxis axisCent10{"axisCent10", {9,0,90},"Centrality axis in 10% bins"};
  ConfigurableAxis axisQ{"axisQ", {100,-2,2},"Q vector (xy) in ZDC"};
  ConfigurableAxis axisVx{"axisVx", {100,0,0.15},"for Pos X of collision"};
  ConfigurableAxis axisVy{"axisVy", {100,0.35,0.41},"for Pos Y of collision"};
  ConfigurableAxis axisVx_mean{"axisVx_mean", {100,-0.01,0.01},"for Pos X of collision"};
  ConfigurableAxis axisVy_mean{"axisVy_mean", {100,-0.01,0.01},"for Pos Y of collision"};
  ConfigurableAxis axisVz{"axisVz", {100,-12,12},"for vz of collision"}; // take 12 because we shift vi - <vi>
  ConfigurableAxis axisRun{"axisRun", {1e6,0,1e6},"for runNumber in ThnSparse"};
  ConfigurableAxis axisPolarity{"axisPolarity", {2,-1,1},"Magnet Polarity"};
  ConfigurableAxis axisMult{"axisMult", {5000,0,5000},"Track Multiplicity"};
  ConfigurableAxis axisMult_mean{"axisMult_mean", {5000,-2500,2500},"Track Multiplicity mean per run"};
  
  O2_DEFINE_CONFIGURABLE(cfgCutVertex,        float,      10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMin,      float,      0.2f, "Minimal.q pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMax,      float,      10.0f, "Maximal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin,         float,      0.2f, "Minimal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax,         float,      3.0f, "Maximal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta,           float,      0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls,  float,      2.5, "Chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgUseNch,           bool,       false, "Use Nch for flow observables")
  O2_DEFINE_CONFIGURABLE(cfgNbootstrap,       int,        10, "Number of subsamples")
  O2_DEFINE_CONFIGURABLE(cfgMagField,         float,      99999, "Configurable magnetic field; default CCDB will be queried")
  O2_DEFINE_CONFIGURABLE(cfgQvecrecent,      std::string, "AnalysisResults_outCalib.root", "path+file from recentring")
  O2_DEFINE_CONFIGURABLE(cfgCalibenergy,      std::string, "../runCalib/AnalysisResults.root", "path+file from energy calibration")
  O2_DEFINE_CONFIGURABLE(cfgOutFitStep3,      std::string, "outFitStep3.root", "path+file from energy calibration")
  
  double centrality;
  double trackmultiplicity;
  double trackmultiplicity_mean;
  int runnumber;
  double polarity;
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //  //define output
  HistogramRegistry registry{"Registry"};
  
  //  Filters
  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtMin) && (aod::track::pt < cfgCutPtMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls);
  
  
  void init(InitContext const&)
  {
    
    //Tower mean energies vs. centrality
    registry.add("meanE/hZNA_mean_t1_cent","", {HistType::kTProfile, {{axisCent}}});
    registry.add("meanE/hZNA_mean_t2_cent","", {HistType::kTProfile, {{axisCent}}});
    registry.add("meanE/hZNA_mean_t3_cent","", {HistType::kTProfile, {{axisCent}}});
    registry.add("meanE/hZNA_mean_t4_cent","", {HistType::kTProfile, {{axisCent}}});
    registry.add("meanE/hZNA_mean_common_cent","", {HistType::kTProfile, {{axisCent}}});
    
    registry.add("meanE/hZNC_mean_t1_cent","", {HistType::kTProfile, {{axisCent}}});
    registry.add("meanE/hZNC_mean_t2_cent","", {HistType::kTProfile, {{axisCent}}});
    registry.add("meanE/hZNC_mean_t3_cent","", {HistType::kTProfile, {{axisCent}}});
    registry.add("meanE/hZNC_mean_t4_cent","", {HistType::kTProfile, {{axisCent}}});
    registry.add("meanE/hZNC_mean_common_cent","", {HistType::kTProfile, {{axisCent}}});
    
    // before
    registry.add("before/hZNA_Qx_vs_Qy","hZNA_Qx_vs_Qy", {HistType::kTH2F, {axisQ,axisQ}});
    registry.add("before/hZNC_Qx_vs_Qy","hZNC_Qx_vs_Qy", {HistType::kTH2F, {axisQ,axisQ}});
    
    
    registry.add("before/hQXA_QXC_vs_cent","hQXA_QXC_vs_cent", {HistType::kTProfile, {axisCent10}});
    registry.add("before/hQYA_QYC_vs_cent","hQYA_QYC_vs_cent", {HistType::kTProfile, {axisCent10}});
    registry.add("before/hQXC_QYA_vs_cent","hQXC_QYA_vs_cent", {HistType::kTProfile, {axisCent10}});
    registry.add("before/hQYC_QXA_vs_cent","hQYC_QXA_vs_cent", {HistType::kTProfile, {axisCent10}});
    
    //step1
    registry.add("step1/hZNA_Qx_vs_Qy","hZNA_Qx_vs_Qy", {HistType::kTH2F, {axisQ,axisQ}});
    registry.add("step1/hZNC_Qx_vs_Qy","hZNC_Qx_vs_Qy", {HistType::kTH2F, {axisQ,axisQ}});
    
    registry.add("step1/hQXA_mean_1perCent_Run","hQXA_mean_1perCent_Run", {HistType::kTProfile2D, {{1,0.,1.} ,axisCent}});
    registry.add("step1/hQYA_mean_1perCent_Run","hQYA_mean_1perCent_Run", {HistType::kTProfile2D, {{1,0.,1.} ,axisCent}});
    registry.add("step1/hQXC_mean_1perCent_Run","hQXC_mean_1perCent_Run", {HistType::kTProfile2D, {{1,0.,1.} ,axisCent}});
    registry.add("step1/hQYC_mean_1perCent_Run","hQYC_mean_1perCent_Run", {HistType::kTProfile2D, {{1,0.,1.} ,axisCent}});

    //step2
    
    registry.add("step2/hZNA_Qx_vs_Qy","hZNA_Qx_vs_Qy", {HistType::kTH2F, {axisQ,axisQ}});
    registry.add("step2/hZNC_Qx_vs_Qy","hZNC_Qx_vs_Qy", {HistType::kTH2F, {axisQ,axisQ}});
    
    registry.add("step2/hQXA_mean_10perCent_v","hQXA_mean_10perCent_v", {HistType::kTHnSparseD, {axisCent10, axisVx_mean, axisVy_mean, axisVz, axisQ}});
    registry.add("step2/hQYA_mean_10perCent_v","hQYA_mean_10perCent_v", {HistType::kTHnSparseD, {axisCent10, axisVx_mean, axisVy_mean, axisVz, axisQ}});
    registry.add("step2/hQXC_mean_10perCent_v","hQXC_mean_10perCent_v", {HistType::kTHnSparseD, {axisCent10, axisVx_mean, axisVy_mean, axisVz, axisQ}});
    registry.add("step2/hQYC_mean_10perCent_v","hQYC_mean_10perCent_v", {HistType::kTHnSparseD, {axisCent10, axisVx_mean, axisVy_mean, axisVz, axisQ}});
    
    registry.add("step2/hvx_mean_Run","hvx_mean_Run", {HistType::kTProfile, {{1,0.,1.}}});
    registry.add("step2/hvy_mean_Run","hvy_mean_Run", {HistType::kTProfile, {{1,0.,1.}}});
    registry.add("step2/hvz_mean_Run","hvz_mean_Run", {HistType::kTProfile, {{1,0.,1.}}});
    
    //step 3
    registry.add("step3/hZNA_Qx_vs_Qy","hZNA_Qx_vs_Qy", {HistType::kTH2F, {axisQ,axisQ}});
    registry.add("step3/hZNC_Qx_vs_Qy","hZNC_Qx_vs_Qy", {HistType::kTH2F, {axisQ,axisQ}});
    
    registry.add("step3/hQXA_mean_10perCent_Run","hQXA_mean_10perCent_Run", {HistType::kTProfile2D, {{1,0.,1.} ,axisCent10}});
    registry.add("step3/hQYA_mean_10perCent_Run","hQYA_mean_10perCent_Run", {HistType::kTProfile2D, {{1,0.,1.} ,axisCent10}});
    registry.add("step3/hQXC_mean_10perCent_Run","hQXC_mean_10perCent_Run", {HistType::kTProfile2D, {{1,0.,1.} ,axisCent10}});
    registry.add("step3/hQYC_mean_10perCent_Run","hQYC_mean_10perCent_Run", {HistType::kTProfile2D, {{1,0.,1.} ,axisCent10}});
    
    // cross terms
    registry.add("step3/hQXA_Run_10perCent_vy","hQXA_Run_10perCent_vy", {HistType::kTHnSparseD, {axisRun,axisCent10,axisVy_mean, axisQ}});
    registry.add("step3/hQXA_Run_10perCent_vz","hQXA_Run_10perCent_vz", {HistType::kTHnSparseD, {axisRun,axisCent10,axisVz, axisQ}});
    
    registry.add("step3/hQYA_Run_10perCent_vx","hQYA_Run_10perCent_vx", {HistType::kTHnSparseD, {axisRun,axisCent10,axisVx_mean, axisQ}});
    registry.add("step3/hQYA_Run_10perCent_vz","hQYA_Run_10perCent_vz", {HistType::kTHnSparseD, {axisRun,axisCent10,axisVz, axisQ}});
    
    registry.add("step3/hQXC_Run_10perCent_vy","hQXC_Run_10perCent_vx", {HistType::kTHnSparseD, {axisRun,axisCent10,axisVy_mean, axisQ}});
    registry.add("step3/hQXC_Run_10perCent_vz","hQXC_Run_10perCent_vz", {HistType::kTHnSparseD, {axisRun,axisCent10,axisVz, axisQ}});
    
    registry.add("step3/hQYC_Run_10perCent_vx","hQYC_Run_10perCent_vx", {HistType::kTHnSparseD, {axisRun,axisCent10,axisVx_mean, axisQ}});
    registry.add("step3/hQYC_Run_10perCent_vz","hQYC_Run_10perCent_vz", {HistType::kTHnSparseD, {axisRun,axisCent10,axisVz, axisQ}});
    

    
    //step 4
    registry.add("step4/hZNA_Qx_vs_Qy","hZNA_Qx_vs_Qy", {HistType::kTH2F, {axisQ,axisQ}});
    registry.add("step4/hZNC_Qx_vs_Qy","hZNC_Qx_vs_Qy", {HistType::kTH2F, {axisQ,axisQ}});
    
    registry.add("step4/hmeanN_1perCent_Run","hmeanN_1perCent_Run", {HistType::kTProfile2D, {{1,0.,1.} ,axisCent}});
    
    registry.add("step4/hQXA_mean_Magnet_10perCent_v","hQXA_mean_Magnet_10perCent_v", {HistType::kTHnSparseD, {axisPolarity, axisCent10, axisVx_mean, axisVy_mean, axisVz, axisQ}});
    registry.add("step4/hQYA_mean_Magnet_10perCent_v","hQYA_mean_Magnet_10perCent_v", {HistType::kTHnSparseD, {axisPolarity, axisCent10, axisVx_mean, axisVy_mean, axisVz, axisQ}});
    registry.add("step4/hQXC_mean_Magnet_10perCent_v","hQXC_mean_Magnet_10perCent_v", {HistType::kTHnSparseD, {axisPolarity, axisCent10, axisVx_mean, axisVy_mean, axisVz, axisQ}});
    registry.add("step4/hQYC_mean_Magnet_10perCent_v","hQYC_mean_Magnet_10perCent_v", {HistType::kTHnSparseD, {axisPolarity, axisCent10, axisVx_mean, axisVy_mean, axisVz, axisQ}});
    
    //step 5
    registry.add("step5/hZNA_Qx_vs_Qy","hZNA_Qx_vs_Qy", {HistType::kTH2F, {axisQ,axisQ}});
    registry.add("step5/hZNC_Qx_vs_Qy","hZNC_Qx_vs_Qy", {HistType::kTH2F, {axisQ,axisQ}});
    
    registry.add("step5/hQXA_mean_run_cent10_Mult","hQXA_mean_run_cent10_Mult", {HistType::kTHnSparseD, {axisRun, axisCent10, axisMult_mean, axisQ}});
    registry.add("step5/hQYA_mean_run_cent10_Mult","hQYA_mean_run_cent10_Mult", {HistType::kTHnSparseD, {axisRun, axisCent10, axisMult_mean, axisQ}});
    registry.add("step5/hQXC_mean_run_cent10_Mult","hQXC_mean_run_cent10_Mult", {HistType::kTHnSparseD, {axisRun, axisCent10, axisMult_mean, axisQ}});
    registry.add("step5/hQYC_mean_run_cent10_Mult","hQYC_mean_run_cent10_Mult", {HistType::kTHnSparseD, {axisRun, axisCent10, axisMult_mean, axisQ}});
    
    registry.add("step5/hQXA_QXC_vs_cent","hQXA_QXC_vs_cent", {HistType::kTProfile, {axisCent10}});
    registry.add("step5/hQYA_QYC_vs_cent","hQYA_QYC_vs_cent", {HistType::kTProfile, {axisCent10}});
    registry.add("step5/hQXC_QYA_vs_cent","hQXC_QYA_vs_cent", {HistType::kTProfile, {axisCent10}});
    registry.add("step5/hQYC_QXA_vs_cent","hQYC_QXA_vs_cent", {HistType::kTProfile, {axisCent10}});
    
    // recentered q-vectors
    registry.add("hStep","hStep", {HistType::kTH1D, {{5,0.,5.}}});
    
  }
  
//  template <typename CollisionObject>
//  inline void fillMeanRegistry(double cent, CollisionObject collision){
//    if constexpr(framework::has_type_v<aod::cent::CentRun2V0M, typename CollisionObject::all_columns>) {
//      registry.fill(HIST("hZNA_mean_common_cent"),int(cent),collision.foundZDC().energyCommonZNA());
//      registry.fill(HIST("hZNA_mean_t1_cent"),int(cent),collision.foundZDC().energySectorZNA()[0]);
//      registry.fill(HIST("hZNA_mean_t2_cent"),int(cent),collision.foundZDC().energySectorZNA()[1]);
//      registry.fill(HIST("hZNA_mean_t3_cent"),int(cent),collision.foundZDC().energySectorZNA()[2]);
//      registry.fill(HIST("hZNA_mean_t4_cent"),int(cent),collision.foundZDC().energySectorZNA()[3]);
//      
//      registry.fill(HIST("hZNC_mean_common_cent"),int(cent),collision.foundZDC().energyCommonZNC());
//      registry.fill(HIST("hZNC_mean_t1_cent"),int(cent),collision.foundZDC().energySectorZNC()[0]);
//      registry.fill(HIST("hZNC_mean_t2_cent"),int(cent),collision.foundZDC().energySectorZNC()[1]);
//      registry.fill(HIST("hZNC_mean_t3_cent"),int(cent),collision.foundZDC().energySectorZNC()[2]);
//      registry.fill(HIST("hZNC_mean_t4_cent"),int(cent),collision.foundZDC().energySectorZNC()[3]);
//      
//    } else if constexpr (framework::has_type_v<aod::cent::CentFT0C, typename CollisionObject::all_columns>){
//      registry.fill(HIST("hZNA_mean_common_cent"),int(cent),collision.foundBC_as<BCsRun3>().zdc().energyCommonZNA());
//      registry.fill(HIST("hZNA_mean_t1_cent"),int(cent),collision.foundBC_as<BCsRun3>().zdc().energySectorZNA()[0]);
//      registry.fill(HIST("hZNA_mean_t2_cent"),int(cent),collision.foundBC_as<BCsRun3>().zdc().energySectorZNA()[1]);
//      registry.fill(HIST("hZNA_mean_t3_cent"),int(cent),collision.foundBC_as<BCsRun3>().zdc().energySectorZNA()[2]);
//      registry.fill(HIST("hZNA_mean_t4_cent"),int(cent),collision.foundBC_as<BCsRun3>().zdc().energySectorZNA()[3]);
//      
//      registry.fill(HIST("hZNC_mean_common_cent"),int(cent),collision.foundBC_as<BCsRun3>().zdc().energyCommonZNC());
//      registry.fill(HIST("hZNC_mean_t1_cent"),int(cent),collision.foundBC_as<BCsRun3>().zdc().energySectorZNC()[0]);
//      registry.fill(HIST("hZNC_mean_t2_cent"),int(cent),collision.foundBC_as<BCsRun3>().zdc().energySectorZNC()[1]);
//      registry.fill(HIST("hZNC_mean_t3_cent"),int(cent),collision.foundBC_as<BCsRun3>().zdc().energySectorZNC()[2]);
//      registry.fill(HIST("hZNC_mean_t4_cent"),int(cent),collision.foundBC_as<BCsRun3>().zdc().energySectorZNC()[3]);
//    }
//  }
  //qxa, qya, qxc, qyc, vx, vy, vz, centrality, trackmultiplicity, runnumber, polarity
  inline void fillRegistry(int step, double qxa, double qya, double qxc, double qyc, double vx, double vy, double vz,  double centrality, double trackmultiplicity, int runnumber, double polarity)
  {
    //    LOGF(info, "centrality = %.2f", centrality);
    //TODO: dit kan korter door vector van histos te maken, kan dan ook in loopje worden gevuld.
    if(step==0){
      registry.get<TProfile>(HIST("before/hQXA_QXC_vs_cent"))->Fill(centrality,qxa*qxc);
      registry.get<TProfile>(HIST("before/hQYA_QYC_vs_cent"))->Fill(centrality,qya*qyc);
      registry.get<TProfile>(HIST("before/hQXC_QYA_vs_cent"))->Fill(centrality,qxc*qya);
      registry.get<TProfile>(HIST("before/hQYC_QXA_vs_cent"))->Fill(centrality,qyc*qxa);
      
      
      registry.fill(HIST("before/hZNA_Qx_vs_Qy"),qxa, qya);
      registry.fill(HIST("before/hZNC_Qx_vs_Qy"),qxc, qyc);
      
      registry.get<TProfile2D>(HIST("step1/hQXA_mean_1perCent_Run"))->Fill(Form("%d", runnumber),centrality,qxa,1);
      registry.get<TProfile2D>(HIST("step1/hQXC_mean_1perCent_Run"))->Fill(Form("%d", runnumber),centrality,qxc,1);
      registry.get<TProfile2D>(HIST("step1/hQYA_mean_1perCent_Run"))->Fill(Form("%d", runnumber),centrality,qya,1);
      registry.get<TProfile2D>(HIST("step1/hQYC_mean_1perCent_Run"))->Fill(Form("%d", runnumber),centrality,qyc,1);
      
    }
    
    if(step==1){
      registry.fill(HIST("step1/hZNA_Qx_vs_Qy"),qxa, qya);
      registry.fill(HIST("step1/hZNC_Qx_vs_Qy"),qxc, qyc);
      
      registry.fill(HIST("step2/hQXA_mean_10perCent_v"), centrality, vx, vy, vz, qxa);
      registry.fill(HIST("step2/hQYA_mean_10perCent_v"), centrality, vx, vy, vz, qya);
      registry.fill(HIST("step2/hQXC_mean_10perCent_v"), centrality, vx, vy, vz, qxc);
      registry.fill(HIST("step2/hQYC_mean_10perCent_v"), centrality, vx, vy, vz, qyc);
      
      registry.get<TProfile>(HIST("step2/hvx_mean_Run"))->Fill(Form("%d",runnumber), vx);
      registry.get<TProfile>(HIST("step2/hvy_mean_Run"))->Fill(Form("%d",runnumber), vy);
      registry.get<TProfile>(HIST("step2/hvz_mean_Run"))->Fill(Form("%d",runnumber), vz);
      
      registry.fill(HIST("hStep"), .5, 1);
    }
    
    if(step==2){
      
      registry.fill(HIST("step2/hZNA_Qx_vs_Qy"),qxa, qya);
      registry.fill(HIST("step2/hZNC_Qx_vs_Qy"),qxc, qyc);
      
      // cross terms -> Denk niet vz!
      registry.fill(HIST("step3/hQXA_Run_10perCent_vy"), runnumber, centrality, vy, qxa);
      registry.fill(HIST("step3/hQXA_Run_10perCent_vz"), runnumber, centrality, vz, qxa);
      registry.fill(HIST("step3/hQYA_Run_10perCent_vx"), runnumber, centrality, vx, qya);
      registry.fill(HIST("step3/hQYA_Run_10perCent_vz"), runnumber, centrality, vz, qya);
      registry.fill(HIST("step3/hQXC_Run_10perCent_vy"), runnumber, centrality, vy, qxc);
      registry.fill(HIST("step3/hQXC_Run_10perCent_vz"), runnumber, centrality, vz, qxc);
      registry.fill(HIST("step3/hQYC_Run_10perCent_vx"), runnumber, centrality, vx, qyc);
      registry.fill(HIST("step3/hQYC_Run_10perCent_vz"), runnumber, centrality, vz, qyc);
      
      // fill with zeros for fitting later.
      registry.get<TProfile2D>(HIST("step3/hQXA_mean_10perCent_Run"))->Fill(Form("%d", runnumber),centrality,0,1);
      registry.get<TProfile2D>(HIST("step3/hQXC_mean_10perCent_Run"))->Fill(Form("%d", runnumber),centrality,0,1);
      registry.get<TProfile2D>(HIST("step3/hQYA_mean_10perCent_Run"))->Fill(Form("%d", runnumber),centrality,0,1);
      registry.get<TProfile2D>(HIST("step3/hQYC_mean_10perCent_Run"))->Fill(Form("%d", runnumber),centrality,0,1);

      
      registry.fill(HIST("hStep"), 1.5, 1);
      
    }
    
    
    if(step==3){
      registry.fill(HIST("step3/hZNA_Qx_vs_Qy"),qxa, qya);
      registry.fill(HIST("step3/hZNC_Qx_vs_Qy"),qxc, qyc);
      
      registry.fill(HIST("step4/hQXA_mean_Magnet_10perCent_v"), polarity, centrality, vx, vy, vz, qxa);
      registry.fill(HIST("step4/hQYA_mean_Magnet_10perCent_v"), polarity, centrality, vx, vy, vz, qya);
      registry.fill(HIST("step4/hQXC_mean_Magnet_10perCent_v"), polarity, centrality, vx, vy, vz, qxc);
      registry.fill(HIST("step4/hQYC_mean_Magnet_10perCent_v"), polarity, centrality, vx, vy, vz, qyc);
      
      registry.get<TProfile2D>(HIST("step4/hmeanN_1perCent_Run"))->Fill(Form("%d",runnumber), centrality, trackmultiplicity,1);

      registry.get<TProfile>(HIST("step3/hQXA_vs_vx"))->Fill(vx, qxa);
      registry.get<TProfile>(HIST("step3/hQYA_vs_vx"))->Fill(vx, qya);
      registry.get<TProfile>(HIST("step3/hQXC_vs_vx"))->Fill(vx, qxc);
      registry.get<TProfile>(HIST("step3/hQYC_vs_vx"))->Fill(vx, qyc);
      
      registry.get<TProfile>(HIST("step3/hQXA_vs_vy"))->Fill(vy, qxa);
      registry.get<TProfile>(HIST("step3/hQYA_vs_vy"))->Fill(vy, qya);
      registry.get<TProfile>(HIST("step3/hQXC_vs_vy"))->Fill(vy, qxc);
      registry.get<TProfile>(HIST("step3/hQYC_vs_vy"))->Fill(vy, qyc);
      
      registry.get<TProfile>(HIST("step3/hQXA_vs_vz"))->Fill(vz, qxa);
      registry.get<TProfile>(HIST("step3/hQYA_vs_vz"))->Fill(vz, qya);
      registry.get<TProfile>(HIST("step3/hQXC_vs_vz"))->Fill(vz, qxc);
      registry.get<TProfile>(HIST("step3/hQYC_vs_vz"))->Fill(vz, qyc);
      
      registry.fill(HIST("hStep"), 2.5, 1);
    }
    
    if(step==4){
      
      registry.fill(HIST("step4/hZNA_Qx_vs_Qy"),qxa, qya);
      registry.fill(HIST("step4/hZNC_Qx_vs_Qy"),qxc, qyc);
      
      registry.fill(HIST("step5/hQXA_mean_run_cent10_Mult"), runnumber, centrality, trackmultiplicity, qxa);
      registry.fill(HIST("step5/hQYA_mean_run_cent10_Mult"), runnumber, centrality, trackmultiplicity, qya);
      registry.fill(HIST("step5/hQXC_mean_run_cent10_Mult"), runnumber, centrality, trackmultiplicity, qxc);
      registry.fill(HIST("step5/hQYC_mean_run_cent10_Mult"), runnumber, centrality, trackmultiplicity, qyc);
      
      registry.fill(HIST("hStep"), 3.5, 1);
    }
    
    //recentered q-vectors
    if(step==5){
      registry.fill(HIST("step5/hZNA_Qx_vs_Qy"),qxa, qya);
      registry.fill(HIST("step5/hZNC_Qx_vs_Qy"),qxc, qyc);
      
      registry.get<TProfile>(HIST("step5/hQXA_QXC_vs_cent"))->Fill(centrality,qxa*qxc);
      registry.get<TProfile>(HIST("step5/hQYA_QYC_vs_cent"))->Fill(centrality,qya*qyc);
      registry.get<TProfile>(HIST("step5/hQXC_QYA_vs_cent"))->Fill(centrality,qxc*qya);
      registry.get<TProfile>(HIST("step5/hQYC_QXA_vs_cent"))->Fill(centrality,qyc*qxa);
      
      registry.fill(HIST("hStep"), 4.5, 1);
    }
  }
  

  // Function to fill the registry
  void fillAllRegistries(double q[6][4], double v[3], double cent, double trackmult, int runnum, int polarity, int start = 0, int end = 5) {
    for (int i = start; i < end; ++i) {
      fillRegistry(i, q[i][0], q[i][1], q[i][2], q[i][3], v[0], v[1], v[2], cent, trackmult, runnum, polarity);
    }
  }
  
  template <typename T>
  bool retrieveAndCheck(std::vector<T*>& vec, const char* names[], const char* filename, const char* fullPath, const char* histname, bool closeFile = false) {
      
      bool file_open = true;
      TFile* file = nullptr;
      
      // Check if the file is already open
      if (gROOT->GetListOfFiles()->FindObject(filename)) {
          file = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
      } else {
          file = TFile::Open(filename, "READ");
      }
      
      if (!file || file->IsZombie() || !file->IsOpen()) {
        LOGF(info, "File %s does not exist!! Abort mission...", filename);
        file_open = false;
      }
      
      if (file_open) {
        for (int i = 0; names[i] != nullptr; ++i) {
            vec[i] = (T*)file->Get(Form("%s%s%s", fullPath, names[i], histname));
            if (!vec[i] || vec[i]->GetEntries() < 1) {
              LOGF(info, "%s%s_%s not found! Abort mission", fullPath, names[i], histname);
              if (closeFile) {
                file->Close();
                delete file;
              }
              return false;
            }
        }
      }
      
      if (closeFile) {
        file->Close();
        delete file;
      }
      
      return file_open;
  }
  
  void setupCalibration()
  {
    const char* names_hQ[4] = {"QXA", "QYA", "QXC", "QYC"};
    const char* names_vcoord[3] = {"vx", "vy", "vz"};
    const char* ZN_towers[10] = {"hZNA_mean_common", "hZNC_mean_common", "hZNA_mean_t1", "hZNA_mean_t2", "hZNA_mean_t3", "hZNA_mean_t4", "hZNC_mean_t1", "hZNC_mean_t2", "hZNC_mean_t3", "hZNC_mean_t4"};
    
    if(retrieveAndCheck(hZN_mean, ZN_towers, cfgCalibenergy->c_str(), "z-d-c-analysis/", "_cent", true)){
      step0_open = true;
    } else return;

    if (retrieveAndCheck(mean_1perCent_Run, names_hQ, cfgQvecrecent->c_str(), "z-d-cqvectors/step1/h","_mean_1perCent_Run", true)) {
      step1_open = true;
    } else return;

    if (retrieveAndCheck(mean_10perCent_v, names_hQ, cfgQvecrecent->c_str(), "z-d-cqvectors/step2/h", "_mean_10perCent_v") &&
        retrieveAndCheck(v_mean_Run, names_vcoord, cfgQvecrecent->c_str(), "z-d-cqvectors/step2/h", "_mean_Run") ){
      step2_open = true;
        } else return;
    
//    if (retrieveAndCheck(hFit_Run_Cent_vx,{names_hQ[1],names_hQ[3]},cfgOutFitStep3->c_str(),"hFit_","_Run_Cent_vx" )
//        && retrieveAndCheck(hFit_Run_Cent_vy,{names_hQ[0],names_hQ[2]},cfgOutFitStep3->c_str(),"hFit_","_Run_Cent_vy" )
//        && retrieveAndCheck(hFit_Run_Cent_vx,names_hQ,cfgOutFitStep3->c_str(),"hFit_","_Run_Cent_vx", true) // close cfgOutFitStep3
//        && retrieveAndCheck(hQXA_Run_10perCent_vy,{names_hQ[0]},cfgQvecrecent->c_str(),"z-d-cqvectors/step3/h","_Run_10perCent_vy")){
//      step3_open = true;
//    } else return;
//    
//    if (retrieveAndCheck(hmeanN_1perCent_Run, {""}, cfgQvecrecent->c_str(),"z-d-cqvectors/step4/","hmeanN_1perCent_Run" )
//        && retrieveAndCheck(mean_Magnet_10perCent_v, names_hQ, cfgQvecrecent->c_str(),"z-d-cqvectors/step4/h","_mean_Magnet_10perCent_v")){
//      step4_open=true;
//    } else return;
//    
//    if (retrieveAndCheck(mean_run_cent10_Mult, names_hQ, cfgQvecrecent->c_str(), "z-d-cqvectors/step5/h", "_mean_run_cent10_Mult", true)){
//      step5_open=true;
//    } else return;
    
    if(step0_open)LOGF(info, "histos energy calibration open! ");
    if(step1_open)LOGF(info, "histos step 1 open! ");
    if(step2_open)LOGF(info, "histos step 2 open! ");
    if(step3_open)LOGF(info, "histos step 3 open! ");
    if(step4_open)LOGF(info, "histos step 4 open! ");
    if(step5_open)LOGF(info, "histos step 5 open! ");

  }
  
  
  int FindBinWithContentAlongAxis(THnSparse *hist, Int_t axis, Double_t targetContent) {
    // to find bin for given axis in THnSparse object
    int bin =1;
    
    TAxis *ax = hist->GetAxis(axis);
    // Check if targetContent is within the range of the axis
    Double_t min = ax->GetXmin();
    Double_t max = ax->GetXmax();
    
    if (targetContent < min || targetContent > max) {
      LOGF(warning,"targetContent is out of range of axis ");
      return 0;
    }
    
    bin = ax->FindBin(targetContent);
    
    return bin;
  }
  
  
  double getQstep2(double q_step1, THnSparseD* mean, int cent, double vxmean, double vymean, double vzmean)
  {
    double mean_q =0.;
    
    
    // do we need to find the bins of v and i
    int bin_cent =   FindBinWithContentAlongAxis(mean, 0, cent);
    int bin_vxmean = FindBinWithContentAlongAxis(mean, 1, vxmean);
    int bin_vymean = FindBinWithContentAlongAxis(mean, 2, vymean);
    int bin_vzmean = FindBinWithContentAlongAxis(mean, 3, vzmean);
    
    mean->GetAxis(0)->SetRange(bin_cent, bin_cent);
    mean->GetAxis(1)->SetRange(bin_vxmean, bin_vxmean);
    mean->GetAxis(2)->SetRange(bin_vymean, bin_vymean);
    mean->GetAxis(3)->SetRange(bin_vzmean, bin_vzmean);
    
    mean_q = mean->Projection(4)->GetMean();
    
    double q_step2 = q_step1 - mean_q;
    
    return q_step2;
  }
  
  
  //  int getMagneticField(uint64_t timestamp)
  //  {
  //    // TODO done only once (and not per run). Will be replaced by CCDBConfigurable
  //    // static o2::parameters::GRPObject* grpo = nullptr;
  //    static o2::parameters::GRPMagField* grpo = nullptr;
  //    if (grpo == nullptr) {
  //      // grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", timestamp);
  //      grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
  //      if (grpo == nullptr) {
  //        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
  //        return 0;
  //      }
  //      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
  //    }
  //    return grpo->getNominalL3Field();
  //  }
  
  double getQstep4(double q_step3, THnSparseD* mean, double cent, double vxmean, double vymean, double vzmean, double pol)
  {
    
    double mean_q =0.;
    
    
    // do we need to find the bins of v and i
    int bin_pol = (pol>0) ? 2 : 1; // only 2 bins in histo, for positive and negative polarity
    int bin_cent =   FindBinWithContentAlongAxis(mean, 0, cent);
    int bin_vxmean = FindBinWithContentAlongAxis(mean, 1, vxmean);
    int bin_vymean = FindBinWithContentAlongAxis(mean, 2, vymean);
    int bin_vzmean = FindBinWithContentAlongAxis(mean, 3, vzmean);
    
    mean->GetAxis(0)->SetRange(bin_pol, bin_pol);
    mean->GetAxis(1)->SetRange(bin_cent, bin_cent);
    mean->GetAxis(2)->SetRange(bin_vxmean, bin_vxmean);
    mean->GetAxis(3)->SetRange(bin_vymean, bin_vymean);
    mean->GetAxis(4)->SetRange(bin_vzmean, bin_vzmean);
    
    mean_q = mean->Projection(5)->GetMean();
    //    LOGF(info, "nEntries proQ step2: %i", proQ->GetEntries());
    //    mean_q = proQ->GetMean();
    
    double q_step4 = q_step3 - mean_q;
    
    //    delete proQ;
    
    return q_step4;
  }
  
  double getQstep5(double q_step4, THnSparseD* mean, int bin_runnmbr, double cent, double mult)
  {
    //
    
    double mean_q =0.;
    
    
    //    int bin_runmbr = FindBinWithContentAlongAxis(mean, 0, runN);
    int bin_cent =   FindBinWithContentAlongAxis(mean, 1, cent);
    int bin_mult =   FindBinWithContentAlongAxis(mean, 2, mult);
    
    mean->GetAxis(0)->SetRange(bin_runnmbr, bin_runnmbr);
    mean->GetAxis(1)->SetRange(bin_cent, bin_cent);
    mean->GetAxis(2)->SetRange(bin_mult, bin_mult);
    
    
    mean_q = mean->Projection(3)->GetMean();
    
    double q_step5 = q_step4 - mean_q;
    
    //    delete proQ;
    
    return q_step5;
  }
  
  inline void PrintMemoryConsumption() {
    ProcInfo_t procInfo;
    gSystem->GetProcInfo(&procInfo);
    
    std::cout << "Memory Usage:" << std::endl;
    std::cout << "  Resident memory: " << (double)procInfo.fMemResident/(1e6) << " MB" << std::endl;
    std::cout << "  Virtual memory: " << (double)procInfo.fMemVirtual/(1e6) << " MB" << std::endl;
  }
  
//  <template T*>
  bool updateBinForRunNumber(int runNumber, int& lastRunNumber, int& binRunNumber, int& numFoundint, TProfile2D* hist) {

      if (runNumber != lastRunNumber) {
          lastRunNumber = runNumber;
          
          for (Int_t i = 0; i < hist->GetXaxis()->GetNbins(); i++) {
              const char* label = hist->GetXaxis()->GetBinLabel(i + 1);
              int labelInt = atoi(label);
              if (labelInt == runNumber) {
                  binRunNumber = i + 1;
                  break; // Exit the loop once the bin is found
              }
          }

        numFoundint = atoi(hist->GetXaxis()->GetBinLabel(binRunNumber));
        return true;
      }

      if (numFoundint != runNumber) {
          std::cerr << "No match found for RUN NUMBER " << runNumber << std::endl;
        return false;
      }
    return true;
  }
  
//TODO: void: make function for step 1! Also for step 2,3,4,5 And apply function in processRun2 and ProcessRun3
  
  void process(myCollisions const& cols,  BCsRun3 const& /*bcs*/, aod::Zdcs const& /*zdcs*/, myTracks const& tracks)
  //Belangrijk dus dat je wel subscribed hier naar aod::Zdcs om toegang te krijgen tot zdc informatie
  {
    setupCalibration();
    
    for(auto& collision : cols){
      // step 0 tm 5 A&C
      double q[6][4]; // 6 steps, each with 4 values
      double v[3]; // vx, vy, vz
      //    double v_mean[3]; // vx_mean, vy_mean, vz_mean
      
      //https://alice-notes.web.cern.ch/system/files/notes/analysis/620/017-May-31-analysis_note-ALICE_analysis_note_v2.pdf
      double ZDC_px[4] = {-1.75, 1.75,-1.75, 1.75};
      double ZDC_py[4] = {-1.75, -1.75,1.75, 1.75};
      double alpha = 0.395;
      //tot hier van die analysis note
      
      // for energy calibration
      double EZN[8]; //uncalibrated energy for the 2x4 towers (a1, a2, a3, a4, c1, c2, c3, c4)
      double meanEZN[10]; //mean energies from calibration histos (common A common C, t1-4 A, t1-4C)
      double e[8]={0., 0., 0., 0.,0., 0., 0., 0.}; // calibrated energies (a1, a2, a3, a4, c1, c2, c3, c4))
      
      //for Q-vector calculation
      // A[0] & C[1]
      double sumZN[2] = {0., 0.};
      double xEnZN[2] = {0., 0.};//"XA", "XC"
      double yEnZN[2] = {0., 0.};//"YA", "YC"
      double Q[4] = {0., 0., 0., 0.}; //"QXA", "QYA", "QXC", "QYC"
      
      
      if (!collision.sel8())
        return;
      auto cent = collision.centFT0C();
      //    LOGF(info, " centrality = %f", cent);
      if(cent<0 || cent>90)
        return;
      
      
      // TODO: check how to pick up polarity for run 2!!!
      //    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      //    auto field = (cfgMagField == 999999) ? getMagneticField(bc.timestamp()) : cfgMagField;
      
      //variable = (condition) ? expressionTrue : expressionFalse;
      //    double polarity = (field<0) ? -0.5 : 0.5;
      
      // for now just take positive polarity
      polarity = .5;
      
      const auto& foundBC = collision.foundBC_as<BCsRun3>();
      if(foundBC.has_zdc()){
        const auto& zdcCol = foundBC.zdc();
        
        // keep track of memory consumption
        if (counter % 1000 == 0) PrintMemoryConsumption();
        
        int runNumber = collision.bc().runNumber();
        
        v[0] = collision.posX();
        v[1] = collision.posY();
        v[2] = collision.posZ();
        centrality = cent;
        trackmultiplicity = tracks.size();
        runnumber = runNumber;
        
        
        for(int i=0; i<4; i++){
          EZN[i] = zdcCol.energySectorZNA()[i];
          EZN[i+4] = zdcCol.energySectorZNC()[i];
        }
        
        // Zeroth step: start energy gain equalisation
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
        if(!step0_open){ if(counter<1) LOGF(info, "cannot continue, file energy calibration not open! ...-> Start energy calibration");
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
        }
        
        if(counter<1) LOGF(info, "files for step 0 (energy Calibraton) are open!");
        
        if(!updateBinForRunNumber(runNumber, lastRunNumber, binRunNumber, numFoundint, hZN_mean[0])) return;
        
        for(int tower=2; tower<10; tower++){
          int common = (tower > 5) ? 1 : 0;
          meanEZN[tower] = hZN_mean[tower]->GetBinContent(int(binRunNumber),int(cent)+1);
          if(meanEZN[tower]>0) e[tower-2] = EZN[tower] * (0.25 * meanEZN[common]) / meanEZN[tower];
        }
        
        // Now calculate Q-vector
        for(int tower=0; tower<8; tower++){
          int side = (tower > 4) ? 1 : 0;
          int sector = tower % 4;
          sumZN[side] += TMath::Power(e[tower],alpha);
          xEnZN[side] += ZDC_px[sector] * TMath::Power(e[tower],alpha);
          yEnZN[side] += ZDC_py[sector] * TMath::Power(e[tower],alpha);
        }
        
        // "QXA", "QYA", "QXC", "QYC"
        for (int i = 0; i < 2; ++i) {
          if (sumZN[i] > 0.0) {
            Q[i*2] = xEnZN[i] / sumZN[i]; // for QXA[0] and QXC[2]
            Q[i*2+1] = yEnZN[i] / sumZN[i]; // for QYA[1] and QYC[3]
          }
        }
        
        // Here define Q-vectors after energy calibration
        for(int i=0; i<5; i++) q[0][i] = Q[i];
        
        
        if(!step1_open) {
          if (counter<1) LOGF(warning, "File for recentring does not exist | Table processed with non-recentered q-vectors!!!!");
          
          fillAllRegistries(q, v, centrality, trackmultiplicity, runnumber, polarity, 0, 1);
          
          counter++;
          lastRunNumber = runNumber;
          
          return;
        }
        
        //
        //
        //              //TODO: dit moet nog in een functie worden gezet!!
        //      // retrieve runnumber and look up in the mean TProfile2D
        //      if (runNumber != lastRunNumber_step1) {
        //        lastRunNumber_step1 = runNumber;
        //
        //        for(Int_t i=0; i<mean_1perCent_Run[0]->GetXaxis()->GetNbins(); i++){
        //          const char* label = mean_1perCent_Run[0]->GetXaxis()->GetBinLabel(i+1);
        //          //          LOGF(info, "RUN NUMBER %i | looking in bin %s", runNumber, label);
        //          int labelInt = atoi(label);
        //          if(labelInt == runNumber){bin=i+1;}
        //          //          LOGF(info, "looking in bin %i", bin);
        //        }
        //
        //        // get runnumber from bin found in previous for loop
        //        numFoundint = atoi(mean_1perCent_Run[0]->GetXaxis()->GetBinLabel(bin));
        //      }
        //
        //      // get mean value for given centrality
        //
        //      double mean1[4]; //"QXA", "QYA", "QXC", "QYC"
        //      for(int i=0; i<4, i++){
        //        mean1[i] = mean_1perCent_Run[i]->GetBinContent(int(bin),int(cent)+1);
        //        q[1][i] = q[0][i] - mean1[i];
        //      }
        //
        //      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        //
        //      if(!step2_open){
        //        if (counter<1)  LOGF(warning, "Histos needed fot step 2 do not exist. Histograms filled up to step 1");
        //        fillAllRegistries(q, v, centrality, trackmultiplicity, runnumber, polarity, 0, 2);
        //
        //        lastRunNumberE = runNumber;
        //        counter++;
        //        return; }
        //
        //      //center vx, vy, vz around 0 on a run-by-run basis
        //      for(int i=0; i<4; i++){
        //        v_mean[i] = v[i] - v_mean_Run[i]->GetBinContent(int(bin));
        //      }
        //
        //
        //      for(int i=0; i<4; i++){
        //        q[2][i] = getQstep2(q[1][i], mean_10perCent_v[i], centrality, v_mean[0], v_mean[1], v_mean[2]);
        //      }
        //
        //      if(!step3_open){
        //        if (counter<1) LOGF(warning, "outFitStep3.root not found! STEP 3 failed: Create table with only first and second step of re-centring");
        //
        //        fillAllRegistries(q, v_mean, centrality, trackmultiplicity, runnumber, polarity, 0, 3);
        //
        //        counter++;
        //        lastRunNumberE = runNumber;
        //        return;
        //      }
        //
        //      //TODO: We zijn tot hier gekomen!! 13-06-2024
        //      int bin_run_step3 = -1;
        //      const char* label;
        //      for(Int_t i=0; i<hFit_QXA_Run_Cent_vy->GetXaxis()->GetNbins(); i++){ // sum over runs
        //        label = hFit_QXA_Run_Cent_vy->GetXaxis()->GetBinLabel(i+1);
        //        int labelInt = atoi(label);
        //        if(labelInt == runNumber){bin_run_step3=i+1; }//LOGF(info, "bin for run is found and set");}
        //      }
        //
        //      int bin_cent = FindBinWithContentAlongAxis(hQXA_Run_10perCent_vy, 1, centrality);
        //
        //      //       LOGF(info, "runnumber is: %i in bin %i | centrality is: %e in bin %i", runnumber, bin_run, centrality, bin_cent); // <- checked, goes well
        //
        //      double q_mean_xa_vy_step3 = hFit_QXA_Run_Cent_vy->GetBinContent(bin_run_step3, bin_cent)*vy_mean;
        //      double q_mean_xa_vz_step3 = hFit_QXA_Run_Cent_vz->GetBinContent(bin_run_step3, bin_cent)*vz_mean;
        //
        //      double q_mean_xc_vy_step3 = hFit_QXC_Run_Cent_vy->GetBinContent(bin_run_step3, bin_cent)*vy_mean;
        //      double q_mean_xc_vz_step3 = hFit_QXC_Run_Cent_vz->GetBinContent(bin_run_step3, bin_cent)*vz_mean;
        //
        //      double q_mean_ya_vx_step3 = hFit_QYA_Run_Cent_vx->GetBinContent(bin_run_step3, bin_cent)*vx_mean;
        //      double q_mean_ya_vz_step3 = hFit_QYA_Run_Cent_vz->GetBinContent(bin_run_step3, bin_cent)*vz_mean;
        //
        //      double q_mean_yc_vx_step3 = hFit_QYC_Run_Cent_vx->GetBinContent(bin_run_step3, bin_cent)*vx_mean;
        //      double q_mean_yc_vz_step3 = hFit_QYC_Run_Cent_vz->GetBinContent(bin_run_step3, bin_cent)*vz_mean;
        //
        //      //      LOGF(info, "q_mean_xa_step3 = %f", q_mean_xa_step3);
        //      //      LOGF(info, "q_mean_xc_step3 = %f", q_mean_xc_step3);
        //      //      LOGF(info, "q_mean_ya_step3 = %f", q_mean_ya_step3);
        //      //      LOGF(info, "q_mean_yc_step3 = %f", q_mean_yc_step3);
        //
        //      qxa3 = qxa2 - q_mean_xa_vz_step3 - q_mean_xa_vy_step3;
        //      qya3 = qya2 - q_mean_ya_vz_step3 - q_mean_ya_vx_step3;
        //      qxc3 = qxc2 - q_mean_xc_vz_step3 - q_mean_xc_vy_step3;
        //      qyc3 = qyc2 - q_mean_yc_vz_step3 - q_mean_yc_vx_step3;
        //
        //
        //      if(!step4_open) {
        //        if (counter<1) LOGF(warning, "Histos needed fot step 4 do not exist. Histograms filled up to step 3");
        //
        //        fillRegistry(0, qxa0, qya0, qxc0, qyc0, vx, vy, vz, centrality, trackmultiplicity, runnumber, polarity);
        //        fillRegistry(1, qxa1, qya1, qxc1, qyc1, vx, vy, vz, centrality, trackmultiplicity, runnumber, polarity);
        //        fillRegistry(2, qxa2, qya2, qxc2, qyc2, vx_mean, vy_mean, vz_mean, centrality, trackmultiplicity, runnumber, polarity);
        //        fillRegistry(3, qxa3, qya3, qxc3, qyc3, vx_mean, vy_mean, vz_mean, centrality, trackmultiplicity, runnumber, polarity);
        //
        //        counter++;
        //        lastRunNumberE = runNumber;
        //        return; }
        //
        //      qxa4 = getQstep4(qxa3, hQXA_mean_Magnet_10perCent_v, polarity, centrality, vx_mean, vy_mean, vz_mean);
        //      qya4 = getQstep4(qya3, hQYA_mean_Magnet_10perCent_v, polarity, centrality, vx_mean, vy_mean, vz_mean);
        //      qxc4 = getQstep4(qxc3, hQXC_mean_Magnet_10perCent_v, polarity, centrality, vx_mean, vy_mean, vz_mean);
        //      qyc4 = getQstep4(qyc3, hQYC_mean_Magnet_10perCent_v, polarity, centrality, vx_mean, vy_mean, vz_mean);
        //
        //      //      LOGF(info,"trak multiplicity before = %.2f", trackmultiplicity);
        //      trackmultiplicity_mean = tracks.size() - hmeanN_1perCent_Run->GetBinContent(bin,int(cent)+1);
        //      //      LOGF(info,"mean trak multiplicity = %.2f", trackmultiplicity_mean);
        //
        //
        //      if(!step5_open){
        //        if (counter<1) LOGF(error, "Histos needed fot step 5 do not exist. Histograms filled up to step 4");
        //        fillRegistry(0, qxa0, qya0, qxc0, qyc0, vx, vy, vz, centrality, trackmultiplicity, runnumber, polarity);
        //        fillRegistry(1, qxa1, qya1, qxc1, qyc1, vx, vy, vz, centrality, trackmultiplicity, runnumber, polarity);
        //        fillRegistry(2, qxa2, qya2, qxc2, qyc2, vx_mean, vy_mean, vz_mean, centrality, trackmultiplicity, runnumber, polarity);
        //        fillRegistry(3, qxa3, qya3, qxc3, qyc3, vx_mean, vy_mean, vz_mean, centrality, trackmultiplicity, runnumber, polarity);
        //        fillRegistry(4, qxa4, qya4, qxc4, qyc4, vx_mean, vy_mean, vz_mean, centrality, trackmultiplicity_mean, runnumber, polarity);
        //        counter++;
        //        lastRunNumberE = runNumber;
        //        return;
        //      }
        //
        //      int bin_run5 = runnumber+1;
        //
        //      qxa5 = getQstep5(qxa4, hQXA_mean_run_cent10_Mult, bin_run5, centrality, trackmultiplicity_mean);
        //      qya5 = getQstep5(qya4, hQYA_mean_run_cent10_Mult, bin_run5, centrality, trackmultiplicity_mean);
        //      qxc5 = getQstep5(qxc4, hQXC_mean_run_cent10_Mult, bin_run5, centrality, trackmultiplicity_mean);
        //      qyc5 = getQstep5(qyc4, hQYC_mean_run_cent10_Mult, bin_run5, centrality, trackmultiplicity_mean);
        //
        //      if(step5_open){
        //        if (counter<1) LOGF(info, "Congrats!!! You are done. Histograms filled up to step 5");
        //
        //        fillRegistry(0, qxa0, qya0, qxc0, qyc0, vx, vy, vz, centrality, trackmultiplicity, runnumber, polarity);
        //        fillRegistry(1, qxa1, qya1, qxc1, qyc1, vx, vy, vz, centrality, trackmultiplicity, runnumber, polarity);
        //        fillRegistry(2, qxa2, qya2, qxc2, qyc2, vx_mean, vy_mean, vz_mean, centrality, trackmultiplicity, runnumber, polarity);
        //        fillRegistry(3, qxa3, qya3, qxc3, qyc3, vx_mean, vy_mean, vz_mean, centrality, trackmultiplicity, runnumber, polarity);
        //        fillRegistry(4, qxa4, qya4, qxc4, qyc4, vx_mean, vy_mean, vz_mean, centrality, trackmultiplicity_mean, runnumber, polarity);
        //        fillRegistry(5, qxa5, qya5, qxc5, qyc5, vx_mean, vy_mean, vz_mean, centrality, trackmultiplicity_mean, runnumber, polarity);
        //
        //        counter++;
        //        lastRunNumberE = runNumber;
        //        return;
        //      }
        
        
      }// end collision found ZDC
    } // end of for (col : cols)
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ZDCqvectors>(cfgc)};
}




//RUN2
//if (collision.foundZDCId() >= 0) {
//  // keep track of memory consumption
//  if (counter % 1000 == 0) PrintMemoryConsumption();
//  
//  int runNumber = collision.bc().runNumber();
//  
//  v[0] = collision.posX();
//  v[1] = collision.posY();
//  v[2] = collision.posZ();
//  centrality = cent;
//  trackmultiplicity = tracks.size();
//  runnumber = runNumber;
//  
//  
//  for(int i=0; i<4; i++){
//    EZN[i] = collision.foundZDC().energySectorZNA()[i];
//    EZN[i+4] = collision.foundZDC().energySectorZNC()[i];
//  }

