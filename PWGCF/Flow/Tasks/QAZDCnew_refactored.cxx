
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
#include <DataFormatsParameters/GRPObject.h>
#include <DataFormatsParameters/GRPMagField.h>
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

// define my.....
using myCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>>;
// run2
// soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>>;
using myTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>>;
using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;

namespace o2::analysis::myanalysistask
{

bool step0_open = false;
bool step1_open = false;
bool step2_open = false;
bool step3_open = false;
bool step4_open = false;
bool step5_open = false;
bool step6_open = false;

int binRunNumber = -1;
int numFoundint = -1;
int lastRunNumber = -1;

int counter = 0;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// here for each step open histos (if possible) and set bool for each step to true. Then in process check bool and proceed or start collecting q-vectors

// step0
std::shared_ptr<TProfile2D> ZN_Energy[10] = {{nullptr}}; // fill for weights
std::shared_ptr<TH2> hQx_vs_Qy[6] = {{nullptr}};         // fill Qx vs Qy for each step in the recentring process.

// <XX> <XY> <YX> and <YY> for step 0 and step 5
std::shared_ptr<TProfile> COORD_correlations[2][4] = {{nullptr}};

std::vector<TProfile2D*> hZN_mean(10, nullptr); // Get from calibration file
// step1
std::vector<TProfile2D*> mean_1perCent_Run(4, nullptr); // hQXA, hQYA, hQXC, hQYC
// step2
std::vector<THnSparseD*> mean_10perCent_v(4, nullptr); // hQXA, hQYA, hQXC, hQYC
std::vector<TProfile*> v_mean_Run(3, nullptr);         // hvx, hvy, hvz
// step3 [fit]
std::vector<TH2D*> hFit_Run_Cent_vx(2, nullptr);            // hQYA, hQYC
std::vector<TH2D*> hFit_Run_Cent_vy(2, nullptr);            // hQXA, hQXC
std::vector<TH2D*> hFit_Run_Cent_vz(4, nullptr);            // hQXA, hQYA, hQXC, hQYC
std::vector<THnSparseD*> hQXA_Run_10perCent_vy(1, nullptr); // for check bin step3
// step4
std::vector<TProfile2D*> hmeanN_1perCent_Run(1, nullptr);     // Only 1
std::vector<THnSparseD*> mean_Magnet_10perCent_v(4, nullptr); // hQXA, hQYA, hQXC, hQYC
// step5
std::vector<THnSparseD*> mean_run_cent10_Mult(4, nullptr); // hQXA, hQYA, hQXC, hQYC

std::vector<const char*> sides = {"A", "C"};
std::vector<const char*> coords = {"x", "y", "z"};
std::vector<const char*> COORDS = {"X", "Y"};

// https://alice-notes.web.cern.ch/system/files/notes/analysis/620/017-May-31-analysis_note-ALICE_analysis_note_v2.pdf
double ZDC_px[4] = {-1.75, 1.75, -1.75, 1.75};
double ZDC_py[4] = {-1.75, -1.75, 1.75, 1.75};
double alphaZDC = 0.395;
// tot hier van die analysis note

} // namespace o2::analysis::myanalysistask

using namespace o2::analysis::myanalysistask;

struct ZDCqvectors {
  ConfigurableAxis axisCent{"axisCent", {90, 0, 90}, "Centrality axis in 1% bins"};
  ConfigurableAxis axisCent10{"axisCent10", {9, 0, 90}, "Centrality axis in 10% bins"};
  ConfigurableAxis axisQ{"axisQ", {100, -2, 2}, "Q vector (xy) in ZDC"};
  ConfigurableAxis axisVx{"axisVx", {100, 0, 0.15}, "for Pos X of collision"};
  ConfigurableAxis axisVy{"axisVy", {100, 0.35, 0.41}, "for Pos Y of collision"};
  ConfigurableAxis axisVx_mean{"axisVx_mean", {100, -0.03, 0.03}, "for Pos X of collision"};
  ConfigurableAxis axisVy_mean{"axisVy_mean", {100, -0.03, 0.03}, "for Pos Y of collision"};
  ConfigurableAxis axisVz{"axisVz", {100, -12, 12}, "for vz of collision"}; // take 12 because we shift vi - <vi>
  ConfigurableAxis axisRun{"axisRun", {1e6, 0, 1e6}, "for runNumber in ThnSparse"};
  ConfigurableAxis axisPolarity{"axisPolarity", {2, -1, 1}, "Magnet Polarity"};
  ConfigurableAxis axisMult{"axisMult", {5000, 0, 5000}, "Track Multiplicity"};
  ConfigurableAxis axisMult_mean{"axisMult_mean", {5000, -2500, 2500}, "Track Multiplicity mean per run"};

  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMin, float, 0.2f, "Minimal.q pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMax, float, 10.0f, "Maximal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "Minimal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 3.0f, "Maximal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5, "Chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgUseNch, bool, false, "Use Nch for flow observables")
  O2_DEFINE_CONFIGURABLE(cfgNbootstrap, int, 10, "Number of subsamples")
  O2_DEFINE_CONFIGURABLE(cfgMagField, float, 99999, "Configurable magnetic field; default CCDB will be queried")
  O2_DEFINE_CONFIGURABLE(cfgQvecrecent, std::string, "AnalysisResults_out.root", "path+file from recentring")
  O2_DEFINE_CONFIGURABLE(cfgCalibenergy, std::string, "energy/AnalysisResults.root", "path+file from energy calibration")
  O2_DEFINE_CONFIGURABLE(cfgOutFitStep3, std::string, "outFitStep3.root", "path+file from energy calibration")

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

  Service<ccdb::BasicCCDBManager> ccdb;

  void init(InitContext const&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    // Tower mean energies vs. centrality ( common=t0 )
    for (int tower = 0; tower < 5; tower++) {
      ZN_Energy[tower] = registry.add<TProfile2D>(Form("Energy/hZNA_mean_t%i_cent", tower), "", kTProfile2D, {{1, 0, 1}, axisCent});
      ZN_Energy[tower + 5] = registry.add<TProfile2D>(Form("Energy/hZNC_mean_t%i_cent", tower), "", kTProfile2D, {{1, 0, 1}, axisCent});
    }

    bool hmeanStep2 = false;
    bool hmeanStep4 = false;
    bool count_step3 = false;

    // Qx_vs_Qy for each step for ZNA and ZNC
    for (int step = 0; step < 6; step++) {
      for (const char* side : sides) {
        hQx_vs_Qy[step] = registry.add<TH2>(Form("step%i/hZN%s_Qx_vs_Qy", step, side), Form("hZN%s_Qx_vs_Qy", side), kTH2F, {axisQ, axisQ});
      }
      if (step == 0 || step == 5) {
        int i = 0;
        for (const char* COORD1 : COORDS) {
          for (const char* COORD2 : COORDS) {
            // Now we get: <XX> <XY> & <YX> <YY> vs. Centrality : do only before (step0) and after (step5) of recentring.
            COORD_correlations[step % 4][i] = registry.add<TProfile>(Form("step%i/hQ%sA_Q%sC_vs_cent", step, COORD1, COORD2), Form("hQ%sA_Q%sC_vs_cent", COORD1, COORD2), kTProfile, {axisCent10});
            i++;
          }
        }
      }
      // Add histograms for each step in the calibration process.
      // Sides is {A,C} and coords is {X,Y}
      for (const char* side : sides) {
        for (const char* coord : COORDS) {
          if (step == 1) {
            registry.add(Form("step%i/hQ%s%s_mean_1perCent_Run", step, coord, side), Form("hQ%s%s_mean_1perCent_Run", coord, side), {HistType::kTProfile2D, {{1, 0., 1.}, axisCent}});
          }
          if (step == 2) {
            registry.add(Form("step%i/hQ%s%s_mean_10perCent_v", step, coord, side), Form("hQ%s%s_mean_10perCent_v", coord, side), {HistType::kTHnSparseD, {axisCent10, axisVx_mean, axisVy_mean, axisVz, axisQ}});
            if (!hmeanStep2) {
              for (const char* c : coords) {
                registry.add(Form("step%i/hv%s_mean_Run", step, c), Form("hv%s_mean_Run", c), {HistType::kTProfile, {{1, 0., 1.}}});
              }
            }
            hmeanStep2 = true;
          }
          if (step == 3) {
            registry.add(Form("step%i/hQ%s%s_mean_10perCent_Run", step, coord, side), Form("hQ%s%s_mean_10perCent_Run", coord, side), {HistType::kTProfile2D, {{1, 0., 1.}, axisCent10}});
            registry.add(Form("step%i/hQ%s%s_Run_10perCent_vx", step, coord, side), Form("hQ%s%s_Run_10perCent_vx", coord, side), {HistType::kTHnSparseD, {axisRun, axisCent10, axisVx_mean, axisQ}});
            registry.add(Form("step%i/hQ%s%s_Run_10perCent_vy", step, coord, side), Form("hQ%s%s_Run_10perCent_vy", coord, side), {HistType::kTHnSparseD, {axisRun, axisCent10, axisVy_mean, axisQ}});
            registry.add(Form("step%i/hQ%s%s_Run_10perCent_vz", step, coord, side), Form("hQ%s%s_Run_10perCent_vz", coord, side), {HistType::kTHnSparseD, {axisRun, axisCent10, axisVz, axisQ}});
          }
          if (step == 4) {
            if (!hmeanStep4) {
              registry.add(Form("step%i/hmeanN_1perCent_Run", step), "hmeanN_1perCent_Run", {HistType::kTProfile2D, {{1, 0., 1.}, axisCent}});
              hmeanStep4 = true;
            }
            registry.add(Form("step%i/hQ%s%s_mean_Magnet_10perCent_v", step, coord, side), Form("hQ%s%s_mean_Magnet_10perCent_v", coord, side), {HistType::kTHnSparseD, {axisPolarity, axisCent10, axisVx_mean, axisVy_mean, axisVz, axisQ}});
          }
          if (step == 5) {
            registry.add(Form("step%i/hQ%s%s_mean_run_cent10_Mult", step, coord, side), Form("hQ%s%s_mean_run_cent10_Mult", coord, side), {HistType::kTHnSparseD, {axisRun, axisCent10, axisMult_mean, axisQ}});
          }
        } // end of COORDS
      } // end of sides

    } // end of sum over steps

    // recentered q-vectors
    registry.add("hStep", "hStep", {HistType::kTH1D, {{6, 0., 6.}}});

    LOGF(info, "................................ ===========> SETUP CALIBRATION <=========== ................................");
    setupCalibration(); // find out what caibration histos are available and at what step to proceed.
    LOGF(info, "................................... ===========> SETUP DONE <=========== ....................................");
  }

  // qxa, qya, qxc, qyc, vx, vy, vz, centrality, trackmultiplicity, runnumber, polarity
  inline void fillRegistry(int step, double qxa, double qya, double qxc, double qyc, double vx, double vy, double vz, double vx_mean, double vy_mean, double vz_mean, double centrality, double trackmultiplicity, double trackmultiplicity_mean, int runnumber, double polarity)
  {
    if (std::isnan(qxa) || std::isnan(qya) || std::isnan(qxc) || std::isnan(qyc)) {
      //  LOGF(info, "qxa = %.4f", qxa);
      //  LOGF(info, "qya = %.4f", qya);
      //  LOGF(info, "qxc = %.4f", qxc);
      //  LOGF(info, "qyc = %.4f", qyc);
      return;
    }

    // LOGF(info, "vx = %.2f", vx);
    // LOGF(info, "vy = %.2f", vy);
    // LOGF(info, "vz = %.2f", vz);
    // LOGF(info, "centrality = %.2f", centrality);
    // LOGF(info, "trackmultiplicity = %.2f", trackmultiplicity);
    // LOGF(info, "runnumber = %i", runnumber);
    // LOGF(info, "polarity = %.2f", polarity);

    // TODO: dit kan korter door vector van histos te maken, kan dan ook in loopje worden gevuld.
    if (step == 0) {
      float hQXA_QXC = qxa * qxc;
      float hQYA_QYC = qya * qyc;
      float hQYA_QXC = qxc * qya;
      float hQXA_QYC = qyc * qxa;

      registry.get<TProfile>(HIST("step0/hQXA_QXC_vs_cent"))->Fill(centrality, hQXA_QXC);
      registry.get<TProfile>(HIST("step0/hQYA_QYC_vs_cent"))->Fill(centrality, hQYA_QYC);
      registry.get<TProfile>(HIST("step0/hQYA_QXC_vs_cent"))->Fill(centrality, hQYA_QXC);
      registry.get<TProfile>(HIST("step0/hQXA_QYC_vs_cent"))->Fill(centrality, hQXA_QYC);

      registry.fill(HIST("step0/hZNA_Qx_vs_Qy"), qxa, qya);
      registry.fill(HIST("step0/hZNC_Qx_vs_Qy"), qxc, qyc);

      registry.get<TProfile2D>(HIST("step1/hQXA_mean_1perCent_Run"))->Fill(Form("%d", runnumber), centrality, qxa, 1);
      registry.get<TProfile2D>(HIST("step1/hQXC_mean_1perCent_Run"))->Fill(Form("%d", runnumber), centrality, qxc, 1);
      registry.get<TProfile2D>(HIST("step1/hQYA_mean_1perCent_Run"))->Fill(Form("%d", runnumber), centrality, qya, 1);
      registry.get<TProfile2D>(HIST("step1/hQYC_mean_1perCent_Run"))->Fill(Form("%d", runnumber), centrality, qyc, 1);

      registry.get<TProfile>(HIST("step2/hvx_mean_Run"))->Fill(Form("%d", runnumber), vx);
      registry.get<TProfile>(HIST("step2/hvy_mean_Run"))->Fill(Form("%d", runnumber), vy);
      registry.get<TProfile>(HIST("step2/hvz_mean_Run"))->Fill(Form("%d", runnumber), vz);

      registry.fill(HIST("hStep"), .5, 1);
    }

    if (step == 1) {
      registry.fill(HIST("step1/hZNA_Qx_vs_Qy"), qxa, qya);
      registry.fill(HIST("step1/hZNC_Qx_vs_Qy"), qxc, qyc);

      registry.fill(HIST("step2/hQYA_mean_10perCent_v"), centrality, vx_mean, vy_mean, vz_mean, qya);
      registry.fill(HIST("step2/hQXC_mean_10perCent_v"), centrality, vx_mean, vy_mean, vz_mean, qxc);
      registry.fill(HIST("step2/hQYC_mean_10perCent_v"), centrality, vx_mean, vy_mean, vz_mean, qyc);
      registry.fill(HIST("step2/hQXA_mean_10perCent_v"), centrality, vx_mean, vy_mean, vz_mean, qxa);

      registry.fill(HIST("hStep"), 1.5, 1);
    }

    if (step == 2) {

      registry.fill(HIST("step2/hZNA_Qx_vs_Qy"), qxa, qya);
      registry.fill(HIST("step2/hZNC_Qx_vs_Qy"), qxc, qyc);

      // cross terms -> Denk niet vz!
      registry.fill(HIST("step3/hQXA_Run_10perCent_vy"), runnumber, centrality, vy_mean, qxa);
      registry.fill(HIST("step3/hQXA_Run_10perCent_vz"), runnumber, centrality, vz_mean, qxa);
      registry.fill(HIST("step3/hQYA_Run_10perCent_vx"), runnumber, centrality, vx_mean, qya);
      registry.fill(HIST("step3/hQYA_Run_10perCent_vz"), runnumber, centrality, vz_mean, qya);
      registry.fill(HIST("step3/hQXC_Run_10perCent_vy"), runnumber, centrality, vy_mean, qxc);
      registry.fill(HIST("step3/hQXC_Run_10perCent_vz"), runnumber, centrality, vz_mean, qxc);
      registry.fill(HIST("step3/hQYC_Run_10perCent_vx"), runnumber, centrality, vx_mean, qyc);
      registry.fill(HIST("step3/hQYC_Run_10perCent_vz"), runnumber, centrality, vz_mean, qyc);

      // fill with zeros for fitting later.
      registry.get<TProfile2D>(HIST("step3/hQXA_mean_10perCent_Run"))->Fill(Form("%d", runnumber), centrality, 0, 1);
      registry.get<TProfile2D>(HIST("step3/hQXC_mean_10perCent_Run"))->Fill(Form("%d", runnumber), centrality, 0, 1);
      registry.get<TProfile2D>(HIST("step3/hQYA_mean_10perCent_Run"))->Fill(Form("%d", runnumber), centrality, 0, 1);
      registry.get<TProfile2D>(HIST("step3/hQYC_mean_10perCent_Run"))->Fill(Form("%d", runnumber), centrality, 0, 1);

      registry.fill(HIST("hStep"), 2.5, 1);
    }

    if (step == 3) {
      registry.fill(HIST("step3/hZNA_Qx_vs_Qy"), qxa, qya);
      registry.fill(HIST("step3/hZNC_Qx_vs_Qy"), qxc, qyc);

      registry.fill(HIST("step4/hQXA_mean_Magnet_10perCent_v"), polarity, centrality, vx_mean, vy_mean, vz_mean, qxa);
      registry.fill(HIST("step4/hQYA_mean_Magnet_10perCent_v"), polarity, centrality, vx_mean, vy_mean, vz_mean, qya);
      registry.fill(HIST("step4/hQXC_mean_Magnet_10perCent_v"), polarity, centrality, vx_mean, vy_mean, vz_mean, qxc);
      registry.fill(HIST("step4/hQYC_mean_Magnet_10perCent_v"), polarity, centrality, vx_mean, vy_mean, vz_mean, qyc);

      registry.get<TProfile2D>(HIST("step4/hmeanN_1perCent_Run"))->Fill(Form("%d", runnumber), centrality, trackmultiplicity, 1);

      //  registry.get<TProfile>(HIST("step3/hQXA_Run_10perCent_vs_vx"))->Fill(vx_mean, qxa);
      //  registry.get<TProfile>(HIST("step3/hQYA_Run_10perCent_vs_vx"))->Fill(vx_mean, qya);
      //  registry.get<TProfile>(HIST("step3/hQXC_Run_10perCent_vs_vx"))->Fill(vx_mean, qxc);
      //  registry.get<TProfile>(HIST("step3/hQYC_Run_10perCent_vs_vx"))->Fill(vx_mean, qyc);

      //  registry.get<TProfile>(HIST("step3/hQXA_Run_10perCent_vs_vy"))->Fill(vy_mean, qxa);
      //  registry.get<TProfile>(HIST("step3/hQYA_Run_10perCent_vs_vy"))->Fill(vy_mean, qya);
      //  registry.get<TProfile>(HIST("step3/hQXC_Run_10perCent_vs_vy"))->Fill(vy_mean, qxc);
      //  registry.get<TProfile>(HIST("step3/hQYC_Run_10perCent_vs_vy"))->Fill(vy_mean, qyc);

      //  registry.get<TProfile>(HIST("step3/hQXA_Run_10perCent_vs_vz"))->Fill(vz_mean, qxa);
      //  registry.get<TProfile>(HIST("step3/hQYA_Run_10perCent_vs_vz"))->Fill(vz_mean, qya);
      //  registry.get<TProfile>(HIST("step3/hQXC_Run_10perCent_vs_vz"))->Fill(vz_mean, qxc);
      //  registry.get<TProfile>(HIST("step3/hQYC_Run_10perCent_vs_vz"))->Fill(vz_mean, qyc);

      registry.fill(HIST("hStep"), 3.5, 1);
    }

    if (step == 4) {

      registry.fill(HIST("step4/hZNA_Qx_vs_Qy"), qxa, qya);
      registry.fill(HIST("step4/hZNC_Qx_vs_Qy"), qxc, qyc);

      registry.fill(HIST("step5/hQXA_mean_run_cent10_Mult"), runnumber, centrality, trackmultiplicity_mean, qxa);
      registry.fill(HIST("step5/hQYA_mean_run_cent10_Mult"), runnumber, centrality, trackmultiplicity_mean, qya);
      registry.fill(HIST("step5/hQXC_mean_run_cent10_Mult"), runnumber, centrality, trackmultiplicity_mean, qxc);
      registry.fill(HIST("step5/hQYC_mean_run_cent10_Mult"), runnumber, centrality, trackmultiplicity_mean, qyc);

      registry.fill(HIST("hStep"), 4.5, 1);
    }

    // recentered q-vectors
    if (step == 5) {
      registry.fill(HIST("step5/hZNA_Qx_vs_Qy"), qxa, qya);
      registry.fill(HIST("step5/hZNC_Qx_vs_Qy"), qxc, qyc);

      registry.get<TProfile>(HIST("step5/hQXA_QXC_vs_cent"))->Fill(centrality, qxa * qxc);
      registry.get<TProfile>(HIST("step5/hQYA_QYC_vs_cent"))->Fill(centrality, qya * qyc);
      registry.get<TProfile>(HIST("step5/hQYA_QXC_vs_cent"))->Fill(centrality, qxc * qya);
      registry.get<TProfile>(HIST("step5/hQXA_QYC_vs_cent"))->Fill(centrality, qyc * qxa);

      registry.fill(HIST("hStep"), 5.5, 1);
    }
  }

  // Function to fill the registry
  void fillAllRegistries(double q[6][4], double v[3], double v_mean[3], double cent, double trackmult, double trackmult_mean, int runnum, int polarity, int start = 0, int end = 6)
  {
    for (int i = start; i < end; ++i) {
      // LOGF(info, "filling step %i", i); 
      fillRegistry(i, q[i][0], q[i][1], q[i][2], q[i][3], v[0], v[1], v[2], v_mean[0], v_mean[1], v_mean[2], cent, trackmult, trackmult_mean, runnum, polarity);
    }
  }

  template <typename T>
  bool retrieveAndCheck(std::vector<T*>& vec,
                        const char* filename,
                        const char* fullPath,
                        std::vector<const char*> names,
                        const char* histname,
                        bool closeFile = false)
  {

    bool file_open = true;
    TFile* file = nullptr;

    // TODO: Find file in CCDB and take TList. In this TList is exactly the same as in file in TDir z-d-cqvectors/ !!
    // CHANGE THIS FUNCTION FOR CCDB 
    
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
      // LOGF(info, "We opened file %s", file->GetName());
      for (int i = 0; i < names.size(); ++i) {
        vec[i] = (T*)file->Get(Form("%s%s%s", fullPath, names[i], histname));
        if (vec[i] == nullptr || vec[i]->GetEntries() < 1) {
          LOGF(info, "%s%s%s not found or empty! Produce calibration file at given step", fullPath, names[i], histname);
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
    std::vector<const char*> names_hQ = {"QXA", "QYA", "QXC", "QYC"};
    std::vector<const char*> names_vcoord = {"vx", "vy", "vz"};
    std::vector<const char*> ZN_towers = {"ZNA_mean_t0", "ZNC_mean_t0", "ZNA_mean_t1", "ZNA_mean_t2", "ZNA_mean_t3", "ZNA_mean_t4", "ZNC_mean_t1", "ZNC_mean_t2", "ZNC_mean_t3", "ZNC_mean_t4"};

    // for(auto name : ZN_towers) LOGF(info, "%s", name);

    if (retrieveAndCheck(hZN_mean, cfgQvecrecent->c_str(), "z-d-cqvectors/Energy/h", ZN_towers, "_cent")) {
      step0_open = true;
    } else {
      LOGF(info, "we return at step 0");
      return;
    }

    if (retrieveAndCheck(mean_1perCent_Run, cfgQvecrecent->c_str(), "z-d-cqvectors/step1/h", names_hQ, "_mean_1perCent_Run", true) &&
        retrieveAndCheck(v_mean_Run, cfgQvecrecent->c_str(), "z-d-cqvectors/step2/h", names_vcoord, "_mean_Run")) {
      step1_open = true;
    } else {
      LOGF(info, "we return at step 1");
      return;
    }

    if (retrieveAndCheck(mean_10perCent_v, cfgQvecrecent->c_str(), "z-d-cqvectors/step2/h", names_hQ, "_mean_10perCent_v")) {
      step2_open = true;
    } else {
      LOGF(info, "we return at step 2");
      return;
    }

    if (retrieveAndCheck(hFit_Run_Cent_vx, cfgOutFitStep3->c_str(), "hFit_", {names_hQ[1], names_hQ[3]}, "_Run_Cent_vx") && retrieveAndCheck(hFit_Run_Cent_vy, cfgOutFitStep3->c_str(), "hFit_", {names_hQ[0], names_hQ[2]}, "_Run_Cent_vy") && retrieveAndCheck(hFit_Run_Cent_vz, cfgOutFitStep3->c_str(), "hFit_", names_hQ, "_Run_Cent_vz", true) // close cfgOutFitStep3
        && retrieveAndCheck(hQXA_Run_10perCent_vy, cfgQvecrecent->c_str(), "z-d-cqvectors/step3/h", {names_hQ[0]}, "_Run_10perCent_vy")) {
      step3_open = true;
    } else
      return;

    if (retrieveAndCheck(hmeanN_1perCent_Run, cfgQvecrecent->c_str(), "z-d-cqvectors/step4/", {""}, "hmeanN_1perCent_Run") && retrieveAndCheck(mean_Magnet_10perCent_v, cfgQvecrecent->c_str(), "z-d-cqvectors/step4/h", names_hQ, "_mean_Magnet_10perCent_v")) {
      step4_open = true;
    } else
      return;

    if (retrieveAndCheck(mean_run_cent10_Mult, cfgQvecrecent->c_str(), "z-d-cqvectors/step5/h", names_hQ, "_mean_run_cent10_Mult", true)) {
      step5_open = true;
    } else
      return;

    if (step0_open)
      LOGF(info, "histos energy calibration open! ");
    if (step1_open)
      LOGF(info, "histos step 1 open! ");
    if (step2_open)
      LOGF(info, "histos step 2 open! ");
    if (step3_open)
      LOGF(info, "histos step 3 open! ");
    if (step4_open)
      LOGF(info, "histos step 4 open! ");
    if (step5_open)
      LOGF(info, "histos step 5 open! ");
  }

  int FindBinWithContentAlongAxis(THnSparse* hist, Int_t axis, Double_t targetContent)
  {
    // to find bin for given axis in THnSparse object
    int bin = 1;

    TAxis* ax = hist->GetAxis(axis);
    // Check if targetContent is within the range of the axis
    Double_t min = ax->GetXmin();
    Double_t max = ax->GetXmax();

    if (targetContent < min || targetContent > max) {
      LOGF(warning, "targetContent (%.2f) is out of range of axis (%.2f,%.2f)", targetContent, min, max);
      return 0;
    }

    bin = ax->FindBin(targetContent);

    return bin;
  }

  int getMagneticField(uint64_t timestamp)
  {
    // TODO done only once (and not per run). Will be replaced by CCDBConfigurable
    // static o2::parameters::GRPObject* grpo = nullptr;
    static o2::parameters::GRPMagField* grpo = nullptr;
    if (grpo == nullptr) {
      // grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", timestamp);
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    return grpo->getNominalL3Field();
  }

  double recenterQfromTHn(double qBef, THnSparseD* meanObject, std::vector<double> variables)
  {
    // Make sure the meanObjects have the Axis filled in the same order as variables with the last axis containing the Q-vectors!
    double mean_q = 0.;
    std::vector<int> bins;
    // LOGF(info, "entries before: %.2f ",meanObject->GetEntries());
    for (int i = 0; i < variables.size(); i++) {
      bins.push_back(FindBinWithContentAlongAxis(meanObject, i, variables[i]));
      meanObject->GetAxis(i)->SetRange(bins[i], bins[i]);
      // LOGF(info, "set range to bin %i ", bins[i] );
    }

    // LOGF(info, "entries after: %.2f ",meanObject->Projection(variables.size())->GetEntries() );
    if (meanObject->Projection(variables.size())->GetEntries() > 1) {
      mean_q = meanObject->Projection(variables.size())->GetMean();
    }

    double qAft = qBef - mean_q;

    return qAft;
  }

  // From step4:
  // double pol, double cent, double vxmean, double vymean, double vzmean)

  //  double getQstep2(double q_step1, THnSparseD* mean, int cent, double vxmean, double vymean, double vzmean)
  //  {
  //    double mean_q =0.;
  //
  //
  //    // do we need to find the bins of v and i
  //    int bin_cent =   FindBinWithContentAlongAxis(mean, 0, cent);
  //    int bin_vxmean = FindBinWithContentAlongAxis(mean, 1, vxmean);
  //    int bin_vymean = FindBinWithContentAlongAxis(mean, 2, vymean);
  //    int bin_vzmean = FindBinWithContentAlongAxis(mean, 3, vzmean);
  //
  //    mean->GetAxis(0)->SetRange(bin_cent, bin_cent);
  //    mean->GetAxis(1)->SetRange(bin_vxmean, bin_vxmean);
  //    mean->GetAxis(2)->SetRange(bin_vymean, bin_vymean);
  //    mean->GetAxis(3)->SetRange(bin_vzmean, bin_vzmean);
  //
  //    mean_q = mean->Projection(4)->GetMean();
  //
  //    double q_step2 = q_step1 - mean_q;
  //
  //    return q_step2;
  //  }

  //  double getQstep5(double q_step4, THnSparseD* mean, double runN, double cent, double mult)
  //  {
  //    //
  //
  //    double mean_q =0.;
  //
  //    int bin_runmbr = FindBinWithContentAlongAxis(mean, 0, runN);
  //    int bin_cent =   FindBinWithContentAlongAxis(mean, 1, cent);
  //    int bin_mult =   FindBinWithContentAlongAxis(mean, 2, mult);
  //
  //    mean->GetAxis(0)->SetRange(bin_runnmbr, bin_runnmbr);
  //    mean->GetAxis(1)->SetRange(bin_cent, bin_cent);
  //    mean->GetAxis(2)->SetRange(bin_mult, bin_mult);
  //
  //
  //    mean_q = mean->Projection(3)->GetMean();
  //
  //    double q_step5 = q_step4 - mean_q;
  //
  //    //    delete proQ;
  //
  //    return q_step5;
  //  }

  inline void PrintMemoryConsumption()
  {
    ProcInfo_t procInfo;
    gSystem->GetProcInfo(&procInfo);

    std::cout << "Memory Usage:" << std::endl;
    std::cout << "  Resident memory: " << (double)procInfo.fMemResident / (1e6) << " MB" << std::endl;
    std::cout << "  Virtual memory: " << (double)procInfo.fMemVirtual / (1e6) << " MB" << std::endl;
  }

  template <typename T>
  bool updateBinForRunNumber(int runNumber, int& lastRunNumber, int& binRunNumber, int& label_int, T* hist)
  // Function that returns true if new runnumber is used and match is found in
  {
    bool found = false;

    if (runNumber != lastRunNumber) {
      // lastRunNumber = runNumber; // dit is alleen handig als je het maar 1x doet

      for (int i = 0; i < hist->GetXaxis()->GetNbins(); i++) {
        const char* label = hist->GetXaxis()->GetBinLabel(i + 1);
        label_int = atoi(label);
        if (label_int == runNumber) {
          binRunNumber = i + 1;
          found = true;
        }
        if (found)
          break; // exit the for loop if the right bin is found!
      }
    }

    // do this to make sure function returns false if runnumber not found.
    if (label_int != runNumber) {
      LOGF(warning, "No match found for RUN NUMBER %i", runNumber);
      found = false;
    } else {
      found = true;
    }

    return found;
  }

  // TODO: void: make function for step 1! Also for step 2,3,4,5 And apply function in processRun2 and ProcessRun3

  void process(myCollisions::iterator const& collision,
               BCsRun3 const& /*bcs*/,
               aod::Zdcs const& /*zdcs*/,
               myTracks const& tracks)
  // Belangrijk dus dat je wel subscribed hier naar aod::Zdcs om toegang te krijgen tot zdc informatie
  {
    // step 0 tm 5 A&C
    double q[6][4];   // 6 steps, each with 4 values
    double v[3];      // vx, vy, vz
    double v_mean[3]; // vx_mean, vy_mean, vz_mean

    // for energy calibration
    double EZN[8];                                  // uncalibrated energy for the 2x4 towers (a1, a2, a3, a4, c1, c2, c3, c4)
    double meanEZN[10];                             // mean energies from calibration histos (common A common C, t1-4 A, t1-4C)
    double e[8] = {0., 0., 0., 0., 0., 0., 0., 0.}; // calibrated energies (a1, a2, a3, a4, c1, c2, c3, c4))

    // for Q-vector calculation
    //  A[0] & C[1]
    double sumZN[2] = {0., 0.};
    double xEnZN[2] = {0., 0.};     //"XA", "XC"
    double yEnZN[2] = {0., 0.};     //"YA", "YC"
    double Q[4] = {0., 0., 0., 0.}; //"QXA", "QYA", "QXC", "QYC"

    if (!collision.sel8())
      return;
    auto cent = collision.centFT0C();
    if (cent < 0 || cent > 90)
      return;

    const auto& foundBC = collision.foundBC_as<BCsRun3>();
    auto field = (cfgMagField == 999999) ? getMagneticField(foundBC.timestamp()) : cfgMagField;
    if (foundBC.has_zdc()) {

      if (counter % 1000 == 0)
        PrintMemoryConsumption();

      v[0] = collision.posX();
      v[1] = collision.posY();
      v[2] = collision.posZ();
      centrality = cent;
      trackmultiplicity = tracks.size();
      polarity = (field < 0) ? -.5 : .5;

      int runnumber = foundBC.runNumber();
      const auto& zdcCol = foundBC.zdc();

      // Get the raw energies EZN[8] (not the mean A,C)
      for (int tower = 0; tower < 4; tower++) {
        EZN[tower] = zdcCol.energySectorZNA()[tower];
        EZN[tower + 4] = zdcCol.energySectorZNC()[tower];
      }

      if (!step0_open) {
        if (counter < 1) {
          LOGF(info, "files for step 0 not found ");
        }
      }
      if (counter < 1) {
        LOGF(info, "=====================> .....Start Energy Calibration..... <=====================");
      }
      // Always create the mean energy per tower histos !
      for (int tower = 0; tower < 5; tower++) {
        if (tower == 0) {
          if (zdcCol.energyCommonZNA() > 0)
            ZN_Energy[tower]->Fill(Form("%d", runnumber), cent, zdcCol.energyCommonZNA(), 1);
          if (zdcCol.energyCommonZNC() > 0)
            ZN_Energy[tower + 5]->Fill(Form("%d", runnumber), cent, zdcCol.energyCommonZNC(), 1);
          // LOGF(info, "Common A tower filled with: %i, %.2f, %.2f", runnumber, cent,zdcCol.energyCommonZNA() );
        } else {
          if (zdcCol.energySectorZNA()[tower - 1] > 0)
            ZN_Energy[tower]->Fill(Form("%d", runnumber), cent, zdcCol.energySectorZNA()[tower - 1], 1);
          if (zdcCol.energySectorZNC()[tower - 1] > 0)
            ZN_Energy[tower + 5]->Fill(Form("%d", runnumber), cent, zdcCol.energySectorZNC()[tower - 1], 1);
          // LOGF(info, "Tower ZNC[%i] filled with: %i, %.2f, %.2f", tower, runnumber, cent,zdcCol.energySectorZNC()[tower-1]);
        }
      }

      if (step0_open) {
        // Only go here if step_0 is open! (Bool set to true or false in setupCalibration())
        if (counter < 1)
          LOGF(info, "files for step 0 (energy Calibraton) are open!");

        if (!updateBinForRunNumber(runnumber, lastRunNumber, binRunNumber, numFoundint, hZN_mean[0])) {
          LOGF(warning, "Run Energy Calibration for run %i", runnumber);
          return;
        }

        if (counter < 1) {
          LOGF(info, "=====================> .....Start Calculating Q-Vectors..... <=====================");
        }

        // Now start gain equalisation!
        for (int tower = 0; tower < 10; tower++) {
          int common = (tower > 5) ? 1 : 0; // 0 and 1st element of meanEZN are the common channels
          if (hZN_mean[tower] == nullptr)
            LOGF(error, "hZN_meam[%i] out of range.. ", tower);
          meanEZN[tower] = hZN_mean[tower]->GetBinContent(int(binRunNumber), int(cent) + 1);
          if (meanEZN[tower] > 0 && tower >= 2) {
            e[tower - 2] = EZN[tower - 2] * (0.25 * meanEZN[common]) / meanEZN[tower];
            // LOGF(info,"e[%i] = %.2f", tower-2, e[tower-2]);  //Dit gaat goed
          }
        }

        // Now calculate Q-vector
        for (int tower = 0; tower < 8; tower++) {
          int side = (tower > 4) ? 1 : 0;
          int sector = tower % 4;
          sumZN[side] += TMath::Power(e[tower], alphaZDC);
          xEnZN[side] += ZDC_px[sector] * TMath::Power(e[tower], alphaZDC);
          yEnZN[side] += ZDC_py[sector] * TMath::Power(e[tower], alphaZDC);
        }

        // "QXA", "QYA", "QXC", "QYC"
        for (int i = 0; i < 2; ++i) {
          if (sumZN[i] > 0) {
            Q[i * 2] = xEnZN[i] / sumZN[i];     // for QXA[0] and QXC[2]
            Q[i * 2 + 1] = yEnZN[i] / sumZN[i]; // for QYA[1] and QYC[3]
          }
        }

        // Here define Q-vectors after energy calibration
        for (int i = 0; i < 5; i++) {
          q[0][i] = Q[i];
        }

        if (!step1_open) {
          if (counter < 1)
            LOGF(warning, "File for recentring does not exist | Table processed with non-recentered q-vectors!!!!");

          fillAllRegistries(q, v, v_mean, centrality, trackmultiplicity, trackmultiplicity_mean, runnumber, polarity, 0, 1);

          counter++;
          lastRunNumber = runnumber;

          return;
        }
      } // end of step0_open

      if (step1_open) {
        // TODO: check where to change runnumber to last runnumbet -> Think all the way at the bottom of process and after all registry fills!
        if (!updateBinForRunNumber(runnumber, lastRunNumber, binRunNumber, numFoundint, mean_1perCent_Run[0])) {
          LOGF(warning, "Run Recentring step 1 for run %i first!", runnumber);
          return;
        }

        // get mean value for given centrality
        double Q_mean_step1[4]; //"QXA", "QYA", "QXC", "QYC"
        for (int i = 0; i < 4; i++) {
          Q_mean_step1[i] = mean_1perCent_Run[i]->GetBinContent(int(binRunNumber), int(cent) + 1);
          // LOGF(info, "Run %s found in bin %i for centrality %.2f", mean_1perCent_Run[i]->GetXaxis()->GetBinLabel(binRunNumber), binRunNumber, cent);
          q[1][i] = q[0][i] - Q_mean_step1[i];
        }

        // center vx, vy, vz around 0 on a run-by-run basis
        // LOGF(info, "START calculating mean_v");
        if (!updateBinForRunNumber(runnumber, lastRunNumber, binRunNumber, numFoundint, v_mean_Run[0])) {
          LOGF(warning, "Run Recentring step 1 for run %i first!", runnumber);
          return;
        }

        for (int i = 0; i < 3; i++) {
          // LOGF(info, "binRunNumber = %i | for bin %s | with content %.5f", binRunNumber,v_mean_Run[i]->GetXaxis()->GetBinLabel(binRunNumber), v_mean_Run[i]->GetBinContent(binRunNumber));
          v_mean[i] = v[i] - v_mean_Run[i]->GetBinContent(binRunNumber);
          //  LOGF(info, "mean v[%i] = %.5f - %.5f", i, v[i],v_mean_Run[i]->GetBinContent(binRunNumber));
        }

        if (!step2_open) {
          if (counter < 1)
            LOGF(warning, "Histos needed fot step 2 do not exist. Histograms filled up to step 1");
          fillAllRegistries(q, v, v_mean, centrality, trackmultiplicity, trackmultiplicity_mean, runnumber, polarity, 0, 2);

          lastRunNumber = runnumber;
          counter++;
          return;
        }
      } // end of step1_open

      if (step2_open) {

        for (int i = 0; i < 4; i++) {
          q[2][i] = recenterQfromTHn(q[1][i], mean_10perCent_v[i], {centrality, v_mean[0], v_mean[1], v_mean[2]});
        }

        if (!step3_open) {
          if (counter < 1)
            LOGF(warning, "outFitStep3.root not found! STEP 3 failed: Create table with only first and second step of re-centring");

          fillAllRegistries(q, v, v_mean, centrality, trackmultiplicity, trackmultiplicity_mean, runnumber, polarity, 0, 3);

          counter++;
          lastRunNumber = runnumber;
          return;
        }
      } // end of step2_open

      if (step3_open) {
        // TODO: We zijn tot hier gekomen!! 13-06-2024
        //  verder op 22-08-2024 Alles hierbovem werkt nu! (Getest op 1 AO2D file)
        //  int bin_run_step3 = -1;
        //  const char* label;
        //  for(Int_t i=0; i<hFit_QXA_Run_Cent_vy->GetXaxis()->GetNbins(); i++){ // sum over runs
        //    label = hFit_QXA_Run_Cent_vy->GetXaxis()->GetBinLabel(i+1);
        //    int labelInt = atoi(label);
        //    if(labelInt == runNumber){bin_run_step3=i+1; }//LOGF(info, "bin for run is found and set");}
        //  }

        if (!updateBinForRunNumber(runnumber, lastRunNumber, binRunNumber, numFoundint, hFit_Run_Cent_vy[0])) {
          LOGF(warning, "Run Recentring step 3 for run %i first!", runnumber);
          return;
        }

        int bin_cent = FindBinWithContentAlongAxis(hQXA_Run_10perCent_vy[0], 1, centrality);

        // LOGF(info, "runnumber is: %i in bin %i | centrality is: %.2f in bin %i", runnumber, binRunNumber, centrality, bin_cent); // <- checked, goes well

        double q_mean_xa_vy_step3 = hFit_Run_Cent_vy[0]->GetBinContent(binRunNumber, bin_cent) * v_mean[1];
        double q_mean_xa_vz_step3 = hFit_Run_Cent_vz[0]->GetBinContent(binRunNumber, bin_cent) * v_mean[2];

        double q_mean_xc_vy_step3 = hFit_Run_Cent_vy[1]->GetBinContent(binRunNumber, bin_cent) * v_mean[1];
        double q_mean_xc_vz_step3 = hFit_Run_Cent_vz[2]->GetBinContent(binRunNumber, bin_cent) * v_mean[2];

        double q_mean_ya_vx_step3 = hFit_Run_Cent_vx[0]->GetBinContent(binRunNumber, bin_cent) * v_mean[0];
        double q_mean_ya_vz_step3 = hFit_Run_Cent_vz[1]->GetBinContent(binRunNumber, bin_cent) * v_mean[2];

        double q_mean_yc_vx_step3 = hFit_Run_Cent_vx[1]->GetBinContent(binRunNumber, bin_cent) * v_mean[0];
        double q_mean_yc_vz_step3 = hFit_Run_Cent_vz[3]->GetBinContent(binRunNumber, bin_cent) * v_mean[2];

        // LOGF(info, "q_mean_xc_vy_step3 = %.3f", q_mean_xc_vy_step3);
        // LOGF(info, "q_mean_xa_vy_step3 = %.3f", q_mean_xa_vy_step3);
        // LOGF(info, "q_mean_ya_vx_step3 = %.3f", q_mean_ya_vx_step3);
        // LOGF(info, "q_mean_yc_vx_step3 = %.3f", q_mean_yc_vx_step3);

        q[3][0] = q[2][0] - q_mean_xa_vz_step3 - q_mean_xa_vy_step3;
        q[3][1] = q[2][1] - q_mean_ya_vz_step3 - q_mean_ya_vx_step3;
        q[3][2] = q[2][2] - q_mean_xc_vz_step3 - q_mean_xc_vy_step3;
        q[3][3] = q[2][3] - q_mean_yc_vz_step3 - q_mean_yc_vx_step3;

        // LOGF(info, "q[2][0] = q[1][0] - q_mean_xa_vz_step3 - q_mean_xa_vy_step3 || %.2f = %.2f - %.2f - %.2f", q[2][0], q[1][0], q_mean_xa_vz_step3, q_mean_xa_vy_step3);
        // LOGF(info, "q[2][1] = q[1][1] - q_mean_ya_vz_step3 - q_mean_ya_vx_step3 || %.2f = %.2f - %.2f - %.2f", q[2][1], q[1][1], q_mean_ya_vz_step3, q_mean_ya_vx_step3);
        // LOGF(info, "q[2][2] = q[1][2] - q_mean_xc_vz_step3 - q_mean_xc_vy_step3 || %.2f = %.2f - %.2f - %.2f", q[2][2], q[1][2], q_mean_xc_vz_step3, q_mean_xc_vy_step3);
        // LOGF(info, "q[2][3] = q[1][3] - q_mean_yc_vz_step3 - q_mean_yc_vx_step3 || %.2f = %.2f - %.2f - %.2f", q[2][3], q[1][3], q_mean_yc_vz_step3, q_mean_yc_vx_step3);

        if (!step4_open) {
          if (counter < 1)
            LOGF(warning, "Histos needed fot step 4 do not exist. Histograms filled up to step 3");

          fillAllRegistries(q, v, v_mean, centrality, trackmultiplicity, trackmultiplicity_mean, runnumber, polarity, 0, 4);

          counter++;
          lastRunNumber = runnumber;
          return;
        }
      }
 
      if (step4_open) {
        trackmultiplicity_mean = trackmultiplicity - hmeanN_1perCent_Run[0]->GetBinContent(binRunNumber, int(centrality) + 1);
        if (trackmultiplicity == trackmultiplicity_mean) {
          LOGF(info, "track multiplicity = %.2f for centrality %.1f", trackmultiplicity, centrality);
          LOGF(info, "mean track multiplicity = %.2f", trackmultiplicity_mean);
        }

        for (int i = 0; i < 4; i++) {
          q[4][i] = recenterQfromTHn(q[3][i], mean_Magnet_10perCent_v[i], {polarity, centrality, v_mean[0], v_mean[1], v_mean[2]});
        }
        if (!step5_open) {
          if (counter < 1)
            LOGF(warning, "Histos for step 5 not found! Create table with step 1-4 of re-centring");

          fillAllRegistries(q, v, v_mean, centrality, trackmultiplicity, trackmultiplicity_mean, runnumber, polarity, 0, 5);

          counter++;
          lastRunNumber = runnumber;
          return;
        }

      } // end of step4_open
      //
      //      qxa4 = getQstep4(qxa3, hQXA_mean_Magnet_10perCent_v, polarity, centrality, vx_mean, vy_mean, vz_mean);
      //      qya4 = getQstep4(qya3, hQYA_mean_Magnet_10perCent_v, polarity, centrality, vx_mean, vy_mean, vz_mean);
      //      qxc4 = getQstep4(qxc3, hQXC_mean_Magnet_10perCent_v, polarity, centrality, vx_mean, vy_mean, vz_mean);
      //      qyc4 = getQstep4(qyc3, hQYC_mean_Magnet_10perCent_v, polarity, centrality, vx_mean, vy_mean, vz_mean);
      //
      //      //      LOGF(info,"trak multiplicity step0 = %.2f", trackmultiplicity);
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

      if (step5_open) {
        for (int i = 0; i < 4; i++) {
          q[5][i] = recenterQfromTHn(q[4][i], mean_run_cent10_Mult[i], {double(runnumber), centrality, trackmultiplicity_mean});
        }

       if(!step6_open){
        if (counter < 1) {  LOGF(info, "Congrats!!! You are done. Histograms filled up to step 5");} 

        fillAllRegistries(q, v, v_mean, centrality, trackmultiplicity, trackmultiplicity_mean, runnumber, polarity, 0, 6);
        counter++;
        lastRunNumber = runnumber;
        return;
       }
      } // end of step5_open

    } // end collision found ZDC
    counter++;
    lastRunNumber = runnumber;
  } // end of process
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ZDCqvectors>(cfgc)};
}

// RUN2
// if (collision.foundZDCId() >= 0) {
//   // keep track of memory consumption
//   if (counter % 1000 == 0) PrintMemoryConsumption();
//
//   int runNumber = collision.bc().runNumber();
//
//   v[0] = collision.posX();
//   v[1] = collision.posY();
//   v[2] = collision.posZ();
//   centrality = cent;
//   trackmultiplicity = tracks.size();
//   runnumber = runNumber;
//
//
//   for(int i=0; i<4; i++){
//     EZN[i] = collision.foundZDC().energySectorZNA()[i];
//     EZN[i+4] = collision.foundZDC().energySectorZNC()[i];
//   }

/// In the init was:

//      //step1
// DONE
//      registry.add("step1/hQXA_mean_1perCent_Run","hQXA_mean_1perCent_Run", {HistType::kTProfile2D, {{1,0.,1.} ,axisCent}});
//      registry.add("step1/hQYA_mean_1perCent_Run","hQYA_mean_1perCent_Run", {HistType::kTProfile2D, {{1,0.,1.} ,axisCent}});
//      registry.add("step1/hQXC_mean_1perCent_Run","hQXC_mean_1perCent_Run", {HistType::kTProfile2D, {{1,0.,1.} ,axisCent}});
//      registry.add("step1/hQYC_mean_1perCent_Run","hQYC_mean_1perCent_Run", {HistType::kTProfile2D, {{1,0.,1.} ,axisCent}});
//
//      //step2
//   DONE
//      registry.add("step2/hQXA_mean_10perCent_v","hQXA_mean_10perCent_v", {HistType::kTHnSparseD, {axisCent10, axisVx_mean, axisVy_mean, axisVz, axisQ}});
//      registry.add("step2/hQYA_mean_10perCent_v","hQYA_mean_10perCent_v", {HistType::kTHnSparseD, {axisCent10, axisVx_mean, axisVy_mean, axisVz, axisQ}});
//      registry.add("step2/hQXC_mean_10perCent_v","hQXC_mean_10perCent_v", {HistType::kTHnSparseD, {axisCent10, axisVx_mean, axisVy_mean, axisVz, axisQ}});
//      registry.add("step2/hQYC_mean_10perCent_v","hQYC_mean_10perCent_v", {HistType::kTHnSparseD, {axisCent10, axisVx_mean, axisVy_mean, axisVz, axisQ}});
//   DONE
//      registry.add("step2/hvx_mean_Run","hvx_mean_Run", {HistType::kTProfile, {{1,0.,1.}}});
//      registry.add("step2/hvy_mean_Run","hvy_mean_Run", {HistType::kTProfile, {{1,0.,1.}}});
//      registry.add("step2/hvz_mean_Run","hvz_mean_Run", {HistType::kTProfile, {{1,0.,1.}}});
//
//      //step 3
//   DONE
//      registry.add("step3/hQXA_mean_10perCent_Run","hQXA_mean_10perCent_Run", {HistType::kTProfile2D, {{1,0.,1.} ,axisCent10}});
//      registry.add("step3/hQYA_mean_10perCent_Run","hQYA_mean_10perCent_Run", {HistType::kTProfile2D, {{1,0.,1.} ,axisCent10}});
//      registry.add("step3/hQXC_mean_10perCent_Run","hQXC_mean_10perCent_Run", {HistType::kTProfile2D, {{1,0.,1.} ,axisCent10}});
//      registry.add("step3/hQYC_mean_10perCent_Run","hQYC_mean_10perCent_Run", {HistType::kTProfile2D, {{1,0.,1.} ,axisCent10}});
//
//      // cross terms
//    NOT DONE
//      registry.add("step3/hQXA_Run_10perCent_vy","hQXA_Run_10perCent_vy", {HistType::kTHnSparseD, {axisRun,axisCent10,axisVy_mean, axisQ}});
//      registry.add("step3/hQXA_Run_10perCent_vz","hQXA_Run_10perCent_vz", {HistType::kTHnSparseD, {axisRun,axisCent10,axisVz, axisQ}});
//
//      registry.add("step3/hQYA_Run_10perCent_vx","hQYA_Run_10perCent_vx", {HistType::kTHnSparseD, {axisRun,axisCent10,axisVx_mean, axisQ}});
//      registry.add("step3/hQYA_Run_10perCent_vz","hQYA_Run_10perCent_vz", {HistType::kTHnSparseD, {axisRun,axisCent10,axisVz, axisQ}});
//
//      registry.add("step3/hQXC_Run_10perCent_vy","hQXC_Run_10perCent_vx", {HistType::kTHnSparseD, {axisRun,axisCent10,axisVy_mean, axisQ}});
//      registry.add("step3/hQXC_Run_10perCent_vz","hQXC_Run_10perCent_vz", {HistType::kTHnSparseD, {axisRun,axisCent10,axisVz, axisQ}});
//
//      registry.add("step3/hQYC_Run_10perCent_vx","hQYC_Run_10perCent_vx", {HistType::kTHnSparseD, {axisRun,axisCent10,axisVx_mean, axisQ}});
//      registry.add("step3/hQYC_Run_10perCent_vz","hQYC_Run_10perCent_vz", {HistType::kTHnSparseD, {axisRun,axisCent10,axisVz, axisQ}});
//
//
//
//      //step 4
//    DONE
//      registry.add("step4/hmeanN_1perCent_Run","hmeanN_1perCent_Run", {HistType::kTProfile2D, {{1,0.,1.} ,axisCent}});
//   DONE
//      registry.add("step4/hQXA_mean_Magnet_10perCent_v","hQXA_mean_Magnet_10perCent_v", {HistType::kTHnSparseD, {axisPolarity, axisCent10, axisVx_mean, axisVy_mean, axisVz, axisQ}});
//      registry.add("step4/hQYA_mean_Magnet_10perCent_v","hQYA_mean_Magnet_10perCent_v", {HistType::kTHnSparseD, {axisPolarity, axisCent10, axisVx_mean, axisVy_mean, axisVz, axisQ}});
//      registry.add("step4/hQXC_mean_Magnet_10perCent_v","hQXC_mean_Magnet_10perCent_v", {HistType::kTHnSparseD, {axisPolarity, axisCent10, axisVx_mean, axisVy_mean, axisVz, axisQ}});
//      registry.add("step4/hQYC_mean_Magnet_10perCent_v","hQYC_mean_Magnet_10perCent_v", {HistType::kTHnSparseD, {axisPolarity, axisCent10, axisVx_mean, axisVy_mean, axisVz, axisQ}});
//
//      //step 5
//    DONE
//      registry.add("step5/hQXA_mean_run_cent10_Mult","hQXA_mean_run_cent10_Mult", {HistType::kTHnSparseD, {axisRun, axisCent10, axisMult_mean, axisQ}});
//      registry.add("step5/hQYA_mean_run_cent10_Mult","hQYA_mean_run_cent10_Mult", {HistType::kTHnSparseD, {axisRun, axisCent10, axisMult_mean, axisQ}});
//      registry.add("step5/hQXC_mean_run_cent10_Mult","hQXC_mean_run_cent10_Mult", {HistType::kTHnSparseD, {axisRun, axisCent10, axisMult_mean, axisQ}});
//      registry.add("step5/hQYC_mean_run_cent10_Mult","hQYC_mean_run_cent10_Mult", {HistType::kTHnSparseD, {axisRun, axisCent10, axisMult_mean, axisQ}});
//