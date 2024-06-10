
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
#include "TF1.h"
#include "TFitResult.h"
#include "TCanvas.h"
#include <iostream>
#include "TSystem.h"

// TODO: skipp steps if not needed anymore! So do not calculate qvec_step1 if we are at step 3!!
//using namespace o2::analysis;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

namespace o2::aod{
namespace myTable
{
DECLARE_SOA_COLUMN(QxA, qxa, float);//, "q-vector x component ZNA");
DECLARE_SOA_COLUMN(QyA, qya, float);//, "q-vector y component ZNA");
DECLARE_SOA_COLUMN(QxC, qxc, float);//, "q-vector x component ZNC");
DECLARE_SOA_COLUMN(QyC, qyc, float);//, "q-vector y component ZNC");
DECLARE_SOA_COLUMN(Vx, vx, float);//, "collision vertex position x");
DECLARE_SOA_COLUMN(Vy, vy, float);//, "collision vertex position y");
DECLARE_SOA_COLUMN(Vz, vz, float);//, "collision vertex position z");
DECLARE_SOA_COLUMN(Centrality, centrality, double);//, "Centrality");
DECLARE_SOA_COLUMN(TrackMultiplicity, trackmultiplicity, double);//, "Track Multiplicity");
DECLARE_SOA_COLUMN(RunNumber, runnumber, int);//, "Run Number");
DECLARE_SOA_COLUMN(Polarity, polarity, double);//, "Magnetic Field"

DECLARE_SOA_COLUMN(QxA_step1, qxa_step1, float);
DECLARE_SOA_COLUMN(QyA_step1, qya_step1, float);
DECLARE_SOA_COLUMN(QxC_step1, qxc_step1, float);
DECLARE_SOA_COLUMN(QyC_step1, qyc_step1, float);
DECLARE_SOA_COLUMN(Step1_done,step1_done, bool);
DECLARE_SOA_COLUMN(QxA_step2, qxa_step2, float);
DECLARE_SOA_COLUMN(QyA_step2, qya_step2, float);
DECLARE_SOA_COLUMN(QxC_step2, qxc_step2, float);
DECLARE_SOA_COLUMN(QyC_step2, qyc_step2, float);
DECLARE_SOA_COLUMN(Step2_done,step2_done, bool);
DECLARE_SOA_COLUMN(QxA_step3, qxa_step3, float);
DECLARE_SOA_COLUMN(QyA_step3, qya_step3, float);
DECLARE_SOA_COLUMN(QxC_step3, qxc_step3, float);
DECLARE_SOA_COLUMN(QyC_step3, qyc_step3, float);
DECLARE_SOA_COLUMN(Step3_done,step3_done, bool);
DECLARE_SOA_COLUMN(QxA_step4, qxa_step4, float);
DECLARE_SOA_COLUMN(QyA_step4, qya_step4, float);
DECLARE_SOA_COLUMN(QxC_step4, qxc_step4, float);
DECLARE_SOA_COLUMN(QyC_step4, qyc_step4, float);
DECLARE_SOA_COLUMN(Step4_done,step4_done, bool);
DECLARE_SOA_COLUMN(QxA_step5, qxa_step5, float);
DECLARE_SOA_COLUMN(QyA_step5, qya_step5, float);
DECLARE_SOA_COLUMN(QxC_step5, qxc_step5, float);
DECLARE_SOA_COLUMN(QyC_step5, qyc_step5, float);
DECLARE_SOA_COLUMN(Step5_done,step5_done, bool);
DECLARE_SOA_COLUMN(Vx_mean, vx_mean, float);
DECLARE_SOA_COLUMN(Vy_mean, vy_mean, float);
DECLARE_SOA_COLUMN(Vz_mean, vz_mean, float);
DECLARE_SOA_COLUMN(TrackMultiplicity_mean, trackmultiplicity_mean, double);


//add bools for each step, set to false but set to true after each calib step to keep track. (dont knos if this is the best way)
}//namespace myTable

DECLARE_SOA_TABLE(MyTableQvec, "AOD", "MYTABLEQVEC",
                  myTable::QxA,
                  myTable::QyA,
                  myTable::QxC,
                  myTable::QyC,
                  myTable::Vx,
                  myTable::Vy,
                  myTable::Vz,
                  myTable::Centrality,
                  myTable::TrackMultiplicity,
                  myTable::RunNumber,
                  myTable::Polarity,
                  myTable::QxA_step1,
                  myTable::QyA_step1,
                  myTable::QxC_step1,
                  myTable::QyC_step1,
                  myTable::Step1_done,
                  myTable::QxA_step2,
                  myTable::QyA_step2,
                  myTable::QxC_step2,
                  myTable::QyC_step2,
                  myTable::Step2_done,
                  myTable::QxA_step3,
                  myTable::QyA_step3,
                  myTable::QxC_step3,
                  myTable::QyC_step3,
                  myTable::Step3_done,
                  myTable::QxA_step4,
                  myTable::QyA_step4,
                  myTable::QxC_step4,
                  myTable::QyC_step4,
                  myTable::Step4_done,
                  myTable::QxA_step5,
                  myTable::QyA_step5,
                  myTable::QxC_step5,
                  myTable::QyC_step5,
                  myTable::Step5_done,
                  myTable::Vx_mean,
                  myTable::Vy_mean,
                  myTable::Vz_mean,
                  myTable::TrackMultiplicity_mean);
} // end of o2::aod

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct ProduceMyTable{
  Produces<o2::aod::MyTableQvec> mytableqvec;
  
  ConfigurableAxis axisCent{"axisCent", {90,0,90},"Centrality axis in 1% bins"};
  ConfigurableAxis axisCent10{"axisCent10", {9,0,90},"Centrality axis in 10% bins"};
  ConfigurableAxis axisQ{"axisQ", {100,-2,2},"Q vector (xy) in ZDC"};
  ConfigurableAxis axisVx{"axisVx", {5e3,0,0.05},"for Pos X of collision"};
  ConfigurableAxis axisVy{"axisVy", {6e3,0.35,0.41},"for Pos Y of collision"};
  ConfigurableAxis axisVx_mean{"axisVx_mean", {4000,-0.02,0.02},"for Pos X of collision"};
  ConfigurableAxis axisVy_mean{"axisVy_mean", {4000,-0.02,0.02},"for Pos Y of collision"};
  ConfigurableAxis axisVz{"axisVz", {20000,-12,12},"for vz of collision"}; // take 12 because we shift vi - <vi>
  ConfigurableAxis axisRun{"axisRun", {1e6,0,1e6},"for runNumber in ThnSparse"};
  
  O2_DEFINE_CONFIGURABLE(cfgCutVertex,        float,      10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMin,      float,      0.2f, "Minimal.q pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMax,      float,      10.0f, "Maximal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin,         float,      0.2f, "Minimal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax,         float,      3.0f, "Maximal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta,           float,      0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls,  float,      2.5, "Chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgUseNch,           bool,       false, "Use Nch for flow observables")
  O2_DEFINE_CONFIGURABLE(cfgNbootstrap,       int,        10, "Number of subsamples")
  O2_DEFINE_CONFIGURABLE(cfgMagField,         float,      99999, "Configurable magnetic field; default CCDB will be queried");
  O2_DEFINE_CONFIGURABLE(cfgCalibenergy,      std::string, "../runCalib/AnalysisResults.root", "path+file from energy calibration")
  O2_DEFINE_CONFIGURABLE(cfgQvecrecent,       std::string, "AnalysisResults_outCalib.root", "path+file from qvec recetring")
  O2_DEFINE_CONFIGURABLE(cfgSaveAs,           std::string, "outstep3", "path+file from qvec recetring")
  
  float qxa;
  float qya;
  float qxc;
  float qyc;
  
  float vx;
  float vy;
  float vz;
  
  float vx_mean;
  float vy_mean;
  float vz_mean;
  
  double centrality;
  double trackmultiplicity;
  double trackmultiplicity_mean;
  int runnumber;
  double polarity;
  
  float qxa_step1;
  float qya_step1;
  float qxc_step1;
  float qyc_step1;
  
  float qxa_step2;
  float qya_step2;
  float qxc_step2;
  float qyc_step2;
  
  float qxa_step3;
  float qya_step3;
  float qxc_step3;
  float qyc_step3;
  
  float qxa_step4;
  float qya_step4;
  float qxc_step4;
  float qyc_step4;
  
  float qxa_step5;
  float qya_step5;
  float qxc_step5;
  float qyc_step5;
  
  bool step1_done = false;
  bool step2_done = false;
  bool step3_done = false;
  bool step4_done = false;
  bool step5_done = false;
  
  int bin=0;
  int bin1=0;
  int lastRunNumberE = -1;
  int lastRunNumber_step1 = -1;
  int numFoundint = -1;
  int counter =0;
  
  // histograms to fill with fit values step 3
  // cross terms! Check if we need vz here..
    
//  TH2D* hFit_hQXA_Cent_vy;
//  TH2D* hFit_hQXA_Cent_vz;
//  
//  TH2D* hFit_hQYA_Cent_vx;
//  TH2D* hFit_hQYA_Cent_vz;
//
//  TH2D* hFit_hQXC_Cent_vy;
//  TH2D* hFit_hQXC_Cent_vz;
//  
//  TH2D* hFit_hQYC_Cent_vx;
//  TH2D* hFit_hQYC_Cent_vz;
  
  
  // open files with data from energy gain equalisation and recentering
     //cfgFitStep3->c_str(), "READ");
  
  //  Filters
  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtMin) && (aod::track::pt < cfgCutPtMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls);
  
  
  //define my.....
  using myCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>>;
  using myTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>>;
  
  inline int FindBinWithContentAlongAxis(THnSparse *hist, Int_t axis, Double_t targetContent) {
    // to find bin for given axis in THnSparse object
    TH1D* projection = hist->Projection(axis);
    int bin;
    Int_t nbins = projection->GetNbinsX();
    double binwidth = projection->GetXaxis()->GetBinWidth(1);
    
    // Iterate over bins along the specified axis
    for (Int_t i = 1; i <= nbins; ++i) {
      Double_t content = projection->GetBinCenter(i);
      
      if (std::abs(content - targetContent) < binwidth) {
        bin = i;
        return bin;
      }
      
    }
    
    // If the target content is not found along the specified axis
    std::cout << "Target content " << targetContent << " not found along axis " << axis << " in the histogram." << std::endl;
    delete projection;
    
    return 0;
  }
  
  inline float getQstep2(float q_step1, THnSparseF* mean, int cent, float vxmean, float vymean, float vzmean)
  {
    int total_counts = 0;
    float q;
    float total_q=0.;
    float mean_q =0.;
    int entries=0;
    
    // do we need to find the bins of v and i
    int bin_cent =   FindBinWithContentAlongAxis(mean, 0, cent);
    int bin_vxmean = FindBinWithContentAlongAxis(mean, 1, vxmean);
    int bin_vymean = FindBinWithContentAlongAxis(mean, 2, vymean);
    int bin_vzmean = FindBinWithContentAlongAxis(mean, 3, vzmean);
    
    mean->GetAxis(0)->SetRange(bin_cent, bin_cent);
    mean->GetAxis(1)->SetRange(bin_vxmean, bin_vxmean);
    mean->GetAxis(2)->SetRange(bin_vymean, bin_vymean);
    mean->GetAxis(3)->SetRange(bin_vzmean, bin_vzmean);
    
    TH1D* proQ = mean->Projection(4);
    
    for(int i=1; i<=proQ->GetNbinsX(); i++){
      entries = proQ->GetBinContent(i);
      q = proQ->GetBinCenter(i);
      total_q += q * entries;
      total_counts += entries;
    }
    
    // for total_counts =1 the mean is the val and you have 0.
    if(total_counts>1) mean_q = total_q/total_counts;
    //
    float q_step2 = q_step1 - mean_q;
    
    delete proQ;
    
    return q_step2;
  }
 inline float fitvalQstep3(THnSparseF* mean, int bin_run, int bin_cent) {
    
    //Make sure that the graph with (cent,run) (x,y) has fitval on z-axis. This has to be made once!! Then use fitvals in this function. But not fit for all events!
    if (!mean) {
      std::cout<< "Null pointer for THnSparseF object"<<std::endl;
      return 0.0;
    }
    
    int entries=0;
    float param=0.;
//    double chi2=0.;
    
    // look only in given centrality and run bin
    mean->GetAxis(0)->SetRange(bin_run, bin_run);
    mean->GetAxis(1)->SetRange(bin_cent, bin_cent);
    
    
    // project vi and Q in a 2d plot (yDim,xDim) : (Q,v)
    
    TH2D* proQv = mean->Projection(3,2);
    if (!proQv) {
      std::cout<<"Projection failed"<<std::endl;
      return 0.0;
    }
    //    LOGF(info, "nEntries total = %e", proQv->GetEntries());
    if(proQv->GetEntries()<1) {
      std::cout<< "No entries is Q(v) TH2D -> NO FIT POSSIBLE"<<std::endl;
      delete proQv;
      return 0.0;
    }
    
  //  TCanvas* c1 = new TCanvas("c1", " ", 800,800);
  //    c1->cd();
  //    proQv->Draw("COLZ");
  //    c1->SaveAs(TString::Format("outstep3/%s_%i_%i.png",mean->GetTitle(),bin_run-1,bin_cent));

    
    //check if x and y are what you expect
    //    LOGF(info, "bins X of ProQv = {%e,%e}", proQv->GetXaxis()->GetXmin(), proQv->GetXaxis()->GetXmax()); // [10:40:42][INFO] bins X of ProQv = {-2.000000e-02,2.000000e-02}
    //    LOGF(info, "bins Y of ProQv = {%e,%e}", proQv->GetYaxis()->GetXmin(), proQv->GetYaxis()->GetXmax()); // [10:40:42][INFO] bins Y of ProQv = {-2.000000e+00,2.000000e+00}
    
    TProfile* Q_v = new TProfile("Q_v", "", proQv->GetNbinsX(), proQv->GetXaxis()->GetXmin(), proQv->GetXaxis()->GetXmax());
    
    for(int i=1; i<=proQv->GetNbinsX(); i++){ // loop over all x-bins : v
      int count_Q_per_v=0;
      for(int j=1; j<=proQv->GetNbinsY(); j++){ // loop over all y-bins : Q
        entries = proQv->GetBinContent(i,j);
        if (entries>0){
          // hier komen we dus niet
          //          LOGF(info,"We fill v=%e with Q=%e %i times", proQv->GetXaxis()->GetBinCenter(i), proQv->GetYaxis()->GetBinCenter(j),entries);
          int counter = 0;
          count_Q_per_v++;
          while(counter < entries){
            //            if (counter==0) LOGF(info,"We fill v=%e with Q=%e %i times", proQv->GetXaxis()->GetBinCenter(i), proQv->GetYaxis()->GetBinCenter(j),entries);
            Q_v->Fill(proQv->GetXaxis()->GetBinCenter(i), proQv->GetYaxis()->GetBinCenter(j));
            counter++;
          }
        }
      }
      //ignore bins with only 1 entry because fit does not ignore.
      if(count_Q_per_v<2) Q_v->SetBinEntries(i,0);
    }
    
  //  TCanvas* c2 = new TCanvas("c2", " ", 800,800);

    if(Q_v->GetEntries()>1){
      TF1* fitfunc = new TF1("fitfunc", "pol1");
      fitfunc->SetParameter(0,0.0);
  //    TCanvas* canvas = new TCanvas("canvas", "Fit Result", 800, 600);
  //    canvas->cd();
      auto fitResult = Q_v->Fit(fitfunc, "QSFB0");
  //    delete canvas;

  //      c2->cd();
  //      Q_v->Draw();
  //      c2->SaveAs(TString::Format("outstep3/%s_fit_%i_%i.png", mean->GetTitle(),bin_run-1,bin_cent));

      if (fitResult->IsValid()) {
        param = fitfunc->GetParameter(1);
//        chi2 = fitfunc->GetChisquare();
        } else {
          std::cout<<"Error: Fit did not converge"<<std::endl;
          param = -0.0;
        }
      delete fitfunc;
    }
    
  //  TString title = mean->GetTitle();
  //  if (chi2>10) std::cout<< "CHi2 > 10!! For" << title.Data()<<" in bin ("<<bin_run-1<<","<< bin_cent*10<<" || Proceed with caution.. "<<std::endl;
    
    delete Q_v;
  //  delete c2;
  //  delete c1;
    delete proQv;
    
    return param;
  }

 inline void fillFitHist_step3(THnSparseF* mean, TH2D* global_hist)
  {
    if (!mean) {
      std::cout<<"Error: Null pointer for THnSparseF object"<<std::endl;
    }
    
    // runnumber bins
    int bins_x = mean->GetAxis(0)->GetNbins();
    // centrality bins
    int bins_y = mean->GetAxis(1)->GetNbins();
    
    TH1D* proRun = mean->Projection(0);
    
    int bin_fill=1;
    // fill histogtam with fitvalues
    for(int i=1; i<=bins_x; i++){ // x = runnumber
      if(proRun->GetBinContent(i)>=1){ // only fill for non-empty bins
        for(int j=1; j<=bins_y; j++){
          global_hist->SetBinContent(bin_fill,j,fitvalQstep3(mean,i,j));
          global_hist->GetXaxis()->SetBinLabel(bin_fill, TString::Format("%i",i-1));
  //          LOGF(info, " fit value: %e  | in bin (%i,%i)",fitvalQstep3(mean,i,j),i,j);
        }
        bin_fill++;
      }
    }
    
    delete proRun;
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
  
  inline float getQstep4(float q_step3, THnSparseF* mean, double cent, float vxmean, float vymean, float vzmean, double pol)
  {
    
    int total_counts = 0;
    float q;
    float total_q=0.;
    float mean_q =0.;
    int entries=0;
    
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
    
    TH1D* proQ = mean->Projection(5);
    //    LOGF(info, "nEntries proQ step2: %i", proQ->GetEntries());
    
    for(int i=1; i<=proQ->GetNbinsX(); i++){
      entries = proQ->GetBinContent(i);
      q = proQ->GetBinCenter(i);
      total_q += q * entries;
      total_counts += entries;
    }
    
    // for total_counts =1 the mean is the val and you have 0.
    if(total_counts>1) mean_q = total_q/total_counts;
    //
    float q_step4 = q_step3 - mean_q;
    
    delete proQ;
    
    return q_step4;
  }

  inline float getQstep5(float q_step4, THnSparseF* mean, int bin_runnmbr, double cent, double mult)
  {
    
    int total_counts = 0;
    float q;
    float total_q=0.;
    float mean_q =0.;
    int entries=0;
    
//    int bin_runmbr = FindBinWithContentAlongAxis(mean, 0, runN);
    int bin_cent =   FindBinWithContentAlongAxis(mean, 1, cent);
    int bin_mult =   FindBinWithContentAlongAxis(mean, 2, mult);
   
    mean->GetAxis(0)->SetRange(bin_runnmbr, bin_runnmbr);
    mean->GetAxis(1)->SetRange(bin_cent, bin_cent);
    mean->GetAxis(2)->SetRange(bin_mult, bin_mult);

    
    TH1D* proQ = mean->Projection(3);
    if (!proQ) {
      LOGF(error,"Projection failed");
      return 0.0;
    }
//    LOGF(info, "nEntries proQ step5: %e", proQ->GetEntries());
    
    for(int i=1; i<=proQ->GetNbinsX(); i++){
      entries = proQ->GetBinContent(i);
      q = proQ->GetBinCenter(i);
      total_q += q * entries;
      total_counts += entries;
    }
    
    // for total_counts =1 the mean is the val and you have 0.
    if(total_counts>1) mean_q = total_q/total_counts;
    //
    float q_step5 = q_step4 - mean_q;
    
    delete proQ;
    
    return q_step5;
  }
  
  inline void PrintMemoryConsumption() {
    ProcInfo_t procInfo;
    gSystem->GetProcInfo(&procInfo);
    
    std::cout << "Memory Usage:" << std::endl;
    std::cout << "  Resident memory: " << (float)procInfo.fMemResident/(1e6) << " MB" << std::endl;
    std::cout << "  Virtual memory: " << (float)procInfo.fMemVirtual/(1e6) << " MB" << std::endl;
}


  
  void process(myCollisions::iterator const& collision, aod::Zdcs const& zdcs, aod::BCs const& bcs, myTracks const& tracks)
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
    if(cent<0 || cent>90)
      return;
    
    
    // TODO: check how to pick up polarity for run 2!!!
//    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
//    auto field = (cfgMagField == 999999) ? getMagneticField(bc.timestamp()) : cfgMagField;
    
    //variable = (condition) ? expressionTrue : expressionFalse;
//    double polarity = (field<0) ? -0.5 : 0.5;
    
    // for now just take positive polarity
    double polarity = .5;
    
    if (collision.foundZDCId() >= 0) {
      // keep track of memory consumption
      if (counter % 200 == 0) PrintMemoryConsumption();
      
      int runNumber = collision.bc().runNumber();
      //      LOGF(info, "RUNNUMBER = %i", runNumber);
      //      float EcomZNA = collision.foundZDC().energyCommonZNA();
      //      float EcomZNC = collision.foundZDC().energyCommonZNC();
      float E1ZNA = collision.foundZDC().energySectorZNA()[0];
      float E2ZNA = collision.foundZDC().energySectorZNA()[1];
      float E3ZNA = collision.foundZDC().energySectorZNA()[2];
      float E4ZNA = collision.foundZDC().energySectorZNA()[3];
      float E1ZNC = collision.foundZDC().energySectorZNC()[0];
      float E2ZNC = collision.foundZDC().energySectorZNC()[1];
      float E3ZNC = collision.foundZDC().energySectorZNC()[2];
      float E4ZNC = collision.foundZDC().energySectorZNC()[3];
      
      
      // First step: start energy gain equalisation
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      TFile *f = TFile::Open(cfgCalibenergy->c_str(), "READ");
      if(!f || f->IsZombie()|| !f->IsOpen()) {
        LOGF(fatal, "Mean file does not exist | run in processCalib mode!!!!");
        return;
      }
      
      TProfile2D* hZNA_centCom = (TProfile2D*)f->Get("z-d-c-analysis/hZNA_mean_common_cent");
      
      if(!hZNA_centCom){
        LOGF(fatal, "Histo not found! Abort mission");
        return; }
      //      else{LOGF(info, "Histo found!");}
      
      if (runNumber != lastRunNumberE) {
        for(Int_t i=0; i<hZNA_centCom->GetXaxis()->GetNbins(); i++){
          const char* label = hZNA_centCom->GetXaxis()->GetBinLabel(i+1);
          //          LOGF(info, "RUN NUMBER %i | looking in bin %s", runNumber, label);
          int labelInt = atoi(label);
          if(labelInt == runNumber){bin=i+1;}
          //          LOGF(info, "looking in bin %i", bin);
        }
        numFoundint = atoi(hZNA_centCom->GetXaxis()->GetBinLabel(bin));
        //        LOGF(info, "RUN NUMBER %i | looking in bin %i : (%i)", runNumber, bin, numFoundint);
      }
      if(numFoundint!=runNumber)LOGF(error, "No match found for RUN NUMBER %e", runNumber);
      
      float ZNA_centCom = hZNA_centCom->GetBinContent(int(bin),int(cent)+1);
      TProfile2D* hZNC_centCom= (TProfile2D*)f->Get(Form("z-d-c-analysis/hZNC_mean_common_cent"));
      float ZNC_centCom = hZNC_centCom->GetBinContent(int(bin),int(cent)+1);
      
      float e1A=0.;
      float e2A=0.;
      float e3A=0.;
      float e4A=0.;
      
      float e1C=0.;
      float e2C=0.;
      float e3C=0.;
      float e4C=0.;
      
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
      
      
      std::vector<float> calE = {e1A, e2A, e3A , e4A, e1C, e2C, e3C , e4C};
      
      for(int i=0; i<4; i++){
        //With calibrated energies!!
        sumZNA += TMath::Power(calE[i],alpha);
        sumZNC += TMath::Power(calE[i+4],alpha);
        
        xEnZNA += ZDC_px[i] * TMath::Power(calE[i],alpha);
        yEnZNA += ZDC_py[i] * TMath::Power(calE[i],alpha);
        
        xEnZNC += ZDC_px[i] * TMath::Power(calE[i+4],alpha);
        yEnZNC += ZDC_py[i] * TMath::Power(calE[i+4],alpha);
        
      } // end of for(int i=0; i<4; i++)
      
      if(sumZNA>0){
        QxA = xEnZNA / sumZNA ;
        QyA = yEnZNA / sumZNA ;}
      
      if(sumZNC>0){
        QxC = xEnZNC / sumZNC ;
        QyC = yEnZNC / sumZNC ;}
      
      if (sumZNC<1e-11 || sumZNA<1e-11){
        LOGF(error, "No amplitudes found in ZNA or ZNC -> Q-vector is 0");
        return;
      }
      
      qxa=QxA;
      qya=QyA;
      qxc=QxC;
      qyc=QyC;
      
      vx = collision.posX();
      vy = collision.posY();
      vz = collision.posZ();
      centrality = cent;
      trackmultiplicity = tracks.size();
      runnumber = runNumber;
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      
      delete hZNC_centCom;
      delete hmeanE1A;
      delete hmeanE2A;
      delete hmeanE3A;
      delete hmeanE4A;
      delete hmeanE1C;
      delete hmeanE2C;
      delete hmeanE3C;
      delete hmeanE4C;
      
      // check if means for step 1 of re-centring exist
      
 

      TFile *f1 = TFile::Open(cfgQvecrecent->c_str(), "READ");
      if(!f1 || f1->IsZombie()|| !f1->IsOpen()) {
        if (counter<1) LOGF(error, "File for recentring does not exist | Table processed with non-recentered q-vectors!!!!");
        mytableqvec(qxa, qya, qxc, qyc, vx, vy, vz, centrality, trackmultiplicity, runnumber, polarity, 0., 0., 0.,0.,step1_done,0.,0.,0.,0.,step2_done,0.,0.,0.,0.,step3_done,0.,0.,0.,0.,step4_done, 0,0,0,0, step5_done, vx_mean, vy_mean, vz_mean, trackmultiplicity_mean);
        counter++;
        lastRunNumberE = runNumber;
        return;
      }
      
      
      
      TProfile2D* hQXA_mean_1perCent_Run = (TProfile2D*)f1->Get("consume-my-table/step1/hQXA_mean_1perCent_Run");
      if(!hQXA_mean_1perCent_Run  || hQXA_mean_1perCent_Run->GetEntries()<1 ){
        if (counter<1) {LOGF(error, "hQXA_mean_1perCent_Run not found or empty! STEP1 failed: Create table without re-centring");}
        mytableqvec(qxa, qya, qxc, qyc, vx, vy, vz, centrality, trackmultiplicity, runnumber, polarity, 0., 0., 0.,0.,step1_done,0.,0.,0.,0.,step2_done,0.,0.,0.,0.,step3_done,0.,0.,0.,0.,step4_done, 0,0,0,0, step5_done, vx_mean, vy_mean, vz_mean, trackmultiplicity_mean);
        counter++;
        lastRunNumberE = runNumber;
        return; }
      
      TProfile2D* hQYA_mean_1perCent_Run = (TProfile2D*)f1->Get("consume-my-table/step1/hQYA_mean_1perCent_Run");
      TProfile2D* hQXC_mean_1perCent_Run = (TProfile2D*)f1->Get("consume-my-table/step1/hQXC_mean_1perCent_Run");
      TProfile2D* hQYC_mean_1perCent_Run = (TProfile2D*)f1->Get("consume-my-table/step1/hQYC_mean_1perCent_Run");
      
      // retrieve runnumber and look up in the mean TProfile2D
      if (runNumber != lastRunNumber_step1) {
        lastRunNumber_step1 = runNumber;
        
        for(Int_t i=0; i<hQXA_mean_1perCent_Run->GetXaxis()->GetNbins(); i++){
          const char* label = hQXA_mean_1perCent_Run->GetXaxis()->GetBinLabel(i+1);
          //          LOGF(info, "RUN NUMBER %i | looking in bin %s", runNumber, label);
          int labelInt = atoi(label);
          if(labelInt == runNumber){bin=i+1;}
          //          LOGF(info, "looking in bin %i", bin);
        }
        
        // get runnumber from bin found in previous for loop
        numFoundint = atoi(hQXA_mean_1perCent_Run->GetXaxis()->GetBinLabel(bin));
      }
      
      // get mean value for given centrality
      float QXA_mean = hQXA_mean_1perCent_Run->GetBinContent(int(bin),int(cent)+1);
      float QYA_mean = hQYA_mean_1perCent_Run->GetBinContent(int(bin),int(cent)+1);
      float QXC_mean = hQXC_mean_1perCent_Run->GetBinContent(int(bin),int(cent)+1);
      float QYC_mean = hQYC_mean_1perCent_Run->GetBinContent(int(bin),int(cent)+1);
      
      
      qxa_step1 = qxa - QXA_mean;
      qya_step1 = qya - QYA_mean;
      qxc_step1 = qxc - QXC_mean;
      qyc_step1 = qyc - QYC_mean;
      
      
      step1_done = true;
      
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      
      THnSparseF* hQXA_mean_10perCent_v = (THnSparseF*)f1->Get("consume-my-table/step2/hQXA_mean_10perCent_v");
      if(!hQXA_mean_10perCent_v || hQXA_mean_10perCent_v->GetEntries()<1){
        if (counter<1) LOGF(error, "hQXA_mean_10perCent_v not found or empty! STEP 2 failed: Create table with only first step of re-centring");
        mytableqvec(qxa, qya, qxc, qyc, vx, vy, vz, centrality, trackmultiplicity, runnumber, polarity, qxa_step1, qya_step1, qxc_step1,qyc_step1,step1_done,0,0,0,0,step2_done,0,0,0,0,step3_done,0,0,0,0,step4_done, 0,0,0,0, step5_done, vx_mean, vy_mean, vz_mean, trackmultiplicity_mean);
        lastRunNumberE = runNumber;
        counter++;
        return; }
      
      THnSparseF* hQYA_mean_10perCent_v = (THnSparseF*)f1->Get("consume-my-table/step2/hQYA_mean_10perCent_v");
      THnSparseF* hQXC_mean_10perCent_v = (THnSparseF*)f1->Get("consume-my-table/step2/hQXC_mean_10perCent_v");
      THnSparseF* hQYC_mean_10perCent_v = (THnSparseF*)f1->Get("consume-my-table/step2/hQYC_mean_10perCent_v");
      
      TProfile* hvx_mean_Run = (TProfile*)f1->Get("consume-my-table/step2/hvx_mean_Run");
      if(!hvx_mean_Run || hvx_mean_Run->GetEntries()<1){
        if (counter<1) LOGF(error, "hQXC_mean_1perCent_Run not found or empty! Create table without re-centring");
        mytableqvec(qxa, qya, qxc, qyc, vx, vy, vz, centrality, trackmultiplicity, runnumber, polarity, 0., 0., 0.,0.,step1_done,0.,0.,0.,0.,step2_done,0.,0.,0.,0.,step3_done,0.,0.,0.,0.,step4_done, 0,0,0,0, step5_done, vx_mean, vy_mean, vz_mean, trackmultiplicity_mean);
        counter++;
        lastRunNumberE = runNumber;
        return; }
      
      TProfile* hvy_mean_Run = (TProfile*)f1->Get("consume-my-table/step2/hvy_mean_Run");
      TProfile* hvz_mean_Run = (TProfile*)f1->Get("consume-my-table/step2/hvz_mean_Run");
      
      //center vx, vy, vz around 0 on a run-by-run basis
      vx_mean = vx - hvx_mean_Run->GetBinContent(int(bin));
      vy_mean = vy - hvy_mean_Run->GetBinContent(int(bin));
      vz_mean = vz - hvz_mean_Run->GetBinContent(int(bin));
      
      
      qxa_step2 = getQstep2(qxa_step1, hQXA_mean_10perCent_v, centrality, vx_mean, vy_mean, vz_mean);
      qya_step2 = getQstep2(qya_step1, hQYA_mean_10perCent_v, centrality, vx_mean, vy_mean, vz_mean);
      qxc_step2 = getQstep2(qxc_step1, hQXC_mean_10perCent_v, centrality, vx_mean, vy_mean, vz_mean);
      qyc_step2 = getQstep2(qyc_step1, hQYC_mean_10perCent_v, centrality, vx_mean, vy_mean, vz_mean);
      
      step2_done=true;
      
      THnSparseF* hQXA_Run_10perCent_vy = (THnSparseF*)f1->Get("consume-my-table/step3/hQXA_Run_10perCent_vy");
      if(!hQXA_Run_10perCent_vy || hQXA_Run_10perCent_vy->GetEntries()<1){
        if (counter<1) LOGF(error, "hQXA_Run_10perCent_vy not found or empty! STEP 3 failed: Create table with only first and second step of re-centring");
        mytableqvec(qxa, qya, qxc, qyc, vx, vy, vz, centrality, trackmultiplicity, runnumber, polarity, qxa_step1, qya_step1, qxc_step1,qyc_step1,step1_done,qxa_step2,qya_step2,qxc_step2,qyc_step2,step2_done,0,0,0,0,step3_done,0,0,0,0,step4_done, 0,0,0,0, step5_done, vx_mean, vy_mean, vz_mean, trackmultiplicity_mean);
        counter++;
        lastRunNumberE = runNumber;
        return; }
      
      //       step 3: centrality 10%, run number and vi : fit <Q(vi)> = a * vi (a=fit parameter)
      //       So we need to do the fit in 10 centrality bins per run. Save mean values in 4dim plot: dim1=run, dim2=cent(10%), dim3=vi, dim4=Q
      //       Project dim4[Q] as function of dim3[vi] and fit with 0th order pol!

      TFile *f2 = TFile::Open("outFitStep3.root", "READ");
      if(!f2 || f2->IsZombie()|| !f2->IsOpen()) {
        if (counter<1) LOGF(error, "outFitStep3.root not found! STEP 3 failed: Create table with only first and second step of re-centring");
        mytableqvec(qxa, qya, qxc, qyc, vx, vy, vz, centrality, trackmultiplicity, runnumber, polarity, qxa_step1, qya_step1, qxc_step1,qyc_step1,step1_done,qxa_step2,qya_step2,qxc_step2,qyc_step2,step2_done,0,0,0,0,step3_done,0,0,0,0,step4_done, 0,0,0,0, step5_done, vx_mean, vy_mean, vz_mean, trackmultiplicity_mean);
        counter++;
        lastRunNumberE = runNumber;
        return;
      }
      
      TH2D* hFit_QXA_Run_Cent_vy = (TH2D*)f2->Get("hFit_QXA_Run_Cent_vy");
      if(!hFit_QXA_Run_Cent_vy || hFit_QXA_Run_Cent_vy->GetEntries()<1){
        if (counter<1) LOGF(error, "hFit_QXA_Run_Cent_vy not found or empty! STEP 3 failed: Create table with only first and second step of re-centring");
        mytableqvec(qxa, qya, qxc, qyc, vx, vy, vz, centrality, trackmultiplicity, runnumber, polarity, qxa_step1, qya_step1, qxc_step1,qyc_step1,step1_done,qxa_step2,qya_step2,qxc_step2,qyc_step2,step2_done,0,0,0,0,step3_done,0,0,0,0,step4_done, 0,0,0,0, step5_done, vx_mean, vy_mean, vz_mean, trackmultiplicity_mean);
        counter++;
        lastRunNumberE = runNumber;
        return; }
      
      TH2D* hFit_QYA_Run_Cent_vx = (TH2D*)f2->Get("hFit_QYA_Run_Cent_vx");
      TH2D* hFit_QXC_Run_Cent_vy = (TH2D*)f2->Get("hFit_QXC_Run_Cent_vy");
      TH2D* hFit_QYC_Run_Cent_vx = (TH2D*)f2->Get("hFit_QYC_Run_Cent_vx");
      
      int bin_run_step3 = -1;
      const char* label;
        for(Int_t i=0; i<hFit_QXA_Run_Cent_vy->GetXaxis()->GetNbins(); i++){ // sum over runs
          label = hFit_QXA_Run_Cent_vy->GetXaxis()->GetBinLabel(i+1);
          int labelInt = atoi(label);
          if(labelInt == runNumber){bin_run_step3=i+1; }//LOGF(info, "bin for run is found and set");}
        }
//      if(bin_run_step3 == -1) { LOGF(warning, "Run %i not found in fit step 3. Continue without..", atoi(label)); return;}
      
      int bin_cent = FindBinWithContentAlongAxis(hQXA_Run_10perCent_vy, 1, centrality);
      
      //       LOGF(info, "runnumber is: %i in bin %i | centrality is: %e in bin %i", runnumber, bin_run, centrality, bin_cent); // <- checked, goes well
      
      float q_mean_xa_step3 = hFit_QXA_Run_Cent_vy->GetBinContent(bin_run_step3, bin_cent)*vy_mean;
      float q_mean_xc_step3 = hFit_QYA_Run_Cent_vx->GetBinContent(bin_run_step3, bin_cent)*vx_mean;
      float q_mean_ya_step3 = hFit_QXC_Run_Cent_vy->GetBinContent(bin_run_step3, bin_cent)*vy_mean;
      float q_mean_yc_step3 = hFit_QYC_Run_Cent_vx->GetBinContent(bin_run_step3, bin_cent)*vx_mean;
      
//      LOGF(info, "q_mean_xa_step3 = %f", q_mean_xa_step3);
//      LOGF(info, "q_mean_xc_step3 = %f", q_mean_xc_step3);
//      LOGF(info, "q_mean_ya_step3 = %f", q_mean_ya_step3);
//      LOGF(info, "q_mean_yc_step3 = %f", q_mean_yc_step3);

      qxa_step3 = qxa_step2 - q_mean_xa_step3;
      qya_step3 = qya_step2 - q_mean_ya_step3;
      qxc_step3 = qxc_step2 - q_mean_xc_step3;
      qyc_step3 = qyc_step2 - q_mean_yc_step3;
      
      step3_done=true;
      
      delete hQXA_mean_1perCent_Run;
      delete hQYA_mean_1perCent_Run;
      delete hQXC_mean_1perCent_Run;
      delete hQYC_mean_1perCent_Run;

      
      TProfile2D* hmeanN_1perCent_Run = (TProfile2D*)f1->Get("consume-my-table/step4/hmeanN_1perCent_Run");
      THnSparseF* hQXA_mean_Magnet_10perCent_v = (THnSparseF*)f1->Get("consume-my-table/step4/hQXA_mean_Magnet_10perCent_v");
      if(!hmeanN_1perCent_Run || hmeanN_1perCent_Run->GetEntries()<1 || !hQXA_mean_Magnet_10perCent_v || hQXA_mean_Magnet_10perCent_v->GetEntries()<1){
        if (counter<1) LOGF(error, "Histos necessary for step 4 not found or empty! STEP 4 failed: Create table with only first, second and third step of re-centring");
        mytableqvec(qxa, qya, qxc, qyc, vx, vy, vz, centrality, trackmultiplicity, runnumber, polarity, qxa_step1, qya_step1, qxc_step1,qyc_step1,step1_done,qxa_step2,qya_step2,qxc_step2,qyc_step2,step2_done,qxa_step3,qya_step3,qxc_step3,qyc_step3,step3_done,0,0,0,0,step4_done, 0,0,0,0, step5_done, vx_mean, vy_mean, vz_mean, trackmultiplicity_mean);
        counter++;
        lastRunNumberE = runNumber;
        return; }
      
      THnSparseF* hQYA_mean_Magnet_10perCent_v = (THnSparseF*)f1->Get("consume-my-table/step4/hQYA_mean_Magnet_10perCent_v");
      THnSparseF* hQXC_mean_Magnet_10perCent_v = (THnSparseF*)f1->Get("consume-my-table/step4/hQXC_mean_Magnet_10perCent_v");
      THnSparseF* hQYC_mean_Magnet_10perCent_v = (THnSparseF*)f1->Get("consume-my-table/step4/hQYC_mean_Magnet_10perCent_v");
      
      qxa_step4 = getQstep4(qxa_step3, hQXA_mean_Magnet_10perCent_v, polarity, centrality, vx_mean, vy_mean, vz_mean);
      qya_step4 = getQstep4(qya_step3, hQYA_mean_Magnet_10perCent_v, polarity, centrality, vx_mean, vy_mean, vz_mean);
      qxc_step4 = getQstep4(qxc_step3, hQXC_mean_Magnet_10perCent_v, polarity, centrality, vx_mean, vy_mean, vz_mean);
      qyc_step4 = getQstep4(qyc_step3, hQYC_mean_Magnet_10perCent_v, polarity, centrality, vx_mean, vy_mean, vz_mean);
      
//      LOGF(info,"trak multiplicity before = %.2f", trackmultiplicity);
      trackmultiplicity_mean = trackmultiplicity - hmeanN_1perCent_Run->GetBinContent(bin,int(cent)+1);
//      LOGF(info,"mean trak multiplicity = %.2f", trackmultiplicity_mean);
      
      step4_done=true;
      
      delete hQXA_mean_Magnet_10perCent_v;
      delete hQYA_mean_Magnet_10perCent_v;
      delete hQXC_mean_Magnet_10perCent_v;
      delete hQYC_mean_Magnet_10perCent_v;
      delete hmeanN_1perCent_Run;
      
      THnSparseF* hQXA_mean_run_cent10_Mult = (THnSparseF*)f1->Get("consume-my-table/step5/hQXA_mean_run_cent10_Mult");
      if(!hQXA_mean_run_cent10_Mult || hQXA_mean_run_cent10_Mult->GetEntries()<1){
        if (counter<1) LOGF(error, "Histos necessary for step 5 not found or empty! STEP 5 failed: Create table with only first, second and third step of re-centring");
        mytableqvec(qxa, qya, qxc, qyc, vx, vy, vz, centrality, trackmultiplicity, runnumber, polarity, qxa_step1, qya_step1, qxc_step1,qyc_step1,step1_done,qxa_step2,qya_step2,qxc_step2,qyc_step2,step2_done,qxa_step3,qya_step3,qxc_step3,qyc_step3,step3_done,qxa_step4,qya_step4,qxc_step4,qyc_step4,step4_done, 0,0,0,0, step5_done, vx_mean, vy_mean, vz_mean, trackmultiplicity_mean);
        counter++;
        lastRunNumberE = runNumber;
        return; }
      
      THnSparseF* hQYA_mean_run_cent10_Mult = (THnSparseF*)f1->Get("consume-my-table/step5/hQYA_mean_run_cent10_Mult");
      THnSparseF* hQXC_mean_run_cent10_Mult = (THnSparseF*)f1->Get("consume-my-table/step5/hQXC_mean_run_cent10_Mult");
      THnSparseF* hQYC_mean_run_cent10_Mult = (THnSparseF*)f1->Get("consume-my-table/step5/hQYC_mean_run_cent10_Mult");
      
      int bin_run5 = FindBinWithContentAlongAxis(hQXA_mean_run_cent10_Mult, 0, runnumber+0.5);
      
      qxa_step5 = getQstep5(qxa_step4, hQXA_mean_run_cent10_Mult, bin_run5, centrality, trackmultiplicity_mean);
      qya_step5 = getQstep5(qya_step4, hQYA_mean_run_cent10_Mult, bin_run5, centrality, trackmultiplicity_mean);
      qxc_step5 = getQstep5(qxc_step4, hQXC_mean_run_cent10_Mult, bin_run5, centrality, trackmultiplicity_mean);
      qyc_step5 = getQstep5(qyc_step4, hQYC_mean_run_cent10_Mult, bin_run5, centrality, trackmultiplicity_mean);
      
      step5_done = true;
      
      delete hQXA_mean_run_cent10_Mult;
      delete hQYA_mean_run_cent10_Mult;
      delete hQXC_mean_run_cent10_Mult;
      delete hQYC_mean_run_cent10_Mult;
    
      
      
      THnSparseF* hQXA_step6 = (THnSparseF*)f1->Get("consume-my-table/step6/hQXA_step6");
      if(!hQXA_step6 || hQXA_step6->GetEntries()<1){
        if (counter<1) LOGF(error, "Histos necessary for step 6 not found or empty! STEP 6 failed: Create table with only first 5 steps of re-centring (this is enough)");
        mytableqvec(qxa, qya, qxc, qyc, vx, vy, vz, centrality, trackmultiplicity, runnumber, polarity, qxa_step1, qya_step1, qxc_step1,qyc_step1,step1_done,qxa_step2,qya_step2,qxc_step2,qyc_step2,step2_done,qxa_step3,qya_step3,qxc_step3,qyc_step3,step3_done,qxa_step4,qya_step4,qxc_step4,qyc_step4,step4_done, qxa_step5,qya_step5,qxc_step5,qyc_step5, step5_done, vx_mean, vy_mean, vz_mean, trackmultiplicity_mean);
        counter++;
        lastRunNumberE = runNumber;
        return;
      }
      
      f->Close();
      f1->Close();
      f2->Close();
      
    }// end collision found ZDC
    
  
    
  }
};

struct ConsumeMyTable{
  
  ConfigurableAxis axisCent{"axisCent", {90,0,90},"Centrality axis in 1% bins"};
  ConfigurableAxis axisCent10{"axisCent10", {9,0,90},"Centrality axis in 10% bins"};
  ConfigurableAxis axisQ{"axisQ", {100,-2,2},"Q vector (xy) in ZDC"};
  ConfigurableAxis axisVx{"axisVx", {5e3,0,0.05},"for Pos X of collision"};
  ConfigurableAxis axisVy{"axisVy", {6e3,0.35,0.41},"for Pos Y of collision"};
  ConfigurableAxis axisVx_mean{"axisVx_mean", {400,-0.01,0.01},"for Pos X of collision"};
  ConfigurableAxis axisVy_mean{"axisVy_mean", {400,-0.01,0.01},"for Pos Y of collision"};
  ConfigurableAxis axisVz{"axisVz", {200,-12,12},"for vz of collision"}; // take 12 because we shift vi - <vi>
  ConfigurableAxis axisRun{"axisRun", {1e6,0,1e6},"for runNumber in ThnSparse"};
  ConfigurableAxis axisPolarity{"axisPolarity", {2,-1,1},"Magnet Polarity"};
  ConfigurableAxis axisMult{"axisMult", {5000,0,5000},"Track Multiplicity"};
  ConfigurableAxis axisMult_mean{"axisMult_mean", {5000,-2500,2500},"Track Multiplicity mean per run"};
  
  
  //  //define output
  HistogramRegistry registry{"Registry"};
  
  void init(InitContext const&)
  {
    // NOTE (todo) : add qx vs qy for each step
    // NOTE (todo) : add flagg to table to see if calib steps are included
    // TODO : add plots to show Q(vars) used in calib before and after to see effect of recentring.
    
    // before recentring
    registry.add("step0/hZNA_Qx_vs_Qy","hZNA_Qx_vs_Qy", {HistType::kTH2F, {axisQ,axisQ}});
    registry.add("step0/hZNC_Qx_vs_Qy","hZNC_Qx_vs_Qy", {HistType::kTH2F, {axisQ,axisQ}});
    // TODO: Hier ook de Q(cent) Q(vi) histos toevoegen!
    
    //step1
    registry.add("step1/hZNA_Qx_vs_Qy","hZNA_Qx_vs_Qy", {HistType::kTH2F, {axisQ,axisQ}});
    registry.add("step1/hZNC_Qx_vs_Qy","hZNC_Qx_vs_Qy", {HistType::kTH2F, {axisQ,axisQ}});
    
    registry.add("step1/hQXA_mean_1perCent_Run","hQXA_mean_1perCent_Run", {HistType::kTProfile2D, {{1,0.,1.} ,axisCent}});
    registry.add("step1/hQYA_mean_1perCent_Run","hQYA_mean_1perCent_Run", {HistType::kTProfile2D, {{1,0.,1.} ,axisCent}});
    registry.add("step1/hQXC_mean_1perCent_Run","hQXC_mean_1perCent_Run", {HistType::kTProfile2D, {{1,0.,1.} ,axisCent}});
    registry.add("step1/hQYC_mean_1perCent_Run","hQYC_mean_1perCent_Run", {HistType::kTProfile2D, {{1,0.,1.} ,axisCent}});
    
    registry.add("step1/hQXA_vs_cent_before","hQXA_vs_cent_before", {HistType::kTProfile, {axisCent}});
    registry.add("step1/hQYA_vs_cent_before","hQYA_vs_cent_before", {HistType::kTProfile, {axisCent}});
    registry.add("step1/hQXC_vs_cent_before","hQXC_vs_cent_before", {HistType::kTProfile, {axisCent}});
    registry.add("step1/hQYC_vs_cent_before","hQYC_vs_cent_before", {HistType::kTProfile, {axisCent}});
    
    registry.add("step1/hQXA_vs_cent_after","hQXA_vs_cent_after", {HistType::kTProfile, {axisCent}});
    registry.add("step1/hQYA_vs_cent_after","hQYA_vs_cent_after", {HistType::kTProfile, {axisCent}});
    registry.add("step1/hQXC_vs_cent_after","hQXC_vs_cent_after", {HistType::kTProfile, {axisCent}});
    registry.add("step1/hQYC_vs_cent_after","hQYC_vs_cent_after", {HistType::kTProfile, {axisCent}});
    
    //step2
    
    registry.add("step2/hZNA_Qx_vs_Qy","hZNA_Qx_vs_Qy", {HistType::kTH2F, {axisQ,axisQ}});
    registry.add("step2/hZNC_Qx_vs_Qy","hZNC_Qx_vs_Qy", {HistType::kTH2F, {axisQ,axisQ}});
    
    registry.add("step2/hQXA_mean_10perCent_v","hQXA_mean_10perCent_v", {HistType::kTHnSparseF, {axisCent10, axisVx_mean, axisVy_mean, axisVz, axisQ}});
    registry.add("step2/hQYA_mean_10perCent_v","hQYA_mean_10perCent_v", {HistType::kTHnSparseF, {axisCent10, axisVx_mean, axisVy_mean, axisVz, axisQ}});
    registry.add("step2/hQXC_mean_10perCent_v","hQXC_mean_10perCent_v", {HistType::kTHnSparseF, {axisCent10, axisVx_mean, axisVy_mean, axisVz, axisQ}});
    registry.add("step2/hQYC_mean_10perCent_v","hQYC_mean_10perCent_v", {HistType::kTHnSparseF, {axisCent10, axisVx_mean, axisVy_mean, axisVz, axisQ}});
    
    registry.add("step2/hvx_mean_Run","hvx_mean_Run", {HistType::kTProfile, {{1,0.,1.}}});
    registry.add("step2/hvy_mean_Run","hvy_mean_Run", {HistType::kTProfile, {{1,0.,1.}}});
    registry.add("step2/hvz_mean_Run","hvz_mean_Run", {HistType::kTProfile, {{1,0.,1.}}});
    
    registry.add("step2/hQXA_vs_vx_before","hQXA_vs_vx_before", {HistType::kTProfile, {axisVx_mean}});
    registry.add("step2/hQYA_vs_vx_before","hQYA_vs_vx_before", {HistType::kTProfile, {axisVx_mean}});
    registry.add("step2/hQXC_vs_vx_before","hQXC_vs_vx_before", {HistType::kTProfile, {axisVx_mean}});
    registry.add("step2/hQYC_vs_vx_before","hQYC_vs_vx_before", {HistType::kTProfile, {axisVx_mean}});
    
    registry.add("step2/hQXA_vs_vy_before","hQXA_vs_vy_before", {HistType::kTProfile, {axisVy_mean}});
    registry.add("step2/hQYA_vs_vy_before","hQYA_vs_vy_before", {HistType::kTProfile, {axisVy_mean}});
    registry.add("step2/hQXC_vs_vy_before","hQXC_vs_vy_before", {HistType::kTProfile, {axisVy_mean}});
    registry.add("step2/hQYC_vs_vy_before","hQYC_vs_vy_before", {HistType::kTProfile, {axisVy_mean}});
    
    registry.add("step2/hQXA_vs_vz_before","hQXA_vs_vz_before", {HistType::kTProfile, {axisVz}});
    registry.add("step2/hQYA_vs_vz_before","hQYA_vs_vz_before", {HistType::kTProfile, {axisVz}});
    registry.add("step2/hQXC_vs_vz_before","hQXC_vs_vz_before", {HistType::kTProfile, {axisVz}});
    registry.add("step2/hQYC_vs_vz_before","hQYC_vs_vz_before", {HistType::kTProfile, {axisVz}});
    
    registry.add("step2/hQXA_vs_vx_after","hQXA_vs_vx_after", {HistType::kTProfile, {axisVx_mean}});
    registry.add("step2/hQYA_vs_vx_after","hQYA_vs_vx_after", {HistType::kTProfile, {axisVx_mean}});
    registry.add("step2/hQXC_vs_vx_after","hQXC_vs_vx_after", {HistType::kTProfile, {axisVx_mean}});
    registry.add("step2/hQYC_vs_vx_after","hQYC_vs_vx_after", {HistType::kTProfile, {axisVx_mean}});
    
    registry.add("step2/hQXA_vs_vy_after","hQXA_vs_vy_after", {HistType::kTProfile, {axisVy_mean}});
    registry.add("step2/hQYA_vs_vy_after","hQYA_vs_vy_after", {HistType::kTProfile, {axisVy_mean}});
    registry.add("step2/hQXC_vs_vy_after","hQXC_vs_vy_after", {HistType::kTProfile, {axisVy_mean}});
    registry.add("step2/hQYC_vs_vy_after","hQYC_vs_vy_after", {HistType::kTProfile, {axisVy_mean}});
    
    registry.add("step2/hQXA_vs_vz_after","hQXA_vs_vz_after", {HistType::kTProfile, {axisVz}});
    registry.add("step2/hQYA_vs_vz_after","hQYA_vs_vz_after", {HistType::kTProfile, {axisVz}});
    registry.add("step2/hQXC_vs_vz_after","hQXC_vs_vz_after", {HistType::kTProfile, {axisVz}});
    registry.add("step2/hQYC_vs_vz_after","hQYC_vs_vz_after", {HistType::kTProfile, {axisVz}});
    
    //step 3
    registry.add("step3/hZNA_Qx_vs_Qy","hZNA_Qx_vs_Qy", {HistType::kTH2F, {axisQ,axisQ}});
    registry.add("step3/hZNC_Qx_vs_Qy","hZNC_Qx_vs_Qy", {HistType::kTH2F, {axisQ,axisQ}});
    
    registry.add("step3/hQXA_mean_10perCent_Run","hQXA_mean_10perCent_Run", {HistType::kTProfile2D, {{1,0.,1.} ,axisCent10}});
    registry.add("step3/hQYA_mean_10perCent_Run","hQYA_mean_10perCent_Run", {HistType::kTProfile2D, {{1,0.,1.} ,axisCent10}});
    registry.add("step3/hQXC_mean_10perCent_Run","hQXC_mean_10perCent_Run", {HistType::kTProfile2D, {{1,0.,1.} ,axisCent10}});
    registry.add("step3/hQYC_mean_10perCent_Run","hQYC_mean_10perCent_Run", {HistType::kTProfile2D, {{1,0.,1.} ,axisCent10}});
    
    // cross terms
    registry.add("step3/hQXA_Run_10perCent_vy","hQXA_Run_10perCent_vy", {HistType::kTHnSparseF, {axisRun,axisCent10,axisVy_mean, axisQ}});
//    registry.add("step3/hQXA_Run_10perCent_vz","hQXA_Run_10perCent_vz", {HistType::kTHnSparseF, {axisRun,axisCent10,axisVz, axisQ}});
    
    registry.add("step3/hQYA_Run_10perCent_vx","hQYA_Run_10perCent_vx", {HistType::kTHnSparseF, {axisRun,axisCent10,axisVx_mean, axisQ}});
//    registry.add("step3/hQYA_Run_10perCent_vz","hQYA_Run_10perCent_vz", {HistType::kTHnSparseF, {axisRun,axisCent10,axisVz, axisQ}});
    
    registry.add("step3/hQXC_Run_10perCent_vy","hQXC_Run_10perCent_vx", {HistType::kTHnSparseF, {axisRun,axisCent10,axisVy_mean, axisQ}});
//    registry.add("step3/hQXC_Run_10perCent_vz","hQXC_Run_10perCent_vz", {HistType::kTHnSparseF, {axisRun,axisCent10,axisVz, axisQ}});
    
    registry.add("step3/hQYC_Run_10perCent_vx","hQYC_Run_10perCent_vx", {HistType::kTHnSparseF, {axisRun,axisCent10,axisVx_mean, axisQ}});
//    registry.add("step3/hQYC_Run_10perCent_vz","hQYC_Run_10perCent_vz", {HistType::kTHnSparseF, {axisRun,axisCent10,axisVz, axisQ}});
    
    //step 4
    registry.add("step4/hZNA_Qx_vs_Qy","hZNA_Qx_vs_Qy", {HistType::kTH2F, {axisQ,axisQ}});
    registry.add("step4/hZNC_Qx_vs_Qy","hZNC_Qx_vs_Qy", {HistType::kTH2F, {axisQ,axisQ}});
    
    registry.add("step4/hmeanN_1perCent_Run","hmeanN_1perCent_Run", {HistType::kTProfile2D, {{1,0.,1.} ,axisCent}});
    
    registry.add("step4/hQXA_mean_Magnet_10perCent_v","hQXA_mean_Magnet_10perCent_v", {HistType::kTHnSparseF, {axisPolarity, axisCent10, axisVx_mean, axisVy_mean, axisVz, axisQ}});
    registry.add("step4/hQYA_mean_Magnet_10perCent_v","hQYA_mean_Magnet_10perCent_v", {HistType::kTHnSparseF, {axisPolarity, axisCent10, axisVx_mean, axisVy_mean, axisVz, axisQ}});
    registry.add("step4/hQXC_mean_Magnet_10perCent_v","hQXC_mean_Magnet_10perCent_v", {HistType::kTHnSparseF, {axisPolarity, axisCent10, axisVx_mean, axisVy_mean, axisVz, axisQ}});
    registry.add("step4/hQYC_mean_Magnet_10perCent_v","hQYC_mean_Magnet_10perCent_v", {HistType::kTHnSparseF, {axisPolarity, axisCent10, axisVx_mean, axisVy_mean, axisVz, axisQ}});
    
    
    
    //step 5
    registry.add("step5/hZNA_Qx_vs_Qy","hZNA_Qx_vs_Qy", {HistType::kTH2F, {axisQ,axisQ}});
    registry.add("step5/hZNC_Qx_vs_Qy","hZNC_Qx_vs_Qy", {HistType::kTH2F, {axisQ,axisQ}});
    
    registry.add("step5/hQXA_mean_run_cent10_Mult","hQXA_mean_run_cent10_Mult", {HistType::kTHnSparseF, {axisRun, axisCent10, axisMult_mean, axisQ}});
    registry.add("step5/hQYA_mean_run_cent10_Mult","hQYA_mean_run_cent10_Mult", {HistType::kTHnSparseF, {axisRun, axisCent10, axisMult_mean, axisQ}});
    registry.add("step5/hQXC_mean_run_cent10_Mult","hQXC_mean_run_cent10_Mult", {HistType::kTHnSparseF, {axisRun, axisCent10, axisMult_mean, axisQ}});
    registry.add("step5/hQYC_mean_run_cent10_Mult","hQYC_mean_run_cent10_Mult", {HistType::kTHnSparseF, {axisRun, axisCent10, axisMult_mean, axisQ}});
    
    registry.add("step5/hQXA_vs_cent","hQXA_vs_cent", {HistType::kTProfile, {axisCent}});
    registry.add("step5/hQYA_vs_cent","hQYA_vs_cent", {HistType::kTProfile, {axisCent}});
    registry.add("step5/hQXC_vs_cent","hQXC_vs_cent", {HistType::kTProfile, {axisCent}});
    registry.add("step5/hQYC_vs_cent","hQYC_vs_cent", {HistType::kTProfile, {axisCent}});
    
//    registry.add("step5/hQXA_vs_vx","hQXA_vs_vx", {HistType::kTProfile, {axisVx_mean}});
    registry.add("step5/hQYA_vs_vx","hQYA_vs_vx", {HistType::kTProfile, {axisVx_mean}});
//    registry.add("step5/hQXC_vs_vx","hQXC_vs_vx", {HistType::kTProfile, {axisVx_mean}});
    registry.add("step5/hQYC_vs_vx","hQYC_vs_vx", {HistType::kTProfile, {axisVx_mean}});
    
    registry.add("step5/hQXA_vs_vy","hQXA_vs_vy", {HistType::kTProfile, {axisVy_mean}});
//    registry.add("step5/hQYA_vs_vy","hQYA_vs_vy", {HistType::kTProfile, {axisVy_mean}});
    registry.add("step5/hQXC_vs_vy","hQXC_vs_vy", {HistType::kTProfile, {axisVy_mean}});
//    registry.add("step5/hQYC_vs_vy","hQYC_vs_vy", {HistType::kTProfile, {axisVy_mean}});
    
    registry.add("step5/hQXA_vs_vz","hQXA_vs_vz", {HistType::kTProfile, {axisVz}});
    registry.add("step5/hQYA_vs_vz","hQYA_vs_vz", {HistType::kTProfile, {axisVz}});
    registry.add("step5/hQXC_vs_vz","hQXC_vs_vz", {HistType::kTProfile, {axisVz}});
    registry.add("step5/hQYC_vs_vz","hQYC_vs_vz", {HistType::kTProfile, {axisVz}});
    
                 
    // recentered q-vectors
    
    registry.add("hZNA_Qx_vs_Qy","hZNA_Qx_vs_Qy", {HistType::kTH2F, {axisQ,axisQ}});
    registry.add("hZNC_Qx_vs_Qy","hZNC_Qx_vs_Qy", {HistType::kTH2F, {axisQ,axisQ}});
    registry.add("hStep","hStep", {HistType::kTH1D, {{5,0.,5.}}});
    
  }
  
  void process(aod::MyTableQvec const& mytableqvecs)
  {
    
    for(auto& qvec : mytableqvecs){
      
      registry.fill(HIST("step0/hZNA_Qx_vs_Qy"),qvec.qxa(), qvec.qya());
      registry.fill(HIST("step0/hZNC_Qx_vs_Qy"),qvec.qxc(), qvec.qyc());
      
      //step1
      registry.get<TProfile2D>(HIST("step1/hQXA_mean_1perCent_Run"))->Fill(Form("%d", qvec.runnumber()),qvec.centrality(),qvec.qxa(),1);
      registry.get<TProfile2D>(HIST("step1/hQXC_mean_1perCent_Run"))->Fill(Form("%d", qvec.runnumber()),qvec.centrality(),qvec.qxc(),1);
      registry.get<TProfile2D>(HIST("step1/hQYA_mean_1perCent_Run"))->Fill(Form("%d", qvec.runnumber()),qvec.centrality(),qvec.qya(),1);
      registry.get<TProfile2D>(HIST("step1/hQYC_mean_1perCent_Run"))->Fill(Form("%d", qvec.runnumber()),qvec.centrality(),qvec.qyc(),1);
      
      registry.get<TProfile>(HIST("step1/hQXA_vs_cent_before"))->Fill(qvec.centrality(), qvec.qxa());
      registry.get<TProfile>(HIST("step1/hQYA_vs_cent_before"))->Fill(qvec.centrality(), qvec.qya());
      registry.get<TProfile>(HIST("step1/hQXC_vs_cent_before"))->Fill(qvec.centrality(), qvec.qxc());
      registry.get<TProfile>(HIST("step1/hQYC_vs_cent_before"))->Fill(qvec.centrality(), qvec.qyc());
      
      //      std::cout<<" - vx: "<<qvec.vx()<<" - vy: "<<qvec.vy()<<std::endl;
      registry.get<TProfile>(HIST("step2/hvx_mean_Run"))->Fill(Form("%d",qvec.runnumber()), qvec.vx());
      registry.get<TProfile>(HIST("step2/hvy_mean_Run"))->Fill(Form("%d",qvec.runnumber()), qvec.vy());
      registry.get<TProfile>(HIST("step2/hvz_mean_Run"))->Fill(Form("%d",qvec.runnumber()), qvec.vz());
      
      if(qvec.step1_done()){
        registry.fill(HIST("step2/hQXA_mean_10perCent_v"), qvec.centrality(), qvec.vx_mean(), qvec.vy_mean(), qvec.vz_mean(), qvec.qxa_step1());
        registry.fill(HIST("step2/hQYA_mean_10perCent_v"), qvec.centrality(), qvec.vx_mean(), qvec.vy_mean(), qvec.vz_mean(), qvec.qya_step1());
        registry.fill(HIST("step2/hQXC_mean_10perCent_v"), qvec.centrality(), qvec.vx_mean(), qvec.vy_mean(), qvec.vz_mean(), qvec.qxc_step1());
        registry.fill(HIST("step2/hQYC_mean_10perCent_v"), qvec.centrality(), qvec.vx_mean(), qvec.vy_mean(), qvec.vz_mean(), qvec.qyc_step1());
        
        registry.fill(HIST("step1/hZNA_Qx_vs_Qy"),qvec.qxa_step1(), qvec.qya_step1());
        registry.fill(HIST("step1/hZNC_Qx_vs_Qy"),qvec.qxc_step1(), qvec.qyc_step1());
        
        registry.fill(HIST("step1/hQXA_vs_cent_after"), qvec.centrality(), qvec.qxa_step1());
        registry.fill(HIST("step1/hQYA_vs_cent_after"), qvec.centrality(), qvec.qya_step1());
        registry.fill(HIST("step1/hQXC_vs_cent_after"), qvec.centrality(), qvec.qxc_step1());
        registry.fill(HIST("step1/hQYC_vs_cent_after"), qvec.centrality(), qvec.qyc_step1());
        
        registry.get<TProfile>(HIST("step2/hQXA_vs_vx_before"))->Fill(qvec.vx_mean(), qvec.qxa_step1());
        registry.get<TProfile>(HIST("step2/hQYA_vs_vx_before"))->Fill(qvec.vx_mean(), qvec.qya_step1());
        registry.get<TProfile>(HIST("step2/hQXC_vs_vx_before"))->Fill(qvec.vx_mean(), qvec.qxc_step1());
        registry.get<TProfile>(HIST("step2/hQYC_vs_vx_before"))->Fill(qvec.vx_mean(), qvec.qyc_step1());
        
        registry.get<TProfile>(HIST("step2/hQXA_vs_vy_before"))->Fill(qvec.vy_mean(), qvec.qxa_step1());
        registry.get<TProfile>(HIST("step2/hQYA_vs_vy_before"))->Fill(qvec.vy_mean(), qvec.qya_step1());
        registry.get<TProfile>(HIST("step2/hQXC_vs_vy_before"))->Fill(qvec.vy_mean(), qvec.qxc_step1());
        registry.get<TProfile>(HIST("step2/hQYC_vs_vy_before"))->Fill(qvec.vy_mean(), qvec.qyc_step1());
        
        registry.get<TProfile>(HIST("step2/hQXA_vs_vz_before"))->Fill(qvec.vz_mean(), qvec.qxa_step1());
        registry.get<TProfile>(HIST("step2/hQYA_vs_vz_before"))->Fill(qvec.vz_mean(), qvec.qya_step1());
        registry.get<TProfile>(HIST("step2/hQXC_vs_vz_before"))->Fill(qvec.vz_mean(), qvec.qxc_step1());
        registry.get<TProfile>(HIST("step2/hQYC_vs_vz_before"))->Fill(qvec.vz_mean(), qvec.qyc_step1());
        
        registry.fill(HIST("hStep"), .5, 1);
        
      }
      if(qvec.step2_done()){
        
        //        std::cout<<" - qxc_step1: "<<qvec.qxc_step1()<<" - qxc_step2: "<<qvec.qxc_step2()<<" - qxc_step1 - qxc_step2: "<<qvec.qxc_step1() - qvec.qxc_step2()<<std::endl;
        
        registry.fill(HIST("step2/hZNA_Qx_vs_Qy"),qvec.qxa_step2(), qvec.qya_step2());
        registry.fill(HIST("step2/hZNC_Qx_vs_Qy"),qvec.qxc_step2(), qvec.qyc_step2());
        
        registry.get<TProfile>(HIST("step2/hQXA_vs_vx_after"))->Fill(qvec.vx_mean(), qvec.qxa_step2());
        registry.get<TProfile>(HIST("step2/hQYA_vs_vx_after"))->Fill(qvec.vx_mean(), qvec.qya_step2());
        registry.get<TProfile>(HIST("step2/hQXC_vs_vx_after"))->Fill(qvec.vx_mean(), qvec.qxc_step2());
        registry.get<TProfile>(HIST("step2/hQYC_vs_vx_after"))->Fill(qvec.vx_mean(), qvec.qyc_step2());
        
        registry.get<TProfile>(HIST("step2/hQXA_vs_vy_after"))->Fill(qvec.vy_mean(), qvec.qxa_step2());
        registry.get<TProfile>(HIST("step2/hQYA_vs_vy_after"))->Fill(qvec.vy_mean(), qvec.qya_step2());
        registry.get<TProfile>(HIST("step2/hQXC_vs_vy_after"))->Fill(qvec.vy_mean(), qvec.qxc_step2());
        registry.get<TProfile>(HIST("step2/hQYC_vs_vy_after"))->Fill(qvec.vy_mean(), qvec.qyc_step2());
        
        registry.get<TProfile>(HIST("step2/hQXA_vs_vz_after"))->Fill(qvec.vz_mean(), qvec.qxa_step2());
        registry.get<TProfile>(HIST("step2/hQYA_vs_vz_after"))->Fill(qvec.vz_mean(), qvec.qya_step2());
        registry.get<TProfile>(HIST("step2/hQXC_vs_vz_after"))->Fill(qvec.vz_mean(), qvec.qxc_step2());
        registry.get<TProfile>(HIST("step2/hQYC_vs_vz_after"))->Fill(qvec.vz_mean(), qvec.qyc_step2());
        
        //        std::cout<<" - vx: "<<qvec.vx()<<" - vy: "<<qvec.vy()<<std::endl;
        
        // cross terms -> Denk niet vz!
        registry.fill(HIST("step3/hQXA_Run_10perCent_vy"), qvec.runnumber(), qvec.centrality(), qvec.vy_mean(), qvec.qxa_step2());
//        registry.fill(HIST("step3/hQXA_Run_10perCent_vz"), qvec.runnumber(), qvec.centrality(), qvec.vz_mean(), qvec.qxa_step2());
        
        registry.fill(HIST("step3/hQYA_Run_10perCent_vx"), qvec.runnumber(), qvec.centrality(), qvec.vx_mean(), qvec.qya_step2());
//        registry.fill(HIST("step3/hQYA_Run_10perCent_vz"), qvec.runnumber(), qvec.centrality(), qvec.vz_mean(), qvec.qya_step2());
        
        registry.fill(HIST("step3/hQXC_Run_10perCent_vy"), qvec.runnumber(), qvec.centrality(), qvec.vy_mean(), qvec.qxc_step2());
//        registry.fill(HIST("step3/hQXC_Run_10perCent_vz"), qvec.runnumber(), qvec.centrality(), qvec.vz_mean(), qvec.qxc_step2());
        
        registry.fill(HIST("step3/hQYC_Run_10perCent_vx"), qvec.runnumber(), qvec.centrality(), qvec.vx_mean(), qvec.qyc_step2());
//        registry.fill(HIST("step3/hQYC_Run_10perCent_vz"), qvec.runnumber(), qvec.centrality(), qvec.vz_mean(), qvec.qyc_step2());
        
        // fill with zeros for fitting later.
        registry.get<TProfile2D>(HIST("step3/hQXA_mean_10perCent_Run"))->Fill(Form("%d", qvec.runnumber()),qvec.centrality(),0,1);
        registry.get<TProfile2D>(HIST("step3/hQXC_mean_10perCent_Run"))->Fill(Form("%d", qvec.runnumber()),qvec.centrality(),0,1);
        registry.get<TProfile2D>(HIST("step3/hQYA_mean_10perCent_Run"))->Fill(Form("%d", qvec.runnumber()),qvec.centrality(),0,1);
        registry.get<TProfile2D>(HIST("step3/hQYC_mean_10perCent_Run"))->Fill(Form("%d", qvec.runnumber()),qvec.centrality(),0,1);
        
        
        registry.fill(HIST("hStep"), 1.5, 1);
        
      }
      
      
      if(qvec.step3_done()){
        
        //        std::cout<<" - qxc_step1: "<<qvec.qxc_step1()<<" - qxc_step3: "<<qvec.qxc_step3()<<" - qxc_step1 - qxc_step3: "<<qvec.qxc_step1() - qvec.qxc_step3()<<std::endl;
        registry.fill(HIST("step3/hZNA_Qx_vs_Qy"),qvec.qxa_step3(), qvec.qya_step3());
        registry.fill(HIST("step3/hZNC_Qx_vs_Qy"),qvec.qxc_step3(), qvec.qyc_step3());
        
        registry.fill(HIST("step4/hQXA_mean_Magnet_10perCent_v"), qvec.polarity(), qvec.centrality(), qvec.vx_mean(), qvec.vy_mean(), qvec.vz_mean(), qvec.qxa_step3());
        registry.fill(HIST("step4/hQYA_mean_Magnet_10perCent_v"), qvec.polarity(), qvec.centrality(), qvec.vx_mean(), qvec.vy_mean(), qvec.vz_mean(), qvec.qya_step3());
        registry.fill(HIST("step4/hQXC_mean_Magnet_10perCent_v"), qvec.polarity(), qvec.centrality(), qvec.vx_mean(), qvec.vy_mean(), qvec.vz_mean(), qvec.qxc_step3());
        registry.fill(HIST("step4/hQYC_mean_Magnet_10perCent_v"), qvec.polarity(), qvec.centrality(), qvec.vx_mean(), qvec.vy_mean(), qvec.vz_mean(), qvec.qyc_step3());
        
        registry.get<TProfile2D>(HIST("step4/hmeanN_1perCent_Run"))->Fill(Form("%d",qvec.runnumber()), qvec.centrality(), qvec.trackmultiplicity(),1);
        
        registry.fill(HIST("hStep"), 2.5, 1);
      }
      
      if(qvec.step4_done()){
        
        registry.fill(HIST("step4/hZNA_Qx_vs_Qy"),qvec.qxa_step4(), qvec.qya_step4());
        registry.fill(HIST("step4/hZNC_Qx_vs_Qy"),qvec.qxc_step4(), qvec.qyc_step4());
        
        registry.fill(HIST("step5/hQXA_mean_run_cent10_Mult"), qvec.runnumber(), qvec.centrality(), qvec.trackmultiplicity_mean(), qvec.qxa_step4());
        registry.fill(HIST("step5/hQYA_mean_run_cent10_Mult"), qvec.runnumber(), qvec.centrality(), qvec.trackmultiplicity_mean(), qvec.qya_step4());
        registry.fill(HIST("step5/hQXC_mean_run_cent10_Mult"), qvec.runnumber(), qvec.centrality(), qvec.trackmultiplicity_mean(), qvec.qxc_step4());
        registry.fill(HIST("step5/hQYC_mean_run_cent10_Mult"), qvec.runnumber(), qvec.centrality(), qvec.trackmultiplicity_mean(), qvec.qyc_step4());
        
        registry.fill(HIST("hStep"), 3.5, 1);
      }
      
      //recentered q-vectors
      if(qvec.step1_done() && qvec.step2_done() && qvec.step3_done() && qvec.step4_done() && qvec.step5_done()){
        
        registry.get<TProfile>(HIST("step5/hQXA_vs_cent"))->Fill(qvec.centrality(), qvec.qxa_step5());
        registry.get<TProfile>(HIST("step5/hQYA_vs_cent"))->Fill(qvec.centrality(), qvec.qya_step5());
        registry.get<TProfile>(HIST("step5/hQXC_vs_cent"))->Fill(qvec.centrality(), qvec.qxc_step5());
        registry.get<TProfile>(HIST("step5/hQYC_vs_cent"))->Fill(qvec.centrality(), qvec.qyc_step5());
        
//        registry.get<TProfile>(HIST("step5/hQXA_vs_vx"))->Fill(qvec.vx_mean(), qvec.qxa_step5());
        registry.get<TProfile>(HIST("step5/hQYA_vs_vx"))->Fill(qvec.vx_mean(), qvec.qya_step5());
//        registry.get<TProfile>(HIST("step5/hQXC_vs_vx"))->Fill(qvec.vx_mean(), qvec.qxc_step5());
        registry.get<TProfile>(HIST("step5/hQYC_vs_vx"))->Fill(qvec.vx_mean(), qvec.qyc_step5());
        
        registry.get<TProfile>(HIST("step5/hQXA_vs_vy"))->Fill(qvec.vy_mean(), qvec.qxa_step5());
//        registry.get<TProfile>(HIST("step5/hQYA_vs_vy"))->Fill(qvec.vy_mean(), qvec.qya_step5());
        registry.get<TProfile>(HIST("step5/hQXC_vs_vy"))->Fill(qvec.vy_mean(), qvec.qxc_step5());
//        registry.get<TProfile>(HIST("step5/hQYC_vs_vy"))->Fill(qvec.vy_mean(), qvec.qyc_step5());
        
        registry.get<TProfile>(HIST("step5/hQXA_vs_vz"))->Fill(qvec.vz_mean(), qvec.qxa_step5());
        registry.get<TProfile>(HIST("step5/hQYA_vs_vz"))->Fill(qvec.vz_mean(), qvec.qya_step5());
        registry.get<TProfile>(HIST("step5/hQXC_vs_vz"))->Fill(qvec.vz_mean(), qvec.qxc_step5());
        registry.get<TProfile>(HIST("step5/hQYC_vs_vz"))->Fill(qvec.vz_mean(), qvec.qyc_step5());
        
        registry.fill(HIST("step5/hZNA_Qx_vs_Qy"),qvec.qxa_step5(), qvec.qya_step5());
        registry.fill(HIST("step5/hZNC_Qx_vs_Qy"),qvec.qxc_step5(), qvec.qyc_step5());
        
        registry.fill(HIST("hZNA_Qx_vs_Qy"),qvec.qxa_step5(), qvec.qya_step5());
        registry.fill(HIST("hZNC_Qx_vs_Qy"),qvec.qxc_step5(), qvec.qyc_step5());
        
        registry.fill(HIST("hStep"), 4.5, 1);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    //    adaptAnalysisTask<ZDCAnalysis>(cfgc)};
    adaptAnalysisTask<ProduceMyTable>(cfgc),
    adaptAnalysisTask<ConsumeMyTable>(cfgc),
  };
  
}


// backup


//  inline float getmeanQstep3(THnSparseF* mean, double runnumber, double cent, bool isY=false)
//  {
//    if(isY){return 0.;} // to include Y histos.. something goes wrong.
//
//    int total_counts = 0;
//    float q=0.;
//    float total_q=0.;
//    int entries=0;
//    float param=0.;
//
//    // do we need to find the bins of v and i
//    int bin_run    = FindBinWithContentAlongAxis(mean, 0, runnumber+.5);
//    int bin_cent   = FindBinWithContentAlongAxis(mean, 1, cent);
////    int bin_v      = FindBinWithContentAlongAxis(mean, axisV, v);
//
//    mean->GetAxis(0)->SetRange(bin_run, bin_run);
//    mean->GetAxis(1)->SetRange(bin_cent, bin_cent);
//
//    // project vi and Q in a 2d plot (yDim,xDim) : (Q,v)
//    TH2D* proQ = mean->Projection(3,2);
//
//    TH1D* hprof  = new TH1D("hprof","Profile of Q vs vi", proQ->GetNbinsX(),proQ->GetXaxis()->GetXmin(),proQ->GetXaxis()->GetXmax());
//
//    for(int i=1; i<=proQ->GetNbinsX(); i++){
//      total_q=0;
//      total_counts=0;
//      for(int j=1; j<=proQ->GetNbinsY(); j++){
//        entries = proQ->GetBinContent(i,j);
//        q = proQ->GetYaxis()->GetBinCenter(j);
//        total_q += q * entries;
//        total_counts += entries;
//      }
//      if(total_counts>0) {
//        hprof->SetBinContent(i, total_q/total_counts);
////        LOGF(info, "x: %e   | y: %e",proQ->GetXaxis()->GetBinCenter(i), total_q/total_counts);
//      }
//
//    }
//
//    TF1* fitfunc = new TF1("fitfunc", "pol0");
//    hprof->Fit(fitfunc,"Q");
//    param = fitfunc->GetParameter(0);
//
//    delete proQ; delete hprof;
//
//    return param;
//  }



// do we need to find the bins of v and i
//      int bin_run    = FindBinWithContentAlongAxis(mean, 0, runnumber+.5);
//      int bin_cent   = FindBinWithContentAlongAxis(mean, 1, cent);
//
//      if (bin_run < 0 || bin_run >= mean->GetAxis(0)->GetNbins() ||
//          bin_cent < 0 || bin_cent >= mean->GetAxis(1)->GetNbins()) {
//          // Error: Axis range check
//          std::cerr << "Error: Invalid bin indices\n";
//          return 0.0;
//      }
