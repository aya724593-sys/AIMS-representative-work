#include <algorithm>
#include <chrono>
#include <cmath>
#include <csignal>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <random>
#include <string>
#include <vector>
//---------------------------------------------------------------
//---------------------------------------------------------------
#include "Pythia8/HeavyIons.h"
#include "Pythia8/Pythia.h"
//---------------------------------------------------------------
//---------------------------------------------------------------
#include "TBrowser.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TClassTable.h"
#include "TComplex.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TUnixSystem.h"
#include "TVector2.h"
#include "TVector3.h"
//---------------------------------------------------------------
#ifdef DEBUG
#define LogDebug(message) std::cout << "[DEBUG] " << message << std::endl
#else
#define LogDebug(message)
#endif
//---------------------------------------------------------------
using namespace Pythia8;
//---------------------------------------------------------------
char hname[500];
char histname[500];
//---------------------------------------------------------------
//-------------Flow Vectors Constructions------------------------
TComplex Q(TComplex Harmonic[4][20][20], Int_t n, Int_t p) {
  Int_t m = fabs(n);
  TComplex Qnp = 0.0;
  if (n >= 0) {
    Qnp = Harmonic[0][p][m];
  } else {
    Qnp = TComplex::Conjugate(Harmonic[0][p][m]);
  }
  return Qnp;
}
//---------------------------------------------------------------
//---------Two Particle Correlations Constructions---------------
double Two(TComplex Harmonic[4][20][20], Int_t n1, Int_t n2) {
  double two =
      (Q(Harmonic, n1, 1) * Q(Harmonic, n2, 1) - Q(Harmonic, n1 + n2, 2)).Re();
  return two;
}
//---------------------------------------------------------------
//-------------------Histograms definitions----------------------
TH1D *histogram[10];
TProfile *TP_d00_c22t0;
//---------------------------------------------------------------
int main() {
  //-------------------------------------------------------------
  //-----------------Histograms definitions----------------------
  //-------------------------------------------------------------
  histogram[0] = new TH1D("Multa", "Multa", 5000, 0.0, 5000.0);
  histogram[0]->Sumw2();

  histogram[1] = new TH1D("Imbac", "Imbac", 1000, 0.0, 50.0);
  histogram[1]->Sumw2();

  histogram[2] = new TH1D("Npart", "Npart", 1000, 0.0, 1000.0);
  histogram[2]->Sumw2();
  //-------------------------------------------------------------
  TP_d00_c22t0 = new TProfile("TP_d00_c22t0", "TP_d00_c22t0", 1000, 0.0, 1000.);
  TP_d00_c22t0->Sumw2();
  //-------------------------------------------------------------
  //-------------------------------------------------------------
  Pythia pythia;
  //-------------------------------------------------------------
  //-------------------------------------------------------------
  pythia.readFile("input.cmnd");
  //-------------------------------------------------------------
  int nError = pythia.mode("Main:timesAllowErrors");
  const bool countErrors = (nError > 0 ? true : false);
  //-------------------------------------------------------------
  //-------------------------------------------------------------
  int nEvent = pythia.mode("Main:numberOfEvents");
  //-------------------------------------------------------------
  // Get the current time as a duration since the epoch
  auto currentTime = std::chrono::system_clock::now().time_since_epoch();
  auto currentTimeMs =
      std::chrono::duration_cast<std::chrono::milliseconds>(currentTime)
          .count();
  // Seed the random number generator with the current time
  std::mt19937_64 generator(currentTimeMs);
  // Generate a random odd number
  std::uniform_int_distribution<int> distribution(1, 50000);
  int randomNumber = distribution(generator);
  if (randomNumber % 2 == 0) {
    randomNumber++; // Make sure it's odd
  }
  std::cout << "Random odd number: " << randomNumber << std::endl;

  // Retrieve beam particle IDs
  int idBeamA = pythia.settings.mode("Beams:idA"); // Beam 1
  int idBeamB = pythia.settings.mode("Beams:idB"); // Beam 2

  pythia.readString("Random:setSeed = on");
  pythia.readString(("Random:seed = " + std::to_string(randomNumber)).c_str());
  //-------------------------------------------------------------
  pythia.init();
  //-------------------------------------------------------------
  // Pythia8ToHepMC toHepMC("my_fifo");
  //-------------------------------------------------------------
  for (int i = 0; i < nEvent; ++i) {
    if (!pythia.next())
      continue;
    //-----------------------------------------------------------
    double Evweight = pythia.info.weight();
    //-----------------------------------------------------------
    double Imb = 0.0;
    int Npart = 0.0;
    if (fabs(idBeamA) != 2212 && fabs(idBeamA) != 2212) {
      Imb = pythia.info.hiInfo->b();
      Npart = pythia.info.hiInfo->nPartTarg() + pythia.info.hiInfo->nPartProj();
    }
    //-----------------------------------------------------------
    //--------------------Charged Particles counting-------------
    int Ncharge = 0;
    for (int i = 0; i < pythia.event.size(); ++i) {
      Particle &p = pythia.event[i];
      //--------------------------------------------------------
      if (p.isFinal()) {
        //--------------------------------------------------------
        if (p.isCharged()) {
          if (p.pT() > 0.2 && p.pT() < 3.0 && fabs(p.eta()) < 0.5)
            Ncharge++;
        }
      }
    }
    //------------------------------------------------------------
    if (Ncharge < 1.0)
      continue;
    //------------------------------------------------------------
    //--------------------Histograms filling----------------------
    histogram[0]->Fill(Ncharge);
    histogram[1]->Fill(Imb);
    histogram[2]->Fill(Npart);
    //------------------------------------------------------------
    //---------------Flow Vectors x,y initialization--------------
    double QxQ[5][20][20] = {0};
    double QyQ[5][20][20] = {0};
    for (int k = 0; k < 5; k++) {
      for (int i = 0; i < 20; i++) {
        for (int j = 0; j < 20; j++) {
          QxQ[k][i][j] = 0.0;
          QyQ[k][i][j] = 0.0;
        }
      }
    }
    //------------------------------------------------------------
    //------------------------------------------------------------
    //------------------------------------------------------------
    for (int i = 0; i < pythia.event.size(); ++i) {
      Particle &p = pythia.event[i];
      //---------------------------------------------------------
      if (p.isFinal()) {
        // Keep track of the sum of waights
        //-------------------------------------------------------
        double pT = p.pT();
        double Eta = p.eta();
        double Rap = p.y();
        double Phi = p.phi();
        int PID = p.id();
        int ch = p.charge();
        //-------------------------------------------------------
        //--------------------Particle selection-----------------
        if (pT < 0.2 || pT > 4.0)
          continue;
        if (fabs(Eta) > 1.0)
          continue;
        //-------------------------------------------------------
        //--------------------Particle ID assignment-------------
        int id = -999;
        if (PID == 211)
          id = 0;
        else if (PID == -211)
          id = 1;
        else if (PID == 321)
          id = 2;
        else if (PID == -321)
          id = 3;
        else if (PID == 2212)
          id = 4;
        else if (PID == -2212)
          id = 5;
        else
          continue;
        if (id < 0)
          continue;
        //-------------------------------------------------------
        //--------------------Flow Vectors x,y filling-----------
        for (int i = 0; i < 20; i++) {
          for (int n = 0; n < 20; n++) {
            QxQ[0][i][n] += pow(1.0, double(i)) * TMath::Cos(double(n) * Phi);
            QyQ[0][i][n] += pow(1.0, double(i)) * TMath::Sin(double(n) * Phi);
          }
        }
        //-------------------------------------------------------
        //--------------------Eta bin assignment-----------------
        int ett2 = -999;
        if (Eta > -1.00 && Eta < -0.20) {
          ett2 = 1;
        } else if (Eta > 0.20 && Eta < 1.00) {
          ett2 = 2;
        } else {
          ett2 = -999;
        }
        if (ett2 < 0)
          continue;
        //-------------------------------------------------------
        //--------------------Flow Vectors x,y filling-----------
        for (int i = 0; i < 20; i++) {
          for (int n = 0; n < 20; n++) {
            QxQ[ett2][i][n] +=
                pow(1.0, double(i)) * TMath::Cos(double(n) * Phi);
            QyQ[ett2][i][n] +=
                pow(1.0, double(i)) * TMath::Sin(double(n) * Phi);
          }
        }
        //-------------------------------------------------------
      }
    }
    //----------------------------------------------------------
    //---------------Flow Vectors x,y conversion to complex-----
    TComplex QnR[4][20][20] = {0};
    for (int e = 0; e < 4; e++) {
      for (int p = 0; p < 20; p++) {
        for (int n = 0; n < 20; n++) {
          QnR[e][p][n] = TComplex(QxQ[e][p][n], QyQ[e][p][n]);
        }
      }
    }
    //----------------------------------------------------------
    //---------------Two Particle Correlations calculations-----
    double d00_w00t0 = Two(QnR, 0, 0);
    double d00_c22t0 = Two(QnR, 2, -2);
    double d00_c33t0 = Two(QnR, 3, -3);
    //----------------------------------------------------------
    //---------------Two Particle Correlations filling-----------
    if (d00_w00t0 > 0.00000001) {
      TP_d00_c22t0->Fill(Ncharge, (d00_c22t0) / (d00_w00t0), (d00_w00t0));
    }
  }
  //-------------------------------------------------------------
  //-------------------------------------------------------------
  //-------------------------------------------------------------
  pythia.stat();
  //-------------------------------------------------------------
  //-------------------------------------------------------------
  std::cout << "The pythia end of the code" << "\n";
  //-------------------------------------------------------------
  //-------------------------------------------------------------
  //-------------------------------------------------------------
  //---------------Output file creation--------------------------
  TFile *output = new TFile("Output_PYTHIA8_file.root", "RECREATE");
  //-------------------------------------------------------------
  //---------------Histograms writing----------------------------
  for (int i = 0; i < 3; i++) {
    histogram[i]->Write();
  }
  TP_d00_c22t0->Write();
  //-------------------------------------------------------------
  output->Close();
  //-------------------------------------------------------------
  //-------------------------------------------------------------
  return 0;
}
