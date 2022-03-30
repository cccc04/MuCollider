
#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesStream.h"
#include "classes/SortableObject.h"
#include "external/ExRootAnalysis/ExRootProgressBar.h"
#include "external/ExRootAnalysis/ExRootTask.h"
#include "external/ExRootAnalysis/ExRootTreeBranch.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootTreeWriter.h"
#include "modules/Delphes.h"
#endif
#include <TBranch.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLeaf.h>
#include <TStopwatch.h>
#include <TTree.h>
#include <cmath>
#include <iostream>
#include <ROOT/TTreeProcessorMT.hxx>
#include <ROOT/TThreadedObject.hxx>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

void nb(Int_t workers, string fileii, TString files, Int_t En)
{
  TStopwatch end;//Timer
  //Input format, e.g.:.x nb.C+(12,"rs/dm","dm",2) for 12 processors, rs/dm directory with dm file name, 2 TeV Ecm
  string filef = ".root";
  string Enn = to_string(En);
  TString Eenn = to_string(En);
  TString fil = files + Eenn;
  string filei = fileii + Enn;
  string file = filei + filef;
  ROOT::EnableImplicitMT(workers);//Mulithread

  // Create one TThreadedObject per histogram to fill during the processing of the tree

  // Create a TTreeProcessorMT: specify the file and the tree in it
  ROOT::TTreeProcessorMT tp(file, "Delphes");

  //Create Histograms

  ROOT::TThreadedObject<TH1F> hmet0("MET0", "MissingET", 150, 0, En * 1000 / 2 + En * 80);
  /**/ ROOT::TThreadedObject<TH1F> hmet("MET1", "MissingET", 150, 0, En * 1000 / 2 + En * 80);
  ROOT::TThreadedObject<TH1F> hmetg("MET", "MissingET", 150, 0, En * 1000 / 2 + En * 80);
  ROOT::TThreadedObject<TH1F> hzmu("Invariant mass un", "Reconstructed Z", 150, 0, 400);
  ROOT::TThreadedObject<TH1F> hzm("Invariant mass", "Reconstructed Z", 150, 0, 400);
  ROOT::TThreadedObject<TH1F> hjpt("PT", "JetPt", 150, 0, En * 1000 / 2 + En * 80);
  ROOT::TThreadedObject<TH1F> hee("E", "Extra Energy", 150, 0, En * 1000 / 2 + En * 80);
  ROOT::TThreadedObject<TH1F> hee2("E2", "Extra Energy2", 150, 0, En * 1000 / 2 + En * 80);
  ROOT::TThreadedObject<TH1F> hlpt1("PT1", "Lepton with higher PT", 150, 0, En * 1000);
  ROOT::TThreadedObject<TH1F> hlpt2("PT2", "Lepton with lower PT", 150, 0, En * 300);
  ROOT::TThreadedObject<TH1F> hm1("Mass", "MissingMass", 150, 0, 1000 * En);
  ROOT::TThreadedObject<TH1F> hm2("Mass", "GenMissingMass", 150, 0, 1000 * En);
  ROOT::TThreadedObject<TH1F> hnop("N", "Number of Particles", 8, 0, 8);
  ROOT::TThreadedObject<TH1F> hdfmz("#Delta#phi", "#Delta#phi", 24, -4, 4);
  ROOT::TThreadedObject<TH1F> hmetph("Phi", "MissingET", 150, -4, 4);
  ROOT::TThreadedObject<TH1F> hdptmz("PT", "MET-ZPT", 150, 0, En * 500);
  ROOT::TThreadedObject<TH1F> hht("HT", "HT", 150, 0, En * 500);
  ROOT::TThreadedObject<TH2F> hmm("Mass", "MissingMass", 150, 0, 1000 * En, 150, 0, En * 1000 / 2 + En * 80);
  ROOT::TThreadedObject<TH2F> hmm1("Mass", "MissingMass", 150, 0, 1000 * En, 150, 0, En * 1000 / 2 + En * 80);

  //Main function

  auto myFunction = [&](TTreeReader &myReader) {
    TTreeReaderArray<float> Met(myReader, "MissingET.MET");
    TTreeReaderArray<int> pid(myReader, "Track.PID");
    TTreeReaderValue<int> Ts(myReader, "Track_size");
    TTreeReaderValue<int> Tws(myReader, "Tower_size");
    TTreeReaderValue<int> Pts(myReader, "Particle_size");
    TTreeReaderValue<int> Jts(myReader, "KTjet_size");
    TTreeReaderValue<int> Fms(myReader, "ForwardMuon_size");
    TTreeReaderValue<int> Phs(myReader, "Photon_size");
    TTreeReaderValue<int> Efts(myReader, "EFlowTrack_size");
    TTreeReaderValue<int> Efphs(myReader, "EFlowPhoton_size");
    TTreeReaderValue<int> Efnhs(myReader, "EFlowNeutralHadron_size");
    TTreeReaderArray<float> Tp(myReader, "Track.PT");
    TTreeReaderArray<float> Te(myReader, "Track.Eta");
    TTreeReaderArray<float> Tph(myReader, "Track.Phi");
    TTreeReaderArray<float> Twe(myReader, "Track.P");
    TTreeReaderArray<float> Top(myReader, "Tower.ET");
    TTreeReaderArray<float> Toe(myReader, "Tower.Eta");
    TTreeReaderArray<float> Toph(myReader, "Tower.Phi");
    TTreeReaderArray<float> Towe(myReader, "Tower.E");
    TTreeReaderArray<float> Mte(myReader, "MissingET.Eta");
    TTreeReaderArray<float> Mtp(myReader, "MissingET.Phi");
    TTreeReaderArray<int> ppid(myReader, "Particle.PID");
    TTreeReaderArray<float> pe(myReader, "Particle.E");
    TTreeReaderArray<float> ppx(myReader, "Particle.Px");
    TTreeReaderArray<float> ppy(myReader, "Particle.Py");
    TTreeReaderArray<float> ppz(myReader, "Particle.Pz");
    TTreeReaderArray<int> sta(myReader, "Particle.Status");
    TTreeReaderArray<float> jtpt(myReader, "KTjet.PT");
    TTreeReaderArray<float> ht(myReader, "ScalarHT.HT");
    TTreeReaderArray<float> fmpt(myReader, "ForwardMuon.PT");
    TTreeReaderArray<float> fme(myReader, "ForwardMuon.Eta");
    TTreeReaderArray<float> fmph(myReader, "ForwardMuon.Phi");
    TTreeReaderArray<float> phe(myReader, "Photon.E");

    while(myReader.Next())
    {
        auto myhmet0 = hmet0.Get();
        auto myhmet = hmet.Get();
        auto myhm1 = hm1.Get();
        auto myhm2 = hm2.Get();
        auto myhht = hht.Get();
        auto myhmetg = hmetg.Get();
        auto myhzmu = hzmu.Get();
        auto myhzm = hzm.Get();
        auto myhjpt = hjpt.Get();
        auto myhee2 = hee2.Get();
        auto myhee = hee.Get();
        auto myhlpt2 = hlpt2.Get();
        auto myhlpt1 = hlpt1.Get();
        auto myhnop = hnop.Get();
        auto myhdfmz = hdfmz.Get();
        auto myhmetph = hmetph.Get();
        auto myhdptmz = hdptmz.Get();
        auto myhmm = hmm.Get();
        auto myhmm1 = hmm1.Get();
        // For performance reasons, a copy of the pointer associated to this thread on the stack is used
        myhmet0->Fill(Met[0]);//Overall missingET
        if((*Ts >= 2 && ((pid[0] == 13 && pid[1] == -13) || (pid[0] == -13 && pid[1] == 13) || (pid[0] == 11 && pid[1] == -11) || (pid[0] == -11 && pid[1] == 11))) || (*Ts < 2 && *Tws == 2))
        // Selection on final particles. Only includes muon pairs or electron pairs identified by "Tracks". Or the cases with 2 hits in "Towers" (sometimes "Tracks" fails to identify but "Towers" always does). 
        {
          Float_t A;
          Float_t B;
          Float_t C;
          Float_t E;
          Float_t F;
          Float_t G;
          Float_t H;
          Float_t J;
          Float_t K;
          Float_t L;
          Float_t M;
          Float_t O;
          Float_t P;
          Float_t Q;
          Float_t N;
          //Float_t R;
          //Float_t S;
          //Float_t T;
          Float_t U;
          Float_t V;
          Float_t W;
          Float_t X;
          Float_t Y;
          Float_t Z;
          Float_t px1;
          Float_t py1;
          Float_t pz1;
          Float_t px2;
          Float_t py2;
          Float_t pz2;
          Float_t Pe1;
          Float_t Px1;
          Float_t Py1;
          Float_t Pz1;
          Float_t Pe2;
          Float_t Px2;
          Float_t Py2;
          Float_t Pz2;
          Float_t mm;
          Float_t wr0;
          Float_t wr1;
            if(*Ts >= 2 && ((pid[0] == 13 && pid[1] == -13) || (pid[0] == -13 && pid[1] == 13) || (pid[0] == 11 && pid[1] == -11) || (pid[0] == -11 && pid[1] == 11)))// Tracks case
            {
              A = Tp[0];
              B = Tp[1];
              E = Te[0];
              F = Te[1];
              G = Tph[0];
              H = Tph[1];
              wr0 = Twe[0];
              wr1 = Twe[1];
            }
            if(*Ts < 2 && *Tws == 2)// Tower case
            {
              A = Top[0];
              B = Top[1];
              E = Toe[0];
              F = Toe[1];
              G = Toph[0];
              H = Toph[1];
              wr0 = Towe[0];
              wr1 = Towe[1];
            }
            if(A > B)// Selection of the leading particle
            {
              K = A; // leading particle pt
              L = B; // sub leading particle pt
              N = E; // K's eta
              O = F; // L's eta
              P = G; // K's phi
              Q = H; // L's phi
            }
            if(A <= B)
            {
              K = B;
              L = A;
              N = F;
              O = E;
              P = H;
              Q = G;
            }
            // momentum components calculation
            px1 = K * cos(P);
            px2 = L * cos(Q);
            py1 = K * sin(P);
            py2 = L * sin(Q);
            pz1 = K * sinh(N);
            pz2 = L * sinh(O);
            mm = sqrt((En * 1000 - wr0 - wr1) * (En * 1000 - wr0 - wr1) - (px1 + px2) * (px1 + px2) - (py1 + py2) * (py1 + py2) - (pz1 + pz2) * (pz1 + pz2));// missing mass
            M = sqrt((E - F) * (E - F) + (G - H) * (G - H));// Missing ET - Z pt
            //R = sqrt((E - Mte[0]) * (E - Mte[0]) + (G - Mtp[0]) * (G - Mtp[0]));
            //S = sqrt((F - Mte[0]) * (F - Mte[0]) + (H - Mtp[0]) * (H - Mtp[0]));
            J = sqrt(2 * A * B * (cosh(E - F) - cos(G - H)));// reconstructed Z invariant mass before cut
            myhdptmz->Fill(M);
            myhmet->Fill(Met[0]);
            myhzmu->Fill(J);
            if(J < 120 && J > 60)// invariant mass cut
            {
              for(Int_t f = 0; f < *Pts; f++)// loop over generator level particles
              {
                U = ppid[f];// loop over particle pid
                V = sta[f];// loop over particle status
                if(V == 1 && (U == 11 || U == -11 || U == 13 || U == -13))// find a muon or electron in final state
                {
                  if(U < 0)// antiparticle, particle 1
                  {
                    Px1 = ppx[f];
                    Py1 = ppy[f];
                    Pz1 = ppz[f];
                    Pe1 = pe[f];
                  }
                  if(U > 0)// particle 2
                  {
                    Px2 = ppx[f];
                    Py2 = ppy[f];
                    Pz2 = ppz[f];
                    Pe2 = pe[f];
                  }
                }
              }
              if(*Jts == 2)// 2 jets detected
              {
                Float_t Jtpt1 = jtpt[0];
                Float_t Jtpt2 = jtpt[1];
                Z = Jtpt1 + Jtpt2;// total jet pt
              }
              else
              {
                Z = 0;
              }
              myhjpt->Fill(Z);
              if(*Fms > 0)// forward muon detected
              {
                W = 0;
                for(Int_t i = 0; i < *Fms; i++)
                {
                  Int_t aa = fmpt[i];
                  Int_t ab = fme[i];
                  Int_t ac = fmph[i];
                  W = W + sqrt(aa * sin(ac) * aa * sin(ac) + aa * cos(ac) * aa * cos(ac) + aa * sinh(ab) * aa * sinh(ab));// total forward muon pt
                }
              }
              else
              {
                W = 0;
              }
              if(*Phs > 0)// photon detected
              {
                X = 0;
                for(Int_t i = 0; i < *Phs; i++)
                {
                  Int_t aa = phe[i];
                  X = X + aa;// total photon energy
                }
              }
              else
              {
                X = 0;
              }
              Y = Mtp[0] - atan2(py1 + py2, px1 + px2);// delta phi between missingET and Z
              Float_t gmm = sqrt((En * 1000 - Pe1 - Pe2) * (En * 1000 - Pe1 - Pe2) - (Px1 + Px2) * (Px1 + Px2) - (Py1 + Py2) * (Py1 + Py2) - (Pz1 + Pz2) * (Pz1 + Pz2));// gen level missing mass
              myhnop->Fill(*Efts + *Efphs + *Efnhs);// number of particle detected
              myhzm->Fill(J);// Z mass distribution after cut
              myhlpt1->Fill(K);// leading particle pt
              if(Z < 1200 - 800 / (En - 1))// jet pt cut
              {
                myhdfmz->Fill(Y);// dphi after cut
                myhee->Fill(X + W);// extra energy 
                if((Y > 3 && Y < 3.28) || (-3.28 < Y && Y < -3))// dphi cut
                {
                  myhee2->Fill(X + W);// extra energy
                  if((X + W) < En * 50)// extra energy cut
                  {
                    if(Met[0] < 12000)// missingET cut
                    {
                      myhm2->Fill(gmm);// gen missingm 
                      myhm1->Fill(mm);// missingm
                      if(mm > 15000)// missing mass cut
                      {
                        myhmm->Fill(mm, Met[0]);// 2D histogram fill(need to change)
                      }
                    }
                    if(mm > 15000)
                    {
                      myhmetg->Fill(Met[0]);// missingET after missingm cut
                    }
                    myhmm1->Fill(mm, Met[0]);
                    myhht->Fill(ht[0]);// ht 
                  }
                }
              }
              myhlpt2->Fill(L);// sub leading 
              myhmetph->Fill(Mtp[0]);// missingET phi
            }



        }
    }
  };



   // Launch the parallel processing of the tree
                    tp.Process(myFunction);

  // Use the TThreadedObject::Merge method to merge the thread private histograms
  // into the final result
  auto ptHistMerged = hmet0.Merge();
  auto ptHistMerged1 = hm1.Merge();
  auto ptHistMerged2 = hm2.Merge();
  auto ptHistMerged3 = hht.Merge();
  auto ptHistMerged4 = hjpt.Merge();
  auto ptHistMerged5 = hlpt1.Merge();
  auto ptHistMerged6 = hlpt2.Merge();
  auto ptHistMerged7 = hmetg.Merge();
  auto ptHistMerged8 = hee.Merge();
  auto ptHistMerged9 = hzmu.Merge();
  auto ptHistMerged10 = hzm.Merge();
  auto ptHistMerged11 = hdptmz.Merge();
  auto ptHistMerged12 = hmetph.Merge();
  auto ptHistMerged13 = hnop.Merge();
  auto ptHistMerged14 = hdfmz.Merge();
  auto ptHistMerged15 = hee2.Merge();
  auto ptHistMerged16 = hmet.Merge();
  auto ptHistMerged17 = hmm.Merge();
  auto ptHistMerged18 = hmm1.Merge();
  end.Print();// Timer

  // write histograms

  TH1F* h = ptHistMerged.get();
  TH1F *h1 = ptHistMerged1.get();
  TH1F *h2 = ptHistMerged2.get();
  TH1F *h3 = ptHistMerged3.get();
  TH1F *h4 = ptHistMerged4.get();
  TH1F *h5 = ptHistMerged5.get();
  TH1F *h6 = ptHistMerged6.get();
  TH1F *h7 = ptHistMerged7.get();
  TH1F *h8 = ptHistMerged8.get();
  TH1F *h9 = ptHistMerged9.get();
  TH1F *h10 = ptHistMerged10.get();
  TH1F *h11 = ptHistMerged11.get();
  TH1F *h12 = ptHistMerged12.get();
  TH1F *h13 = ptHistMerged13.get();
  TH1F *h14 = ptHistMerged14.get();
  TH1F *h15 = ptHistMerged15.get();
  TH1F *h16 = ptHistMerged16.get();
  TH2F *h17 = ptHistMerged17.get();
  TH2F *h18 = ptHistMerged18.get();
  TFile ff("rs/h/mchi1.root", "UPDATE");
  h->Write(fil);
  h1->Write(fil + "1");
  h2->Write(fil + "2");
  h3->Write(fil + "3");
  h16->Write(fil + "4");
  h4->Write(fil + "5");
  h5->Write(fil + "6");
  h6->Write(fil + "7");
  h7->Write(fil + "8");
  h8->Write(fil + "9");
  h9->Write(fil + "10");
  h10->Write(fil + "11");
  h11->Write(fil + "12");
  h12->Write(fil + "13");
  h13->Write(fil + "14");
  h14->Write(fil + "15");
  h15->Write(fil + "16");
  h17->Write(fil + "17");
  h18->Write(fil + "18");
  ff.Close();
}
