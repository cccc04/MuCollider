
{
  TH1F work(TString vfname, Double_t N); // 1d histogram scaling
  Int_t gn(Int_t Ev);   // function to correlate Ecm with the right column of cross section in the file
  TH1F work1(TString vfname, Double_t N); // unused function
  TH2F work2(TString vfname, Double_t N); // 2d histogram scaling


  Int_t E = 30; // Ecm
  Int_t et = gn(E); 
  Float_t lumi = gnn(E) * 1000000; // luminosity
  TString as = "25";
  TString aas = "2_5";
  TString ee =  to_string(E);
  TString st =  to_string(E) + "18;1"; // the name of histograms to be read; After cut
  TString stt = to_string(E) + "17;1"; // the name of histograms to be read; Before cut
  TTree *t = new TTree("t", "Book");
  // read cross section data
  t->ReadFile("Book"+aas+".csv", "ZZ:DM:WWZ:WWZL:ZZA:ZZAL:ZZZ:WWA");
  TLeaf *zz = t->GetLeaf("ZZ");
  TLeaf *wwz = t->GetLeaf("WWZ");
  TLeaf *dm = t->GetLeaf("DM");
  TLeaf *wwzl = t->GetLeaf("WWZL");
  TLeaf *zza = t->GetLeaf("ZZA");
  TLeaf *zzal = t->GetLeaf("ZZAL");
  TLeaf *zzz = t->GetLeaf("ZZZ");
  TLeaf *wwa = t->GetLeaf("WWA");
  t->GetEntry(et);// get the correct column 
  auto a1 = zz->GetValue();
  auto a2 = wwz->GetValue();
  auto a3 = dm->GetValue();
  auto a4 = wwzl->GetValue();
  auto a5 = zza->GetValue();
  auto a6 = zzal->GetValue();
  auto a7 = zzz->GetValue();
  auto a8 = wwa->GetValue();
  //Scale factor calculation
  auto N1 = a1 * lumi / 1000000 ;
  auto N2 = a2 *lumi / 1000000;
  auto N3 = a3 *lumi / 10000 ;
  auto N4 = a4 *lumi / 1000000;
  auto N5 = a5 *lumi / 1000000;
  auto N6 = a6 *lumi / 1000000;
  auto N7 = a7 *lumi / 1000000;
  auto N8 = a8 *lumi / 1000000;
  TCanvas *c1;
  c1 = new TCanvas("c1", "Plots");
  c1->SetLogz(1);
  TH2F *h1 = new TH2F();
  *h1 = work2("dm"+as+st, N3);
  TH2F *h2 = new TH2F();
  *h2 = work2("zvv" + st, N1);
  TH2F *h3 = new TH2F();
  *h3 = work2("wwz" + st, N4);
  TH2F *h6 = new TH2F();
  //*h6 = work("zzal" + st, N6);
  TH2F *h4 = new TH2F();
  //*h4 = work("wwzl" + st, N4);
  TH2F *h5 = new TH2F();
  *h5 = work2("zza" + st, N6);
  TH2F *h7 = new TH2F();
  *h7 = work2("zzz" + st, N7);
  TH2F *h8 = new TH2F();
  *h8 = work2("wwa" + st, N8);
  TH2F *h9 = new TH2F();
  *h9 = work2("dm"+ as + stt, N3);
  TH2F *h10 = new TH2F();
  *h10 = work2("zvv" + stt, N1);
  TH2F *h11 = new TH2F();
  *h11 = work2("wwz" + stt, N4);
  TH2F *h12 = new TH2F();
  *h12 = work2("zza" + stt, N6);
  TH2F *h13 = new TH2F();
  *h13 = work2("zzz" + stt, N7);
  TH2F *h14 = new TH2F();
  *h14 = work2("wwa" + stt, N8);
  // histogram style settings
  h1->SetName("h1");
  h2->SetName("h2");
  h3->SetName("h3");
  h4->SetName("h4");
  h5->SetName("h5");
  h6->SetName("h6");
  h7->SetName("h7");
  h8->SetName("h8");
  h1->SetMarkerStyle(21);
  h2->SetMarkerStyle(21);
  h3->SetMarkerStyle(21);
  h4->SetMarkerStyle(21);
  h5->SetMarkerStyle(21);
  h6->SetMarkerStyle(21);
  h7->SetMarkerStyle(21);
  h8->SetMarkerStyle(21);
  h1->SetLineColorAlpha(kGray + 2, 1);
  h2->SetLineColorAlpha(kGray + 2, 1);
  h3->SetLineColorAlpha(kGray + 2, 1);
  h4->SetLineColorAlpha(kGray + 2, 1);
  h5->SetLineColorAlpha(kGray + 2, 1);
  h6->SetLineColorAlpha(kGray + 2, 1);
  h7->SetLineColorAlpha(kGray + 2, 1);
  h8->SetLineColorAlpha(kGray + 2, 1);
  h1->SetFillColorAlpha(kGreen,0.7);
  h2->SetFillColorAlpha(kPink, 0.7);
  h3->SetFillColorAlpha(kYellow, 0.7);
  h4->SetFillColorAlpha(kYellow, 0.7);
  h5->SetFillColorAlpha(920, 0.7);
  h6->SetFillColorAlpha(920, 0.7);
  h7->SetFillColorAlpha(kAzure, 0.7);
  h8->SetFillColorAlpha(kViolet, 0.7);
  /* h1->SetFillStyle(3001);
  h2->SetFillStyle(3001);
  h3->SetFillStyle(3001);
  h4->SetFillStyle(3001);
  h5->SetFillStyle(3001);
  h6->SetFillStyle(3001);
  h7->SetFillStyle(3001);
  h8->SetFillStyle(3001);*/
  THStack *hs = new THStack("hs", "#sqrt{s} = "+ee+"TeV");// stack histogram creation (not used in this case for the bug that creates colored bins even if they are empty)
  hs->Add(h2);
  hs->Add(h3);
  //hs->Add(h4);
  hs->Add(h5);
  //hs->Add(h6);
  hs->Add(h7);
  hs->Add(h8);
  // histogram range
  Float_t min = 0.001;
  Float_t max = 100000;
  hs->SetMinimum(min);
  h1->SetMinimum(min);
  /**/h2->SetMinimum(min);
  h3->SetMinimum(min);
  h4->SetMinimum(min);
  h5->SetMinimum(min);
  h6->SetMinimum(min);
  h7->SetMinimum(min);
  h8->SetMinimum(min);
  hs->SetMaximum(max);
  h1->SetMaximum(max);
  /**/h2->SetMaximum(max);
  h3->SetMaximum(max);
  h4->SetMaximum(max);
  h5->SetMaximum(max);
  h6->SetMaximum(max);
  h7->SetMaximum(max);
  h8->SetMaximum(max);
  // legend
  TLegend *l = new TLegend(0.1, 0.1, 0.9, 0.9);
  gStyle->SetLegendTextSize(.05);
  l->AddEntry(h1, "#it{#chi#bar{#chi}Z : #mu^{+}#mu^{-} > #chi_{n}#bar{#chi}_{n}Z , Z > l^{+}l^{-}}", "f");
  l->AddEntry(h2, "#it{Z#nu#nu : #mu^{+}#mu^{-} > Z#nu_{l}#bar{#nu}_{l} , Z > l^{+}l^{-}}", "f");
  //l->AddEntry(h3, "wwz :", "f");
  l->AddEntry(h3, "#it{WWZ : #mu^{+}#mu^{-} > W^{+}W^{-}Z , W^{+} > l^{+}#nu_{l} , W^{-} > l^{-}#bar{#nu}_{l} , Z > #nu_{l}#bar{#nu}_{l}}", "f");
  l->AddEntry(h5, "#it{ZZ#gamma : #mu^{+}#mu^{-} > ZZ#gamma , Z > l^{+}l^{-} , Z > #nu_{l}#bar{#nu}_{l}}", "f");
 // l->AddEntry(h5, "zza'", "f");
  l->AddEntry(h7, "#it{ZZZ : #mu^{+}#mu^{-} > ZZZ, Z > l^{+}l^{-} , Z > #nu_{l}#bar{#nu}_{l} , Z > #nu_{l}#bar{#nu}_{l}}", "f");
  l->AddEntry(h8, "#it{ww#gamma : #mu^{+}#mu^{-} > W^{+}W^{-}#gamma ,  W^{+} > l^{+}#nu_{l} , W^{-} > l^{-}#bar{#nu}_{l}}", "f");
  
  // draw histograms
  
  //hs->Draw("SAME LEGO4");
  TString opt = "LEGO40";
  /**/h2->Draw("SAME "+opt);
  h3->Draw("SAME " + opt);
  h5->Draw("SAME " + opt);
  h7->Draw("SAME " + opt);
  h8->Draw("SAME " + opt);
  h2->SetTitle("#sqrt{s} = " + ee + "TeV, m_{#chi} = 0.8 TeV, L = " + lumi / 1000000 + "ab^{-1}"); 
  h2->SetTitleSize(0.03); 
  h1->Draw("SAME LEGO40");
  h2->GetXaxis()->SetTitle("m_{missing}[GeV]");
  h2->GetYaxis()->SetTitle("MET[GeV]");
  h2->GetXaxis()->SetTitleOffset(1.5);
  h2->GetYaxis()->SetTitleOffset(1.5);
  h2->GetXaxis()->SetLabelSize(0.03);
  h2->GetYaxis()->SetLabelSize(0.03);
  h2->GetXaxis()->SetTitleSize(0.03);
  h2->GetYaxis()->SetTitleSize(0.03);
  h2->GetXaxis()->SetTickSize(0.02);
  h2->GetYaxis()->SetTickSize(0.02);
  Int_t a = h1->GetEntries();
  Int_t b = h2->GetEntries();
  Int_t c = h3->GetEntries();
  //h4->GetEntries();
  Int_t d = h5->GetEntries();
  //h6->GetEntries();
  Int_t e = h7->GetEntries();
  Int_t f = h8->GetEntries();
  Int_t g = h9->GetEntries();
  Int_t h = h10->GetEntries();
  Int_t j = h11->GetEntries();
  Int_t k = h12->GetEntries();
  Int_t n = h13->GetEntries();
  Int_t m = h14->GetEntries();
  // signal strength calculation
  Int_t ss = a *N3;
  Int_t bb = b * N1 + c * N4 + d * N6 + e * N7 + f * N8;
  Int_t ss1 = g * N3;
  Int_t bb1 = h * N1 + j * N4 + k * N6 + n * N7 + m * N8;
  Float_t ss2 = ss/sqrt(bb);
  Float_t bb2 = ss1 / sqrt(bb1);
  //auto sb = ss / bb ^ (0.5);
  TString tt1 = "Before: s = " + to_string(ss);
  TString tt2 = "b = " + to_string(bb);
  TString tt3 = "#frac{s}{#sqrt{b}} = " + to_string(ss2);
  TString tt4 = "After: s = " + to_string(ss1);
  TString tt5 = "b = " + to_string(bb1);
  TString tt6 = " #frac{s}{#sqrt{b}} = " + to_string(bb2);
  TLatex *t1 = new TLatex(0.36, 0.28, tt1);
  TLatex *t2 = new TLatex(0.43, 0.23, tt2);
  TLatex *t3 = new TLatex(0.41, 0.18, tt3);
  TLatex *t4 = new TLatex(0.38, 0.11, tt4);
  TLatex *t5 = new TLatex(0.43, 0.06, tt5);
  TLatex *t6 = new TLatex(0.41, 0.01, tt6);
  t1->SetTextFont(2);
  t1->SetTextSize(0.03);
  t2->SetTextFont(2);
  t2->SetTextSize(0.03);
  t3->SetTextFont(2);
  t3->SetTextSize(0.03);
  t4->SetTextFont(2);
  t4->SetTextSize(0.03);
  t5->SetTextFont(2);
  t5->SetTextSize(0.03);
  t6->SetTextFont(2);
  t6->SetTextSize(0.03);
  t1->Draw();
  t2->Draw();
  t3->Draw();
  t4->Draw();
  t5->Draw();
  t6->Draw();
  c1->Modified();
  //TCanvas c2;
  //l->Draw();
}

TH1F work(TString vfname, Double_t N)
{
  TH1F *h;
  TFile *f = new TFile("h/hist12.root");
  f->GetObject(vfname, h);
  h->Scale(N);
  return *h;
}

TH2F work2(TString vfname, Double_t N)
{
  TH2F *h;
  TFile *f = new TFile("h/mchi"+aas+".root");
  f->GetObject(vfname, h);
  h->Scale(N);
  return *h;
}

Int_t gn(Int_t Ev)
{
  Int_t n;
  if(Ev == 1)
  {
    n = 0;
    return n;
  }
  if(Ev == 2)
  {
    n = 1;
    return n;
  }
  if(Ev == 3)
  {
    n = 2;
    return n;
  }
  if(Ev == 6)
  {
    n = 3;
    return n;
  }
  if(Ev == 10)
  {
    n = 4;
    return n;
  }
  if(Ev == 15)
  {
    n = 5;
    return n;
  }
  if(Ev == 30)
  {
    n = 6;
    return n;
  }
}
Int_t gnn(Int_t Ev)
{
  Int_t n;
  if(Ev == 3)
  {
    n = 1;
    return n;
  }
  if(Ev == 6)
  {
    n = 4;
    return n;
  }
  if(Ev == 10)
  {
    n = 10;
    return n;
  }

  if(Ev == 30)
  {
    n = 10;
    return n;
  }
}
TH1F work1(TString vfname, Double_t N)
{
  TH1F *h;
  TFile *f = new TFile(vfname);
  TTree *a = (TTree *)f->Get("Delphes");
  a->Draw("MissingET.MET");
  h = (TH1F *)gPad->GetPrimitive("htemp");
  h->Scale(N);
  return *h;
}