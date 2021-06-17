
{

  TString filei = "dm3";
  Int_t workers = 11;

  auto start = std::chrono::system_clock::now();
  TString filef = ".root";
  TString file = filei + filef;
  ROOT::TProcessExecutor pool;
  std::vector<int> gh;
  for(Int_t v = 0; v < workers; v++)
  {
    gh.push_back(v);
  }
  auto work = [workers, file](Int_t h) {
    TFile *f = new TFile(file);
    TTree *a = (TTree *)f->Get("Delphes");
    Int_t nEn = a->GetEntries("Events");
    Float_t A;
    Float_t B;
    Float_t C;
    Float_t E;
    Float_t F;
    Float_t G;
    Float_t H;
    Float_t J;
    Int_t D = floor(nEn/ workers);
    TH1F* hmet = new TH1F("MET", "MissingET", 250, 0, 1500);
    TH1F* hmte = new TH1F("Eta", "MissingET", 250, -8, 8);
    TH1F* hmtp = new TH1F("Phi", "MissingET", 250, -4, 4);
    TH1F *hzpt = new TH1F("PT", "Reconstructed Z", 250, 0, 1500);
    TH1F *hzm = new TH1F("Mass", "Reconstructed Z", 250, 0, 500);
    for(Int_t e = h * D; e < D * (h + 1); e++)
    {
        a->GetEntry(e);
        TLeaf* ep = a->GetLeaf("Electron.PT");
        TLeaf* ee = a->GetLeaf("Electron.Eta");
        TLeaf* eph = a->GetLeaf("Electron.Phi");
        TLeaf* mp = a->GetLeaf("Muon.PT");
        TLeaf* me = a->GetLeaf("Muon.Eta");
        TLeaf* mph = a->GetLeaf("Muon.Phi");
        TLeaf* fp = a->GetLeaf("ForwardMuon.PT");
        TLeaf* fe = a->GetLeaf("ForwardMuon.Eta");
        TLeaf* fph = a->GetLeaf("ForwardMuon.Phi");
        TLeaf* es = a->GetLeaf("Electron_size");
        TLeaf* ms = a->GetLeaf("Muon_size");
        TLeaf* fs = a->GetLeaf("ForwardMuon_size");
        TLeaf* met = a->GetLeaf("MissingET.MET");
        TLeaf* mte = a->GetLeaf("MissingET.Eta");
        TLeaf* mtp = a->GetLeaf("MissingET.Phi");
        Int_t Es = es->GetValue();
        Int_t Ms = ms->GetValue();
        Int_t Fs = fs->GetValue();
        Float_t Met;
        Float_t Mte;
        Float_t Mtp;
        if (Es == 2)
        {
            A = ep->GetValue(0);
            B = ep->GetValue(1);
            C = A + B;
            E = ee->GetValue(0);
            F = ee->GetValue(1);
            G = eph->GetValue(0);
            H = eph->GetValue(1);
            J = sqrt(2 * A * B * (cosh(E - F) - cos(G - H)));
            if (J < 110 && J > 70) {
                Met = met->GetValue();
                Mte = mte->GetValue();
                Mtp = mtp->GetValue();
                hmte->Fill(Mte);
                hmtp->Fill(Mtp);
                hmet->Fill(Met);
                hzpt->Fill(C);
                hzm->Fill(J);
            }
        }
        else if (Ms == 2)
        {
            A = mp->GetValue(0);
            B = mp->GetValue(1);
            C = A + B;
            E = me->GetValue(0);
            F = me->GetValue(1);
            G = mph->GetValue(0);
            H = mph->GetValue(1);
            J = sqrt(2 * A * B * (cosh(E - F) - cos(G - H)));
            if (J < 110 && J > 70) {
                Met = met->GetValue();
                Mte = mte->GetValue();
                Mtp = mtp->GetValue();
                hmte->Fill(Mte);
                hmtp->Fill(Mtp);
                hmet->Fill(Met);
                hzpt->Fill(C);
                hzm->Fill(J);
            }
        }
        else if (Fs == 2)
        {
            A = fp->GetValue(0);
            B = fp->GetValue(1);
            C = A + B;
            E = fe->GetValue(0);
            F = fe->GetValue(1);
            G = fph->GetValue(0);
            H = fph->GetValue(1);
            J = sqrt(2 * A * B * (cosh(E - F) - cos(G - H)));
            if (J < 110 && J > 70) {
                Met = met->GetValue();
                Mte = mte->GetValue();
                Mtp = mtp->GetValue();
                hmte->Fill(Mte);
                hmtp->Fill(Mtp);
                hmet->Fill(Met);
                hzpt->Fill(C);
                hzm->Fill(J);
            }
        }
        else if (Ms == 1 && Fs == 1)
        {
            A = mp->GetValue();
            B = fp->GetValue();
            C = A + B;
            E = me->GetValue();
            F = fe->GetValue();
            G = mph->GetValue();
            H = fph->GetValue();
            J = sqrt(2 * A * B * (cosh(E - F) - cos(G - H)));
            if (J < 110 && J > 70) {
                Met = met->GetValue();
                Mte = mte->GetValue();
                Mtp = mtp->GetValue();
                hmte->Fill(Mte);
                hmtp->Fill(Mtp);
                hmet->Fill(Met);
                hzpt->Fill(C);
                hzm->Fill(J);
            }
        }
    }
    std::vector<TH1F*, std::allocator<TH1F*> > gg;
    gg = { hmet, hzpt, hzm,hmte,hmtp };
    return gg;
  };
  auto squares = pool.Map(work, gh);
  TFile *f = new TFile(file);
  TTree *a = (TTree *)f->Get("Delphes");
  Int_t nEn = a->GetEntries("Events");
  Float_t A;
  Float_t B;
  Float_t C;
  Float_t E;
  Float_t F;
  Float_t G;
  Float_t H;
  Float_t J;
  Int_t D = floor(nEn / workers);
  TH1F *hmet = new TH1F("MET", "MissingET", 250, 0, 1500);
  TH1F *hmte = new TH1F("Eta", "MissingET", 250, -8, 8);
  TH1F *hmtp = new TH1F("Phi", "MissingET", 250, -4,4);
  TH1F *hzpt = new TH1F("PT", "Reconstructed Z", 250, 0, 1500);
  TH1F *hzm = new TH1F("Mass", "Reconstructed Z", 250, 0, 500);
  for(Int_t e = D * workers; e < nEn; e++)
  {
    a->GetEntry(e);
    TLeaf *ep = a->GetLeaf("Electron.PT");
    TLeaf *ee = a->GetLeaf("Electron.Eta");
    TLeaf *eph = a->GetLeaf("Electron.Phi");
    TLeaf *mp = a->GetLeaf("Muon.PT");
    TLeaf *me = a->GetLeaf("Muon.Eta");
    TLeaf *mph = a->GetLeaf("Muon.Phi");
    TLeaf *fp = a->GetLeaf("ForwardMuon.PT");
    TLeaf *fe = a->GetLeaf("ForwardMuon.Eta");
    TLeaf *fph = a->GetLeaf("ForwardMuon.Phi");
    TLeaf *es = a->GetLeaf("Electron_size");
    TLeaf *ms = a->GetLeaf("Muon_size");
    TLeaf *fs = a->GetLeaf("ForwardMuon_size");
    TLeaf *met = a->GetLeaf("MissingET.MET");
    TLeaf *mte = a->GetLeaf("MissingET.Eta");
    TLeaf *mtp = a->GetLeaf("MissingET.Phi");
    Int_t Es = es->GetValue();
    Int_t Ms = ms->GetValue();
    Int_t Fs = fs->GetValue();
    Float_t Met;
    Float_t Mte;
    Float_t Mtp;
    if (Es == 2)
    {
        A = ep->GetValue(0);
        B = ep->GetValue(1);
        C = A + B;
        E = ee->GetValue(0);
        F = ee->GetValue(1);
        G = eph->GetValue(0);
        H = eph->GetValue(1);
        J = sqrt(2 * A * B * (cosh(E - F) - cos(G - H)));
        if (J < 110 && J > 70) {
            Met = met->GetValue();
            Mte = mte->GetValue();
            Mtp = mtp->GetValue();
            hmte->Fill(Mte);
            hmtp->Fill(Mtp);
            hmet->Fill(Met);
            hzpt->Fill(C);
            hzm->Fill(J);
        }
    }
    else if(Ms == 2)
    {
      A = mp->GetValue(0);
      B = mp->GetValue(1);
      C = A + B;
      E = me->GetValue(0);
      F = me->GetValue(1);
      G = mph->GetValue(0);
      H = mph->GetValue(1);
      J = sqrt(2 * A * B * (cosh(E - F) - cos(G - H)));
      if (J < 110 && J > 70) {
          Met = met->GetValue();
          Mte = mte->GetValue();
          Mtp = mtp->GetValue();
          hmte->Fill(Mte);
          hmtp->Fill(Mtp);
          hmet->Fill(Met);
          hzpt->Fill(C);
          hzm->Fill(J);
      }
    }
    else if(Fs == 2)
    {
      A = fp->GetValue(0);
      B = fp->GetValue(1);
      C = A + B;
      E = fe->GetValue(0);
      F = fe->GetValue(1);
      G = fph->GetValue(0);
      H = fph->GetValue(1);
      J = sqrt(2 * A * B * (cosh(E - F) - cos(G - H)));
      if (J < 110 && J > 70) {
          Met = met->GetValue();
          Mte = mte->GetValue();
          Mtp = mtp->GetValue();
          hmte->Fill(Mte);
          hmtp->Fill(Mtp);
          hmet->Fill(Met);
          hzpt->Fill(C);
          hzm->Fill(J);
      }
    }
    else if(Ms == 1 && Fs == 1)
    {
      A = mp->GetValue();
      B = fp->GetValue();
      C = A + B;
      E = me->GetValue();
      F = fe->GetValue();
      G = mph->GetValue();
      H = fph->GetValue();
      J = sqrt(2 * A * B * (cosh(E - F) - cos(G - H)));
      if (J < 110 && J > 70) {
          Met = met->GetValue();
          Mte = mte->GetValue();
          Mtp = mtp->GetValue();
          hmte->Fill(Mte);
          hmtp->Fill(Mtp);
          hmet->Fill(Met);
          hzpt->Fill(C);
          hzm->Fill(J);
      }
    }
  }
  std::vector<TH1F *, std::allocator<TH1F *> > gg;
  gg = {hmet, hzpt, hzm,hmte,hmtp};

  c1 = new TCanvas("c1", "Plots", 1920, 1080);
  c1->Divide(3, 3);

  for(Int_t i = 0; i < 5; i++)
  {
    squares[0][i]->Add(gg[i]);
    for(Int_t ii = 1; ii < workers; ii++)
    {
      squares[0][i]->Add(squares[ii][i]);
    }
    c1->cd(i + 1);
    squares[0][i]->Draw();
  }
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> diff = end - start;
  std::cout << "Time = " << diff.count() << " s\n";
  f->Close();
  TFile ff("h/his.root", "UPDATE");
  TH1F *h = squares[0][0];
  TH1F *h1 = squares[0][1];
  TH1F *h2 = squares[0][2];
  TH1F* h3 = squares[0][3];
  TH1F* h4 = squares[0][4];
  h->Write(filei);
  h1->Write(filei + "1");
  h2->Write(filei + "2");
  h3->Write(filei + "3");
  h4->Write(filei + "4");
  ff.Close();

}

