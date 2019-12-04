#define M4bCut 1
inline double m4bTilde(double M4b, double MA0, double Mh, double MA0slice) {
    return M4b - MA0 + MA0slice - Mh + 125.;
}
void sliceFitTest() 
{
    gStyle->SetOptFit(1111);
    gStyle->SetHistLineWidth(2);
    gStyle->SetOptStat(0);
    gStyle->SetMarkerStyle(0);
    TCanvas *c1 = new TCanvas("c1","c1",800,600);
    TFile *f1 = new TFile("rpf2D_alphabet_M4bcut.root");
    //TFile *f1 = new TFile("rpf2D_data_HYto4b.root");
    TFile *ftree = new TFile("ALLBTagCSV-Run2016_tree.root");
    TTree *t1 = (TTree*)ftree->Get("tree");
    TH2F* rpf_MhMA0 = (TH2F*)f1->Get("h_data_MhMA0_rpf");
    TH2F* pass_MhMA0 = (TH2F*)f1->Get("h_data_MhMA0_pass");
    TH2F* fail_MhMA0 = (TH2F*)f1->Get("h_data_MhMA0_fail");
    TF1 *fpol = new TF1("f1","pol3(0)",rpf_MhMA0->GetXaxis()->GetXmin(),rpf_MhMA0->GetXaxis()->GetXmax());
    int nSliceBin = 2;
    int nBinY = rpf_MhMA0->GetNbinsY(); 
    TH1F *h_sliceProj_pass, *h_sliceProj_fail;
    TH1F *h_sliceProj_rpf, *h_srRpf, *h_Rpf, *h_cr, *h_sb, *h_sr;
    TPaveStats *st;
    float minMA0 = 300;
    float maxMA0 = 400;
    string pdfName = Form("sliceTestFit_%dbins_MA0%dto%d.pdf",nSliceBin,(int)minMA0,(int)maxMA0);
    double npar[10];
    float mh,ma0,m4b,hpt;
    bool isTag, isAntiTag;
    t1->SetBranchAddress("Mh",&mh);
    t1->SetBranchAddress("MA0",&ma0);
    t1->SetBranchAddress("MZp",&m4b);
    t1->SetBranchAddress("hPt",&hpt);
    t1->SetBranchAddress("isTag",&isTag);
    t1->SetBranchAddress("isAntiTag",&isAntiTag);


    c1->Print((pdfName+"[").data());

        int i = 1;
        //string namep = Form("h_pass%d",i);
        //string namef = Form("h_fail%d",i);
        string namep = Form("%s_MA0%dto%d",pass_MhMA0->GetName(),(int)rpf_MhMA0->GetYaxis()->GetBinLowEdge(i),(int)rpf_MhMA0->GetYaxis()->GetBinLowEdge(i+nSliceBin));
        string namef = Form("Alphabet_MA0%dto%d",(int)rpf_MhMA0->GetYaxis()->GetBinLowEdge(i),(int)rpf_MhMA0->GetYaxis()->GetBinLowEdge(i+nSliceBin));
        //string name = Form("%s_MA0%dto%d",rpf_MhMA0->GetName(),(int)rpf_MhMA0->GetYaxis()->GetBinLowEdge(i),(int)rpf_MhMA0->GetYaxis()->GetBinLowEdge(i+nSliceBin));
        string name = Form("%s_MA0%dto%d",rpf_MhMA0->GetName(),(int)minMA0,(int)maxMA0);
        h_sliceProj_pass = (TH1F*)pass_MhMA0->ProjectionX(namep.data(),i,i+nSliceBin-1,"eo");
        h_sliceProj_fail = (TH1F*)fail_MhMA0->ProjectionX(namef.data(),i,i+nSliceBin-1,"eo");
        //pass_MhMA0->Draw();
        
        
        h_sliceProj_rpf = (TH1F*)h_sliceProj_pass->Clone(); 
        h_sliceProj_rpf->Divide(h_sliceProj_pass,h_sliceProj_fail);
        h_sliceProj_rpf->SetTitle(name.data());
        h_sliceProj_rpf->SetName(name.data());
        h_sliceProj_rpf->SetMaximum(h_sliceProj_rpf->GetMaximum()*1.8);
        string sCut = Form("isTag&&(MA0<%f)&&(MA0>%f)&&(Mh<136)&&(Mh>114)&MZp>400",maxMA0,minMA0);
        string sbCut = Form("isTag&&(MA0<%f)&&(MA0>%f)&&(Mh>136||Mh<114)&MZp>400",maxMA0,minMA0);
        t1->Draw("Mh>>h_sb(25,100,150)",sbCut.data(),"e");
        if (M4bCut) t1->Draw("Mh>>h_sr(25,100,150)",sCut.data(),"e");
        else t1->Draw("Mh>>h_sr(25,100,150)","isTag&&(MA0<300)&&(MA0>200)&&(Mh<136)&&(Mh>114)","e"); 
        //string sfCut = Form("isAntiTag&&(MA0<%f)&&(MA0>%f)&&(Mh<136)&&(Mh>114)&MZp>400",maxMA0,minMA0);
        string sfCut = Form("isAntiTag&&(MA0<%f)&&(MA0>%f)&&MZp>400",maxMA0,minMA0);
        if (M4bCut) t1->Draw("Mh>>h_cr(25,100,150)",sfCut.data(),"e");
        else t1->Draw("Mh>>h_cr(25,100,150)","isAntiTag&&(MA0<300)&&(MA0>200)&&(Mh<136)&&(Mh>114)","e"); 
        h_sr = (TH1F*)gDirectory->Get("h_sr");
        h_cr = (TH1F*)gDirectory->Get("h_cr");
        h_sb = (TH1F*)gDirectory->Get("h_sb");
        h_sr->Sumw2();
        h_sr->SetTitle("h_sr");
        h_sr->SetLineColor(kGray);
        h_sr->SetMarkerSize(0);
        h_srRpf = (TH1F*)h_sr->Clone("h_srRpf");
        h_Rpf = (TH1F*)h_sb->Clone("h_Rpf");
        h_srRpf->SetTitle("h_srRpf");
        h_Rpf->SetTitle("h_Rpf");
        h_srRpf->Sumw2();
        h_Rpf->Sumw2();
        h_srRpf->Divide(h_cr);
        h_Rpf->Divide(h_cr);
        //h_srRpf->Divide(h_sliceProj_fail);
        h_srRpf->SetMarkerSize(0);
        h_Rpf->SetMarkerSize(0);


        
        
        h_srRpf->SetLineColor(kGray+2);
        //h_sliceProj_rpf->Fit(fpol,"IMF");
        h_Rpf->SetMaximum(h_Rpf->GetMaximum()*1.5);
        h_Rpf->Fit(fpol,"IMF");
        h_srRpf->Draw("esame");
        c1->Print(pdfName.data());
        //h_srRpf->Draw("e");
        //c1->Print(pdfName.data());
        //h_sr->Draw("e");
        //c1->Print(pdfName.data());
        //h_sliceProj_pass->Draw();
        gPad->Update();
        st = (TPaveStats*)h_sliceProj_rpf->GetListOfFunctions()->FindObject("stats");
        //st->SetY2NDC(0.5);

        for (int j=0;j<fpol->GetNpar();j++) npar[j] = fpol->GetParameter(j);
        i = 0; 
        //h_sliceProj_pass->SetMaximum(h_sliceProj_pass->GetMaximum()*1.5);
        //h_sliceProj_pass->Draw("e");
        h_sb->SetMaximum(h_sb->GetMaximum()*1.5);
        h_sb->SetTitle("h_sideband");
        h_sb->Draw("e");
        h_sr->Draw("esame");
        fpol->SetParameters(npar);
        h_cr->Multiply(fpol);
        h_cr->SetTitle("h_failXrpf");
        //h_sliceProj_pass->SetLineColor(kBlack);
        h_sb->SetLineColor(kBlack);
        h_cr->SetLineColor(kBlue);
        h_cr->Draw("esame");
        c1->BuildLegend();
        c1->Print(pdfName.data());
        for (int j=0;j<fpol->GetNpar();j++) cout << npar[j] << "\t";
        cout << endl;
    
    TH1D *h_reduceMass_pass[3];
    TH1D *h_reduceMass_alphabet[3];
    TH1D *h_4bMass_pass[3];
    TH1D *h_4bMass_alphabet[3];
    cout << "nBin = " << nBinY << endl;
    string name2[2] = {"h_4bMass_pass_sb","h_4bMass_pass_sb+sr"};
    string name2Al[2] = {"h_4bMass_failXrpf_sb","h_4bMass_failXrpf_sb+sr"};
    string name3[2] = {"h_reduceMass_pass_sb","h_reduceMass_pass_sb+sr"};
    string name3Al[2] = {"h_reduceMass_failXrpf_sb","h_reduceMass_failXrpf_sb+sr"};
    for (int i=0;i<2;i++) {
        int ma0p = (int)minMA0;
        int width = 100;
        //int ma0p = (int)h_MhMA0_rpf->GetYaxis()->GetBinLowEdge(i+1);
        string namep = Form("h_reduceMass_pass_MA0%dto%d_sig%d",ma0p,(int)maxMA0,i);
        string nameAl = Form("h_reduceMass_Alphabet_MA0%dto%d_sig%d",ma0p,(int)maxMA0,i);
        string name4bp = Form("h_4bMass_pass_MA0%dto%d_sig%d",ma0p,(int)maxMA0,i);
        string name4bAl = Form("h_4bMass_Alphabet_MA0%dto%d_sig%d",ma0p,(int)maxMA0,i);
        h_reduceMass_pass[i] = new TH1D(namep.data(),name3[i].data(),48,400,1600);
        h_reduceMass_alphabet[i] = new TH1D(nameAl.data(),name3Al[i].data(),48,400,1600);
        h_reduceMass_alphabet[i]->SetLineColor(kRed+3*i);
        h_reduceMass_pass[i]->SetLineColor(kBlue+3*i);
        h_reduceMass_pass[i]->Sumw2();
        h_reduceMass_alphabet[i]->Sumw2();
        h_4bMass_pass[i] = new TH1D(name4bp.data(),name2[i].data(),48,400,1600);
        h_4bMass_alphabet[i] = new TH1D(name4bAl.data(),name2Al[i].data(),48,400,1600);
        h_4bMass_alphabet[i]->SetLineColor(kRed+3*i);
        h_4bMass_pass[i]->SetLineColor(kBlue+3*i);
        h_4bMass_pass[i]->Sumw2();
        h_4bMass_alphabet[i]->Sumw2();

    }
    /*
    h_4bMass_pass[1]->SetLineStyle(10);
    h_4bMass_alphabet[1]->SetLineStyle(10);
    h_reduceMass_pass[1]->SetLineStyle(10);
    h_reduceMass_alphabet[1]->SetLineStyle(10);
    */
    for (int i=0;i<t1->GetEntries();i++) {
        if (i%100000==0) cout << i << "/" << t1->GetEntries() << endl;
        t1->GetEntry(i);
        if (ma0>maxMA0||ma0<minMA0) continue; //[200,300]
        if (m4b<400) continue;
        if (mh>150||mh<100) continue;
        //if (mh>114&&mh<136) continue;
        //if (isTag&&(mh>114&&mh<136)) continue;
        double mtilde = m4bTilde(m4b,ma0,mh,300.);
        if (isTag) {
            h_reduceMass_pass[1]->Fill(mtilde);
            h_4bMass_pass[1]->Fill(m4b);
            if (mh<114||mh>136) {
                h_reduceMass_pass[0]->Fill(mtilde);
                h_4bMass_pass[0]->Fill(m4b);
            }
        }
        else if (isAntiTag){
            //fpol->SetParameters(nparList[j]); // set MA0 slice
            h_reduceMass_alphabet[1]->Fill(mtilde,fpol->Eval(mh));
            h_4bMass_alphabet[1]->Fill(m4b,fpol->Eval(mh));
            if (mh<114||mh>136) {
                h_reduceMass_alphabet[0]->Fill(mtilde,fpol->Eval(mh));
                h_4bMass_alphabet[0]->Fill(m4b,fpol->Eval(mh));
            }
        }
    }
    h_reduceMass_pass[0]->GetXaxis()->SetTitle("reduce mass (GeV)");
    h_reduceMass_pass[0]->GetYaxis()->SetTitle("Number of Events");
    h_reduceMass_pass[0]->Draw("e");
    h_reduceMass_pass[1]->Draw("esame");
    //for (int i=0;i<2;i++) h_reduceMass_pass[i]->Draw("esame");
    for (int i=0;i<2;i++) h_reduceMass_alphabet[i]->Draw("esame");
    c1->BuildLegend(0.5,0.6,0.9,0.9);
    c1->Print(pdfName.data());
    h_4bMass_pass[0]->GetXaxis()->SetTitle("reduce mass (GeV)");
    h_4bMass_pass[0]->GetYaxis()->SetTitle("Number of Events");
    h_4bMass_pass[0]->Draw("e");
    h_4bMass_pass[1]->Draw("esame");
    //for (int i=0;i<2;i++) h_4bMass_pass[i]->Draw("esame");
    for (int i=0;i<2;i++) h_4bMass_alphabet[i]->Draw("esame");
    c1->BuildLegend(0.5,0.6,0.9,0.9);
    c1->Print(pdfName.data());
    c1->Print((pdfName+"]").data());
    cout << "4bM pass: " << h_4bMass_pass[0]->Integral() << "\t4bM alphabet: " << h_4bMass_alphabet[0]->Integral() << endl;  
    cout << "reduce pass: " << h_reduceMass_pass[0]->Integral() << "\treduce alphabet: " << h_reduceMass_alphabet[0]->Integral() << endl;
    TFile *fwrt = new TFile("fwrt.root","recreate");
    h_srRpf->Write();
    fwrt->Close();

}
