void sliceFit() 
{
    gStyle->SetOptFit(1111);
    gStyle->SetHistLineWidth(2);
    gStyle->SetOptStat(0);
    TCanvas *c1 = new TCanvas("c1","c1",800,600);
    TFile *f1 = new TFile("rpf2D_alphabet_M4bcut.root");
    //TFile *f1 = new TFile("rpf2D_data_HYto4b.root");
    TH2F* rpf_MhMA0 = (TH2F*)f1->Get("h_data_MhMA0_rpf");
    TH2F* pass_MhMA0 = (TH2F*)f1->Get("h_data_MhMA0_pass");
    TH2F* fail_MhMA0 = (TH2F*)f1->Get("h_data_MhMA0_fail");
    TF1 *fpol = new TF1("f1","pol3(0)",rpf_MhMA0->GetXaxis()->GetXmin(),rpf_MhMA0->GetXaxis()->GetXmax());
    int nSliceBin = 2;
    int nBinY = rpf_MhMA0->GetNbinsY(); 
    TH1F *h_sliceProj_pass[10], *h_sliceProj_fail[10];
    TH1F *h_sliceProj_rpf[10];
    TPaveStats *st;
    string pdfName = Form("sliceFit_%dbins.pdf",nSliceBin);
    double npar[20][10];
    c1->Print((pdfName+"[").data());
    for (int i=1;i<=nBinY-nSliceBin+1;i++) {
        //string namep = Form("h_pass%d",i);
        //string namef = Form("h_fail%d",i);
        string namep = Form("%s_MA0%dto%d",pass_MhMA0->GetName(),(int)rpf_MhMA0->GetYaxis()->GetBinLowEdge(i),(int)rpf_MhMA0->GetYaxis()->GetBinLowEdge(i+nSliceBin));
        string namef = Form("Alphabet_MA0%dto%d",(int)rpf_MhMA0->GetYaxis()->GetBinLowEdge(i),(int)rpf_MhMA0->GetYaxis()->GetBinLowEdge(i+nSliceBin));
        string name = Form("%s_MA0%dto%d",rpf_MhMA0->GetName(),(int)rpf_MhMA0->GetYaxis()->GetBinLowEdge(i),(int)rpf_MhMA0->GetYaxis()->GetBinLowEdge(i+nSliceBin));
        h_sliceProj_pass[i-1] = (TH1F*)pass_MhMA0->ProjectionX(namep.data(),i,i+nSliceBin-1,"eo");
        h_sliceProj_fail[i-1] = (TH1F*)fail_MhMA0->ProjectionX(namef.data(),i,i+nSliceBin-1,"eo");
        h_sliceProj_rpf[i-1] = (TH1F*)h_sliceProj_pass[i-1]->Clone(); 
        h_sliceProj_rpf[i-1]->Divide(h_sliceProj_pass[i-1],h_sliceProj_fail[i-1]);
        h_sliceProj_rpf[i-1]->SetTitle(name.data());
        h_sliceProj_rpf[i-1]->SetName(name.data());
        h_sliceProj_rpf[i-1]->Fit(fpol,"IMF");
        h_sliceProj_rpf[i-1]->SetMaximum(h_sliceProj_rpf[i-1]->GetMaximum()*1.8);
        /*
        h_sliceProj_pass[i-1]->Divide(h_sliceProj_fail[i-1]);
        h_sliceProj_pass[i-1]->SetTitle(name.data());
        h_sliceProj_pass[i-1]->SetName(name.data());
        h_sliceProj_pass[i-1]->Fit(fpol);
        h_sliceProj_pass[i-1]->SetMaximum(h_sliceProj_pass[i-1]->GetMaximum()*1.8);
        */
        //h_sliceProj_pass->Draw();
        gPad->Update();
        st = (TPaveStats*)h_sliceProj_rpf[i-1]->GetListOfFunctions()->FindObject("stats");
        //st->SetY2NDC(0.5);

        c1->Print(pdfName.data());
        for (int j=0;j<fpol->GetNpar();j++) npar[i-1][j] = fpol->GetParameter(j);
    }
    for (int i=0;i<=nBinY-nSliceBin;i++) {
        h_sliceProj_pass[i]->Draw("e");
        fpol->SetParameters(npar[i]);
        h_sliceProj_fail[i]->Multiply(fpol);
        h_sliceProj_pass[i]->SetLineColor(kBlack);
        h_sliceProj_fail[i]->SetLineColor(kBlue);
        h_sliceProj_fail[i]->Draw("esame");
        c1->BuildLegend();
        c1->Print(pdfName.data());
    }
    c1->Print((pdfName+"]").data());
    for (int i=0;i<=nBinY-nSliceBin;i++) {
        for (int j=0;j<fpol->GetNpar();j++) cout << npar[i][j] << "\t";
        cout << endl;
    }
    TObjArray *aSlices = new TObjArray();
    aSlices->SetOwner(kTRUE);
    //rpf_MhMA0->FitSlicesY(fpol,0,-1,0,"QNSG2",aSlices);
    printf("\n");
    //aSlices->Print();

}
