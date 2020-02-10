#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TTree.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TPaveStats.h"
#define runShort 0
inline double m4bTilde(double M4b, double MA0, double Mh, double MA0slice) {
    return M4b - MA0 + MA0slice - Mh + 125.;
}
void sliceFitLoop(float hptStart = 0) 
{
    TRandom3 *gRandom = new TRandom3(0);
    gStyle->SetOptFit(1111);
    gStyle->SetHistLineWidth(2);
    gStyle->SetOptStat(0);
    gStyle->SetMarkerStyle(0);
    TH1::SetDefaultSumw2(1);
    TCanvas *c1 = new TCanvas("c1","c1",800,600);
    TFile *f1 = new TFile("../diff/rpf2D_alphabet_M4bcut.root");
    TFile *ftree = new TFile("BTagCSV-Run2016Full.root");
    TTree *t1 = (TTree*)ftree->Get("tree");
    //TTree *t2 = (TTree*)ftree->Get("tree2");
    //t1->AddFriend("tree2");
    TH2F* rpf_MhMA0 = (TH2F*)f1->Get("h_data_MhMA0_rpf");
    TH2F* pass_MhMA0 = (TH2F*)f1->Get("h_data_MhMA0_pass");
    TH2F* fail_MhMA0 = (TH2F*)f1->Get("h_data_MhMA0_fail");
    TF1 *fpol = new TF1("f1","pol3(0)",rpf_MhMA0->GetXaxis()->GetXmin(),rpf_MhMA0->GetXaxis()->GetXmax());
    int nSliceBin = 2;
    int nBinY = rpf_MhMA0->GetNbinsY(); 
    TH1F *h_srRpf, *h_Rpf, *h_cr, *h_sb, *h_sr;
    TH1F *h_hptPass, *h_hptFail, *h_ma0Pass, *h_ma0Fail;
    TH1F *h_dphiPass = new TH1F("h_dphi_pass","h_dphi_pass",32,0,3.14);
    TH1F *h_dphiFail = new TH1F("h_dphi_fail","h_dphi_fail",32,0,3.14);
    TH1F *h_dphiWeight = new TH1F("h_dphi_failXrpf","h_dphi_failXrpf",32,0,3.14);
    TH1F *h_dRPass = new TH1F("h_dR_pass","h_dR_pass",30,0,6);
    TH1F *h_dRFail = new TH1F("h_dR_fail","h_dR_fail",30,0,6);
    TH1F *h_dRWeight = new TH1F("h_dR_failXrpf","h_dR_failXrpf",30,0,6);
    TH1F *h_a0ptPass = new TH1F("h_A0Pt_pass","h_A0Pt_pass",50,0,500);
    TH1F *h_a0ptFail = new TH1F("h_A0Pt_fail","h_A0Pt_fail",50,0,500);
    TH1F *h_a0ptWeight = new TH1F("h_A0Pt_failXrpf","h_A0Pt_failXrpf",50,0,500);
    TLorentzVector* v1, *v2, *vh1, *vh2;
    TLorentzVector *va01, *va02;
    vh1 = new TLorentzVector();
    vh2 = new TLorentzVector();
    va01 = new TLorentzVector();
    va02 = new TLorentzVector();
    TPaveStats *st;
    double npar[10];
    float mh,ma0,m4b,hpt, m4b_loop, ma0_loop, mh_loop;
    bool isTag, isAntiTag, isHfail, isHpass; 
    t1->SetBranchAddress("Mh",&mh);
    t1->SetBranchAddress("MA0",&ma0);
    t1->SetBranchAddress("MZp",&m4b);
    t1->SetBranchAddress("hPt",&hpt);
    t1->SetBranchAddress("isTag",&isTag);
    t1->SetBranchAddress("isAntiTag",&isAntiTag);
    t1->SetBranchAddress("isHpass",&isHpass);
    t1->SetBranchAddress("isHfail",&isHfail);
    int nFailJetCom, nFailHCom;
    int failInd1[500];
    int failInd2[500];
    int HfailInd1[500];
    int HfailInd2[500];
    TClonesArray *arr = new TClonesArray("TLorentzVector");
    t1->SetBranchAddress("FailJetArr",&arr);
    t1->SetBranchAddress("higgs jet1",&vh1);
    t1->SetBranchAddress("higgs jet2",&vh2);
    t1->SetBranchAddress("A0 jet1",&va01);
    t1->SetBranchAddress("A0 jet2",&va02);
    t1->SetBranchAddress("nFailJetCom",&nFailJetCom);
    t1->SetBranchAddress("failInd1",failInd1);
    t1->SetBranchAddress("failInd2",failInd2);
    t1->SetBranchAddress("HfailInd1",HfailInd1);
    t1->SetBranchAddress("HfailInd2",HfailInd2);
    t1->GetBranch("FailJetArr")->SetAutoDelete(kFALSE);
    t1->SetBranchAddress("nFailHCom",&nFailHCom);

    //for (float hPtSlice=0;hPtSlice<350;hPtSlice+=50) {
        float hPtSlice = hptStart;
        float minhPt = hPtSlice;
        float maxhPt = hPtSlice+50;
        string pdfName = Form("sliceLoopRandFit_hPt%dto%d.pdf",(int)minhPt,(int)maxhPt);
        c1->Print((pdfName+"[").data());

        //string namep = Form("h_pass%d",i);
        int i=1;
        //pass_MhMA0->Draw();
        
        //  Draw histogram 
        //string sCut = Form("isTag&&(hPt<%f)&&(hPt>%f)&&(Mh<136)&&(Mh>114)&MZp>400",maxhPt,minhPt);
        string sbCut = Form("isTag&&(hPt<%f)&&(hPt>%f)&&(Mh>136||Mh<114)&MZp>400",maxhPt,minhPt);
        t1->Draw("Mh>>h_sb(25,100,150)",sbCut.data(),"e");
        t1->Draw(Form("hPt>>h_hptPass(35,%f,%f)",minhPt-10,maxhPt+10),sbCut.data(),"e");
        t1->Draw("MA0>>h_ma0Pass(60,200,800)",sbCut.data(),"e");
        
        t1->Draw("hPt","isTag&&(Mh>136||Mh<114)&MZp>400","e");
        c1->Print(pdfName.data());

        //t1->Draw("Mh>>h_sr(25,100,150)",sCut.data(),"e");
        //string sfCut = Form("isAntiTag&&(MA0<%f)&&(MA0>%f)&&(Mh<136)&&(Mh>114)&MZp>400",maxhPt,minhPt);
        //string sfCut = Form("isAntiTag&&(hPt<%f)&&(hPt>%f)&&MZp>400",maxhPt,minhPt);
        //string sfCut_nosig = Form("isAntiTag&&(hPt<%f)&&(hPt>%f)&&MZp>400&&(Mh>136||Mh<114)",maxhPt,minhPt);
        string sfCut_nosig = "99<1";
        string sfCut = "99<1";  // retrun false
        t1->Draw("Mh>>h_cr(25,100,150)",sfCut.data(),"e");
        t1->Draw("MA0>>h_ma0Fail(60,200,800)",sfCut_nosig.data(),"e");
        t1->Draw(Form("hPt>>h_hptFail(35,%f,%f)",minhPt-10,maxhPt+10),sfCut_nosig.data(),"e");
        //h_sr = (TH1F*)gDirectory->Get("h_sr");
        h_cr = (TH1F*)gDirectory->Get("h_cr");
        h_sb = (TH1F*)gDirectory->Get("h_sb");
        h_hptPass = (TH1F*)gDirectory->Get("h_hptPass");
        h_hptFail = (TH1F*)gDirectory->Get("h_hptFail");
        h_ma0Pass = (TH1F*)gDirectory->Get("h_ma0Pass");
        h_ma0Fail = (TH1F*)gDirectory->Get("h_ma0Fail");
        
        ///////////////////////////////////////////////////////
        cout << "start fill hist" << endl;
        for (int k=0;k<t1->GetEntries();k++) {
            if (k%100000==0) cout << k << "/" << t1->GetEntries() << endl;
            if (runShort&&k>5000) break;
            t1->GetEntry(k);
            
            
            //if (mh>150||mh<100) continue; // Mh boundary
            TLorentzVector vh = *vh1 +*vh2;
            if (isTag&&0) { 
            	if (hpt>maxhPt||hpt<minhPt) continue; //[200,300]
                TLorentzVector va0 = *va01 +*va02;
                if (m4b<400) continue;
                if (ma0<200) continue;
                if (mh>114&&mh<136) continue;
                //h_dphiPass->Fill(abs(vh.DeltaPhi(va0)));
            }
            if (isHfail) {
                if (nFailHCom==0) continue;
                bool find = false;
                static int iseed;
                iseed = gRandom->Integer(nFailHCom); 
                // mh
                v1 = (TLorentzVector*)arr->At(HfailInd1[iseed]);
                v2 = (TLorentzVector*)arr->At(HfailInd2[iseed]);
                float higgsranpt = (*v1+*v2).Pt();
                m4b_loop = (*v1+*v2+*vh1+*vh2).M();
                mh_loop = (*v1+*v2).M();
                if (higgsranpt>maxhPt||higgsranpt<minhPt) continue;
                if (m4b_loop<400) continue;
                if ((*vh1+*vh2).M()<200.) continue; // MA0
                
                if (mh_loop>114&&mh_loop<136) continue;
                if (mh_loop<100||mh_loop>150) continue;
                h_cr->Fill(mh_loop);
                h_dphiFail->Fill(abs(vh.DeltaPhi(*v1+*v2)));
                h_dRFail->Fill(vh.DeltaR(*v1+*v2));
                h_ma0Fail->Fill((*vh1+*vh2).M());
                h_hptFail->Fill(hpt);
            }
            //if (isTag&&m4b<400) continue;
        }
        cout << "end of fill hist" << endl;
        //////////////////////////////////////////////////////
        //h_sr->SetTitle("h_sr");
        //h_sr->SetLineColor(kGray+2);
        //h_sr->SetMarkerSize(0);
        //h_srRpf = (TH1F*)h_sr->Clone("h_srRpf");
        h_Rpf = (TH1F*)h_sb->Clone("h_Rpf");
        //h_srRpf->SetTitle("h_srRpf");
        h_Rpf->SetTitle(Form("h_Rpf_hPt%d-%d",(int)minhPt,(int)maxhPt));
        //h_srRpf->Sumw2();
        h_Rpf->Sumw2();
        //h_srRpf->Divide(h_cr);
        h_Rpf->Divide(h_cr);
        //h_srRpf->SetMarkerSize(0);
        h_Rpf->SetMarkerSize(0);


        
        
        //h_srRpf->SetLineColor(kGray+2);
        h_Rpf->SetMaximum(h_Rpf->GetMaximum()*1.6);
        h_Rpf->GetXaxis()->SetTitle("M_{h} (GeV)");
        h_Rpf->Fit(fpol,"IMF");
        //h_srRpf->Draw("esame");
        c1->Print(pdfName.data());
        //h_srRpf->Draw("e");
        //c1->Print(pdfName.data());
        //h_sr->Draw("e");
        //c1->Print(pdfName.data());
        gPad->Update();

        for (int j=0;j<fpol->GetNpar();j++) npar[j] = fpol->GetParameter(j);
        i = 0; 
        h_sb->SetMaximum(h_sb->GetMaximum()*1.4);
        h_sb->SetTitle(Form("h_sideband_hpt%d-%d",(int)minhPt,(int)maxhPt));
        h_sb->GetXaxis()->SetTitle("M_{h} (GeV)");
        h_sb->GetYaxis()->SetTitle("Number of Events");

        h_sb->Draw("e");
        //h_sr->Draw("esame");
        fpol->SetParameters(npar);
        h_cr->Multiply(fpol);
        h_cr->SetTitle("h_failXrpf");
        h_sb->SetLineColor(kBlack);
        h_cr->SetLineColor(kBlue);
        h_cr->SetMarkerColor(kBlue);
        h_cr->Draw("esame");
        c1->BuildLegend();
        c1->Print(pdfName.data());
        for (int j=0;j<fpol->GetNpar();j++) cout << npar[j] << "\t";
        cout << endl;
    //}
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
        int hPtp = (int)minhPt;
        int width = 100;
        //int hPtp = (int)h_MhMA0_rpf->GetYaxis()->GetBinLowEdge(i+1);
        string namep = Form("h_reduceMass_pass_hPt%dto%d_sig%d",hPtp,(int)maxhPt,i);
        string nameAl = Form("h_reduceMass_Alphabet_hPt%dto%d_sig%d",hPtp,(int)maxhPt,i);
        string name4bp = Form("h_4bMass_pass_hPt%dto%d_sig%d",hPtp,(int)maxhPt,i);
        string name4bAl = Form("h_4bMass_Alphabet_hPt%dto%d_sig%d",hPtp,(int)maxhPt,i);
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
    h_ma0Fail->Reset();
    h_hptFail->Reset();
    for (int i=0;i<t1->GetEntries();i++) {
        //.if (i>500000) break;
        if (runShort&&i>5000) break;
        if (i%100000==0) cout << i << "/" << t1->GetEntries() << endl;
        t1->GetEntry(i);
        TLorentzVector vh = *vh1 +*vh2;
        
        //if (isTag&&(mh>114&&mh<136)) continue;
        //double mtilde = m4bTilde(m4b,ma0,mh,ma0);
        if (isTag) {
            if (m4b<400) continue;
            if (hpt>maxhPt||hpt<minhPt) continue; //[200,300]
        	if (mh>150||mh<100) continue; // Mh boundary
        	if (mh>114&&mh<136) continue;
            double mtilde = m4bTilde(m4b,ma0,mh,300.);
            h_reduceMass_pass[0]->Fill(mtilde);
            h_4bMass_pass[0]->Fill(m4b);
            h_dphiPass->Fill(abs(vh.DeltaPhi(*va02+*va01)));
            h_dRPass->Fill(vh.DeltaR(*va01+*va02));
            h_a0ptPass->Fill((*va01+*va02).Pt());

        }
        if (isHfail){
        	
            if (nFailHCom==0) continue;
            static int iseed;
            iseed = gRandom->Integer(nFailHCom); 
            v1 = (TLorentzVector*)arr->At(HfailInd1[iseed]);
            v2 = (TLorentzVector*)arr->At(HfailInd2[iseed]);
            m4b_loop = (*v1+*v2+*vh1+*vh2).M();
            mh_loop = (*v1+*v2).M();
            float higgsranpt = (*v1+*v2).Pt();
            
            if (higgsranpt>maxhPt||higgsranpt<minhPt) continue;
            if (m4b_loop<400) continue;
            if (mh_loop<100||mh_loop>150) continue;
            if (mh_loop>114&&mh_loop<136) continue;
            if (vh.M()<200.) continue;
            
            
            double mtilde = m4bTilde(m4b_loop,vh.M(),mh_loop,300.);
            double weight = fpol->Eval(mh_loop);
            h_dphiWeight->Fill(abs(vh.DeltaPhi(*v1+*v2)),weight);
            h_dphiFail->Fill(abs(vh.DeltaPhi(*v1+*v2)));
            h_dRWeight->Fill(vh.DeltaR(*v1+*v2),weight);
            h_ma0Fail->Fill(vh.M(),weight);

            h_hptFail->Fill(higgsranpt,weight);
            h_reduceMass_alphabet[0]->Fill(mtilde,weight);
            h_4bMass_alphabet[0]->Fill(m4b_loop,weight);
            h_a0ptFail->Fill(vh.Pt());
            h_a0ptWeight->Fill(vh.Pt(),weight);
        }
    }
    TLegend *leg;
    h_hptFail->SetLineColor(kBlue);
    h_hptFail->SetMarkerColor(kBlue);
    h_ma0Fail->SetLineColor(kBlue);
    h_ma0Fail->SetMarkerColor(kBlue);

    h_hptPass->SetTitle("h_hptPass");
    h_hptFail->SetTitle("h_hptFailXrpf");
    h_hptPass->GetXaxis()->SetTitle("higgs P_{T} (GeV)");
    h_hptPass->GetYaxis()->SetTitle("Number of Events");
    h_hptPass->SetMaximum(h_hptPass->GetMaximum()*2);
    h_hptPass->Draw("e");
    h_hptFail->Draw("esame");
    leg = c1->BuildLegend(0.1,0.6,0.5,0.9);
    leg->SetHeader(Form("hPt%d-%d",(int)minhPt,(int)maxhPt),"C");
    leg->AddEntry((TObject*)0,Form("pass normalization = %.2f",h_hptPass->Integral()),"");
    leg->AddEntry((TObject*)0,Form("RpfXfail normalization = %.2f",h_hptFail->Integral()),"");
    c1->Print(pdfName.data());
    h_a0ptPass->SetLineColor(kBlue);
    h_a0ptFail->SetLineColor(kBlack);
    h_a0ptWeight->SetLineColor(kRed);
    h_a0ptPass->GetYaxis()->SetTitle("Number of Events");
    h_a0ptPass->GetXaxis()->SetTitle("A0 P_{T} (GeV)");
    h_a0ptPass->Draw("e");
    h_a0ptWeight->Draw("esame");
    leg = c1->BuildLegend(0.5,0.6,0.9,0.9);
    leg->SetHeader(Form("hPt%d-%d",(int)minhPt,(int)maxhPt),"C");
    leg->AddEntry((TObject*)0,Form("pass normalization = %.2f",h_a0ptPass->Integral()),"");
    leg->AddEntry((TObject*)0,Form("RpfXfail normalization = %.2f",h_a0ptWeight->Integral()),"");    
    c1->Print(pdfName.data());
    h_ma0Pass->SetTitle("h_ma0Pass");
    h_ma0Fail->SetTitle("h_ma0FailXrpf");
    h_ma0Pass->GetXaxis()->SetTitle("M_{A0} (GeV)");
    h_ma0Pass->GetYaxis()->SetTitle("Number of Events");
    h_ma0Pass->SetMaximum(h_ma0Pass->GetMaximum()*1.4);
    h_ma0Pass->Draw("e");
    h_ma0Fail->Draw("esame");
    leg = c1->BuildLegend(0.5,0.6,0.9,0.9);
    leg->SetHeader(Form("hPt%d-%d",(int)minhPt,(int)maxhPt),"C");
    leg->AddEntry((TObject*)0,Form("pass normalization = %.2f",h_ma0Pass->Integral()),"");
    leg->AddEntry((TObject*)0,Form("RpfXfail normalization = %.2f",h_ma0Fail->Integral()),"");
    c1->Print(pdfName.data());

    h_reduceMass_pass[0]->GetXaxis()->SetTitle("reduce mass (GeV)");
    h_reduceMass_pass[0]->GetYaxis()->SetTitle("Number of Events");
    //h_reduceMass_pass[0]->SetMaximum(h_reduceMass_pass[1]->GetMaximum()*1.2);
    h_reduceMass_pass[0]->GetXaxis()->SetLimits(400,1200.);
    //h_reduceMass_pass[0]->SetTitle(Form("h_reduceMass_hPt%d-%d",(int)minhPt,(int)maxhPt));
    h_reduceMass_pass[0]->Draw("e");
    h_reduceMass_alphabet[0]->Draw("esame");
    c1->Update();
    leg = c1->BuildLegend(0.5,0.6,0.9,0.9);
    leg->SetHeader(Form("hPt%d-%d",(int)minhPt,(int)maxhPt),"C");
    leg->AddEntry((TObject*)0,Form("pass normalization = %.2f",h_reduceMass_pass[0]->Integral()),"");
    leg->AddEntry((TObject*)0,Form("RpfXfail normalization = %.2f",h_reduceMass_alphabet[0]->Integral()),"");
    c1->Print(pdfName.data());
    h_4bMass_pass[0]->GetXaxis()->SetTitle("4b mass (GeV)");
    h_4bMass_pass[0]->GetYaxis()->SetTitle("Number of Events");
    //h_4bMass_pass[0]->SetTitle(Form("h_4bMass_hPt%d-%d",(int)minhPt,(int)maxhPt));
    //h_4bMass_pass[0]->SetMaximum(h_4bMass_pass[1]->GetMaximum()*1.2);
    h_4bMass_pass[0]->GetXaxis()->SetLimits(400.,1200.);
    h_4bMass_pass[0]->Draw("e");
    h_4bMass_alphabet[0]->Draw("esame");
    c1->Update();
    leg = c1->BuildLegend(0.5,0.6,0.9,0.9);
    leg->SetHeader(Form("hPt%d-%d",(int)minhPt,(int)maxhPt),"C");
    leg->AddEntry((TObject*)0,Form("pass normalization = %.2f",h_4bMass_pass[0]->Integral()),"");
    leg->AddEntry((TObject*)0,Form("RpfXfail normalization = %.2f",h_4bMass_alphabet[0]->Integral()),"");
    c1->Print(pdfName.data());
    h_dphiPass->SetLineColor(kBlue);
    h_dphiFail->SetLineColor(kBlack);
    h_dphiPass->GetYaxis()->SetTitle("Normalized to 1");
    h_dphiPass->GetXaxis()->SetTitle("#Delta#phi(h,A0)");
    h_dphiPass->DrawNormalized("e");
    h_dphiFail->DrawNormalized("esame");
    c1->Update();
    leg = c1->BuildLegend(0.1,0.6,0.5,0.9);
    leg->SetHeader(Form("hPt%d-%d",(int)minhPt,(int)maxhPt),"C");
    c1->Print(pdfName.data());
    
    h_dphiWeight->SetLineColor(kRed);
    h_dphiPass->GetXaxis()->SetTitle("#Delta#phi(h,A0)");
    h_dphiPass->GetYaxis()->SetTitle("Number of Events");
    h_dphiPass->Draw("e");
    h_dphiWeight->Draw("esame");
    c1->Update();
    leg = c1->BuildLegend(0.1,0.6,0.5,0.9);
    leg->SetHeader(Form("hPt%d-%d",(int)minhPt,(int)maxhPt),"C");
    leg->AddEntry((TObject*)0,Form("pass normalization = %.2f",h_dphiPass->Integral()),"");
    leg->AddEntry((TObject*)0,Form("RpfXfail normalization = %.2f",h_dphiWeight->Integral()),"");
    c1->Print(pdfName.data());
    h_dRPass->SetLineColor(kBlue);
    h_dRFail->SetLineColor(kBlack);
    h_dRPass->GetYaxis()->SetTitle("Normalized to 1");
    h_dRPass->GetXaxis()->SetTitle("#DeltaR(h,A0)");
    h_dRPass->DrawNormalized("e");
    h_dRFail->DrawNormalized("esame");
    c1->Update();
    leg = c1->BuildLegend(0.1,0.6,0.4,0.9);
    leg->SetHeader(Form("hPt%d-%d",(int)minhPt,(int)maxhPt),"C");
    c1->Print(pdfName.data());
    h_dRWeight->SetLineColor(kRed);
    h_dRPass->GetYaxis()->SetTitle("Number of Events");
    h_dRPass->Draw("e");
    h_dRWeight->Draw("esame");
    c1->Update();
    leg = c1->BuildLegend(0.1,0.6,0.4,0.9);
    leg->SetHeader(Form("hPt%d-%d",(int)minhPt,(int)maxhPt),"C");
    leg->AddEntry((TObject*)0,Form("pass normalization = %.2f",h_dRPass->Integral()),"");
    leg->AddEntry((TObject*)0,Form("RpfXfail normalization = %.2f",h_dRWeight->Integral()),"");
    c1->Print(pdfName.data());
    c1->Print((pdfName+"]").data());

}
