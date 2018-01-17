#include <iostream>
#include "TROOT.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TString.h"

int process_data(TString name){
  std::cout<<"processing: "<<name<<std::endl;

  TFile *f = new TFile(name,"read");
  TNtuple *t =(TNtuple *)f->Get("osc");
  TFile *fout = new TFile("out_" + name,"recreate");
  TString varList="q1:q2:q3:q4:ev:min1:min2:min3:min4:max1:max2:max3:max4:tm1:tm2:tm3:tm4:tmin1:tmin2:tmin3:tmin4:p1:p2:p3:p4";
  TTree *tout = new TTree("data","data processed. time (ns), C1... (V), q1... (nC)");

  Float_t time,tm1,tm2,tm3,tm4,tmin1,tmin2,tmin3,tmin4,c1,c2,c3,c4,p1,p2,p3,p4,evt,evt_prev,min1,min2,min3,min4,max1,max2,max3,max4,q1,q2,q3,q4,dt,hl,ll,phl,pll;
  Float_t  *data = new Float_t[varList.CountChar(':')+1];

  t->SetMaxEntryLoop(1e5);

  Int_t Nmeas=  t->Draw("C1","event==0","goff");

  t->SetMaxEntryLoop();

  Float_t *C1 = new Float_t[Nmeas];
  Float_t *C2 = new Float_t[Nmeas];
  Float_t *C3 = new Float_t[Nmeas];
  Float_t *C4 = new Float_t[Nmeas];
  Float_t *TIME = new Float_t[Nmeas];

  tout->Branch("C1",C1,Form("C1[%d]/F",Nmeas));
  tout->Branch("C2",C2,Form("C2[%d]/F",Nmeas));
  tout->Branch("C3",C3,Form("C3[%d]/F",Nmeas));
  tout->Branch("C4",C4,Form("C4[%d]/F",Nmeas));
  tout->Branch("time",TIME,Form("time[%d]/F",Nmeas));
  tout->Branch("measured",data,varList);

  ////////// GATE AND PEDESTAL LIMIT DEFINITION ////////////////
  ll = -10e-9;// low time signal edge
  hl = 200e-9;// high time signal edge
  pll = -80e-9;// low time pedestal edge 
  phl = -20e-9;// high time pedestal edge

  //////////////////////////////////////////////////////////////
  t->SetBranchAddress("time",&time);
  t->SetBranchAddress("C1",&c1);
  t->SetBranchAddress("C2",&c2);
  t->SetBranchAddress("C3",&c3);
  t->SetBranchAddress("C4",&c4);
  t->SetBranchAddress("event",&evt);

  Long_t Ne=t->GetEntries();
  t->GetEntry(0);
  evt_prev=evt;
  min1 = min2 = min3 = min4 = 100000;
  max1 = max2 = max3 = max4 = -100000;
  q1=q2=q3=q4=p1=p2=p3=p4=0;
  tm1=tm2=tm3=tm4=tmin1=tmin2=tmin3=tmin4=0;
  Int_t  count=0;
  dt = time;
  t->GetEntry(1);
  dt =time-dt;

  for (int i=0;i<Ne;i++)
  {
    t->GetEntry(i);
    if (ll<time&&time<hl)
    {
      q1+=-c1;q2+=-c2;q3+=-c3;q4+=-c4;
      tm1=-c1*time*1e9;tm2=-c2*time*1e9;tm3=-c3*time*1e9;tm4=-c4*time*1e9;
    }
    else if(pll<time&&time<phl)
    {
      p1+=c1;p2+=c2;p3+=c3;p4+=c4;
    }
    if(min1>c1){ min1 = c1;tmin1=time*1e9;}
    if(min2>c2){ min2 = c2;tmin2=time*1e9;}
    if(min3>c3){ min3 = c3;tmin3=time*1e9;}
    if(min4>c4){ min4 = c4;tmin4=time*1e9;}

    if(max1<c1) max1 = c1;
    if(max2<c2) max2 = c2;
    if(max3<c3) max3 = c3;
    if(max4<c4) max4 = c4;

    if (evt != evt_prev)
    {
      
      tm1=tm1/q1;
      tm2=tm2/q2;
      tm3=tm3/q3;
      tm4=tm4/q4;
      q1=q1*dt*1e9/50.;
      q2=q2*dt*1e9/50.;
      q3=q3*dt*1e9/50.;
      q4=q4*dt*1e9/50.;
      data[0]=q1;
      data[1]=q2;
      data[2]=q3;
      data[3]=q4;
      data[4]=evt;
      data[5]=min1;
      data[6]=min2;
      data[7]=min3;
      data[8]=min4;
      data[9]=max1;
      data[10]=max2;
      data[11]=max3;
      data[12]=max4;
      data[13]=tm1;
      data[14]=tm2;
      data[15]=tm3;
      data[16]=tm4;

      data[17]=tmin1;
      data[18]=tmin2;
      data[19]=tmin3;
      data[20]=tmin4;
      p1=p1*dt*1e9/50.*(hl-ll)/(phl-pll);
      p2=p2*dt*1e9/50.*(hl-ll)/(phl-pll);
      p3=p3*dt*1e9/50.*(hl-ll)/(phl-pll);
      p4=p4*dt*1e9/50.*(hl-ll)/(phl-pll);
      data[21]=p1;
      data[22]=p2;
      data[23]=p3;
      data[24]=p4;

      tout->Fill();
      evt_prev = evt;
      //reset variables
      min1 = min2 = min3 = min4 = 100000;
      max1 = max2 = max3 = max4 = -100000;
      q1=q2=q3=q4=p1=p2=p3=p4=0;
      tm1=tm2=tm3=tm4=tmin1=tmin2=tmin3=tmin4=0;
      ///////////////////
      count=0;
      C1[count]=c1; C2[count]=c2; C3[count]=c3; C4[count]=c4; TIME[count++]=time*1e9;
    }
    else
    {
      C1[count]=c1; C2[count]=c2; C3[count]=c3; C4[count]=c4; TIME[count++]=time*1e9;
    }
   
  }
  tout->Write("",TObject::kOverwrite);
  f->Close();
  fout->Close();

  return 0;
}
