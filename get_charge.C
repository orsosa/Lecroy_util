#include "TROOT.h"
#include "TNtuple.h"
#include "TH1F.h"
#include "TFile.h"

int get_charge(char* name="Todo_Soldado.root"){
  cout<<"updating: "<<name<<endl;

  TFile *f = new TFile(name,"update");
  TNtuple * t =(TNtuple *)f->Get("osc");
  char varList[100]="q:ev:min";
  TNtuple *tq = new TNtuple("tq","charge tuple q (nC), min (V)",varList);
  Float_t time,c1,evt,evt_prev,min,q,dt,hl,ll;
  // gate definition.
  ll = -10e-9;// low time edge
  hl = 300e-9;// high time edge
  t->SetBranchAddress("time",&time);
  t->SetBranchAddress("C1",&c1);
  t->SetBranchAddress("event",&evt);
  Long_t Ne=t->GetEntries();
  t->GetEntry(0);
  evt_prev=evt;
  min=100;
  q=0;
  dt = time;
  t->GetEntry(1);
  dt =time-dt;
  for (int i=0;i<Ne;i++)
  {
    t->GetEntry(i);
    if (ll<time&&time<hl) q+=-c1;
    if(min>c1) min = c1;
    if (evt != evt_prev)
    {
      q=q*dt*1e9/50;
      tq->Fill(q,evt_prev,min);
      evt_prev = evt;
      q=0;
      min=100;
    }
    
  }
  tq->Write("",TObject::kOverwrite);
  f->Close();
  return 0;
}
