#include <vector>
#include <iostream>
#include <string>

using namespace std;

//#define test___
#ifdef test___
int total_tracks =0;
int total_events =0;
int total_tracks2=0;
int total_events2=0;
#endif

class Track{
  public:
  Track(float pt1, float eta1, float phi01, int charge1, float eff1,int ptbin11, int ptbin21);

  Track(){}

  float get_pt()    {return pt;}
  float get_eta()   {return eta;}
  float get_phi0()  {return phi0;}
  int   get_charge(){return charge;}
  float get_eff()   {return eff;}
  float get_ptbin1(){return ptbin1;}
  float get_ptbin2(){return ptbin2;}

  private:
  float pt,eta,phi0,eff;
  int charge,ptbin1,ptbin2;
};


Track::Track(float pt1, float eta1, float phi01, int charge1, float eff1,int ptbin11, int ptbin21){
  pt=pt1; eta=eta1; phi0=phi01;charge=charge1;eff=eff1;ptbin1=ptbin11;ptbin2=ptbin21;
}


class Event{
  private:
  int id,cent;
  float zvtx;
  vector<Track*> Tracks;

  public:
  Event(int id1, int cent1, float zvtx1);
  Event(){}
 ~Event();
  void AddTrack( float pt1, float eta1, float phi01, int charge1, float eff1,int ptbin11, int ptbin21);
  int get_id();
  int get_npart();
  Track* GetTrack(int i);
};

int Event::get_id(){return id;}
int Event::get_npart(){return Tracks.size();}

Track* Event::GetTrack(int i){return Tracks.at(i);}

void Event::AddTrack( float pt1, float eta1, float phi01, int charge1, float eff1,int ptbin11, int ptbin21){
  Track *new_track=new Track(pt1,eta1,phi01,charge1,eff1,ptbin11,ptbin21);
  Tracks.push_back(new_track);
  #ifdef test___
  total_tracks++;
  #endif
}

Event::Event::Event(int id1, int cent1, float zvtx1){
  id=id1; cent=cent1; zvtx=zvtx1;
  #ifdef test___
  total_events++;
  #endif
}



Event::~Event(){
  int ntrk=Tracks.size();
  for(int itrk=0;itrk<ntrk;itrk++){
    delete Tracks.at(itrk);
    #ifdef test___
    total_tracks--;
    #endif
  }
  #ifdef test___
  total_events--;
  #endif
  Tracks.clear();  
}


typedef Event* EVENT_PTR;


