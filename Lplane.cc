// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/DirectFinalState.hh"
#include "Rivet/Projections/HeavyHadrons.hh"

#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/LundGenerator.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class Lplane : public Analysis {
  public:
  //counter
  //int jetCounter = 0;
  int slice = 20;
  double kTcut = 1*GeV;
  double jetR = 0.4;
  double XMin = 0.0;
  double XMax = 4;
  double YMin = -5;
  double YMIN = 0.0;
  double ZMIN = 0;
  double ZMAX = 6;
  double YMax = 5;
  double sdzg = 0.5;
  double ZMin = log(1/0.5);
  double ZMax = 8.6*log(1/0.5);
  double Erad = 0*GeV;
  double MinJetPt = 55*GeV;
  double MaxJetPt = 100*GeV;
  int histoSlice = 20;

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(Lplane);


    /// @name Analysis methods
    /// @{
    void init() {
      const FinalState fs(Cuts::abseta < 4.5);
      FastJets jets(fs, JetAlg::ANTIKT, 0.4, JetMuons::NONE, JetInvisibles::NONE);
      declare(jets, "Jets");  
      ChargedFinalState tracks(Cuts::pT > 0.5*GeV && Cuts::abseta < 2.5);
      declare(tracks, "tracks");
      declare(HeavyHadrons(Cuts::pT > 5*GeV), "BHadrons");
     
      //histograms booking
      book(_h_2Dbjets, "bjets", histoSlice, XMin, XMax, histoSlice, Erad, MaxJetPt); 
      book(_h_2Dlightjets, "lightjets", histoSlice , XMin, XMax, histoSlice, Erad, MaxJetPt); 
      book(_h_hadroPt, "hadron pT", histoSlice, XMin, MaxJetPt); 
      book(_h_jetPt, "Jet pT", histoSlice, MinJetPt, MaxJetPt); 
      book(_h_sdjetPt, "sd Jet pT", histoSlice, MinJetPt, MaxJetPt); 
      book(_h_delta, "delta", histoSlice, XMin, XMax);
      book(_h_rg, "rg", histoSlice, XMin, XMax);
      book(_h_zg, "zg", histoSlice, ZMIN, sdzg);
      book(_h_2Dlund, "lund", histoSlice, XMin, XMax, histoSlice, ZMIN, ZMAX); //Testing Leticia's method.
      book(_njets, "njets");
 
      _h_vs.resize(slice);
      for (size_t i = 0; i < _h_vs.size(); ++i) {
        book(_h_vs[i], "vs" + std::to_string(i), slice, ZMin, ZMax);
      }
      _h_hs.resize(slice);
      for (size_t j = 0; j < _h_hs.size(); ++j) {
        book(_h_hs[j], "hs" + std::to_string(j), slice, XMin, XMax);
      }
    }

    //per event analysis
    void analyze(const Event& event) {

      //find bhadrons in events
      Particles bhadrons;

      //loop over particles in event record and check if they are final state bhadrons or not. Place fs bhadrons inside the bhadrons container.
      for(ConstGenParticlePtr p: HepMCUtils::particles(event.genEvent())) {
        if (!( PID::isHadron( p->pdg_id() ) && PID::hasBottom( p->pdg_id() )) ) continue;
        //std::cout<< "hadron id: " << p->pdg_id()<<std::endl;
        ConstGenVertexPtr dv = p->end_vertex();
        bool hasBdaughter = false;
        if ( PID::isHadron( p->pdg_id() ) && PID::hasBottom( p->pdg_id() )) { // b-hadron selection
          if (dv) {
            for(ConstGenParticlePtr pp: HepMCUtils::particles(dv, Relatives::CHILDREN)){
              if ( PID::isHadron( pp->pdg_id() ) && PID::hasBottom( pp->pdg_id()) ) {
                hasBdaughter = true;
              }
            }
          }
          if (hasBdaughter) continue;
          bhadrons += Particle(*p);    
        }
      }
      //std::cout<< "bhadrons.size(): "<< bhadrons.size()<<std::endl;

      // Retrieve clustered jets, sorted by pT, with a minimum pT & eta cut
      const Jets jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 30*GeV && Cuts::abseta < 2.1);
       //separate heavy and light jets 
      Jets bjets, lightjets;

      //loop over every jet in jet container
      for (const Jet& jet : jets){
        if (jet.pT() < MinJetPt) continue;
        if (jet.pT() > MaxJetPt) continue;  //upper limit because DC's dependance on energy
        bool Bjet = false;
        //tag
        for (const Particle& p: bhadrons) {
          double hadroPt = p.pT();
          if (hadroPt > 100*GeV) continue;
          if (deltaR(jet,p, PSEUDORAPIDITY) < jetR) {
            //std::cout<<"jet hadron pid: "<< abs(p.pid())<<std::endl;
            //std::cout<<"hadroPt: "<< hadroPt<< std::endl;
            _h_hadroPt->fill(hadroPt / GeV);
            bjets.push_back(jet);
            Bjet = true;
          }
        }
        if (!Bjet) {
          lightjets.push_back(jet);
        }
      }
      //std::cout<<"bjets.size(): "<< bjets.size()<<std::endl;    
      //std::cout<<"lightjets.size(): "<< lightjets.size()<<std::endl;
      
      //make sure you have a bjet
      if (bjets.size()< 1) vetoEvent;

      //select leading jet
      const Jet& j1 = bjets[0];
      _njets->fill(1); //count your 1 jet
      double jetPt = j1.pT();
      //std::cout<<"jetPt: "<< jetPt<<std::endl; 
      _h_jetPt->fill(jetPt/ GeV);   

      //define sd parameters
      double z_cut = 0.10;
      double beta  = 1.0;
      fjcontrib::SoftDrop sd(beta, z_cut);

      //groom
      PseudoJet sd_j1 = sd(j1);
      double sdpT = sd_j1.pt();

      //because soft drop is a groomer (not a tagger), it should always return a soft-dropped jet
      assert(sd_j1 != 0); 
      _h_sdjetPt->fill(sdpT/ GeV);
  
      //tracks 
      /*const Particles& tracks = apply<ChargedFinalState>(event, "tracks").particlesByPt();
      Particles intracks;  //declare tracks container     

      //analyse jet

      //find tracks in jet
      for (const Particle& p: tracks){
        const double dr = deltaR(sd_j1, p, PSEUDORAPIDITY);
        if (dr > jetR) continue;
        intracks.push_back(p);
        //std::cout<<" Track PID: "<< p.pid()<<std::endl;
      }
      //std::cout<<"intracks.size(): "<<intracks.size()<<std::endl;
      //Particles constituents = sd_j1.constituents;  //declare tracks container  


      //re-cluster tracks with CA algorithm
      JetDefinition tjet_def(fastjet::cambridge_algorithm, 10);
      ClusterSequence tjet_cs(intracks, tjet_def);
      vector<PseudoJet> tjets = fastjet::sorted_by_pt (tjet_cs.inclusive_jets(0.0));

      if(tjets.size() < 1) vetoEvent;  //CA must return at least 1 jet
      //std::cout<<"tjets.size(): " <<tjets.size()<< std::endl;*/
      double rg = sd_j1.structure_of<fjcontrib::SoftDrop>().delta_R();
      double zg = sd_j1.structure_of<fjcontrib::SoftDrop>().symmetry();
      double rg_p = -log(rg);

      _h_rg->fill(rg_p);
      _h_zg->fill(zg);

      //decluster jet with the lund generator
      fjcontrib::LundGenerator lund; //declare lund generator
      vector<fjcontrib::LundDeclustering> declusts = lund(sd_j1);  //decluster first pseudojet.
      //std::cout<<"declusts.size(): "<< declusts.size()<<std::endl;
      for (size_t idecl = 0; idecl < declusts.size(); ++idecl) {  //continue declustering jet till you reach core
        pair<double,double> coords = declusts[idecl].lund_coordinates(); //find lund coordinates for each declustering step
        double X = coords.first; //ln(1/theta))
        double Y = coords.second; //this is actually ln(kt)
        double Z = - log(declusts[idecl].z()); 
        //double kT = declusts[idecl].kt();
        double E = exp(X + Y + Z);  //radiator energy
        
        if (X > XMin && X < XMax && E > Erad && E < MaxJetPt) { 
          //std::cout<<"E: "<<E<<std::endl;
          //std::cout<<"kT: "<<kT<<std::endl;
          _h_2Dbjets->fill(X,E); //fill bjets lund plane
          _h_delta->fill(X); //fill bjets lund plane
          //_h_2Dlightjets->fill(X,E); //fill lightjets lund plane

          double hdiv = (double)XMax/(double)slice;
          int i = floor(X/hdiv);
         // std::cout<<"i: "<<i<<std::endl;
          _h_vs[i]->fill(E);

          double vdiv = (double)(MaxJetPt - Erad)/(double)slice;
          int j = floor((E - Erad)/vdiv);
          //std::cout<<"j: "<<j<<std::endl;
          _h_hs[j]->fill(X);
        }  //end if statement
        if (X > XMin && X < XMax && Z > ZMin && Z < ZMax) {
          _h_2Dlund->fill(X,Z); //fill angular sep histo
        }
      } //end of declust for loop
    }  //end analyze() function

    void finalize() {
      //std::cout<<"_njets: "<<_njets<<std::endl;
      //const double jetCounterW = _njets->sumW();
      //YODA::Counter jetCounter = *_njets;
      //std::cout<<"jetCounterW: "<<jetCounterW<<std::endl;
      //std::cout<<"sumW(): "<<sumW()<<std::endl;

      //normalize the histograms using x section
      //const double scaling = crossSection()/picobarn/sumW();
      //const double scaling = 1/jetCounterW;
      const double scaling = 1/sumW();
      scale(_h_2Dbjets, scaling);
      scale(_h_2Dlightjets, scaling);
      scale(_h_2Dlund, scaling);
      scale(_h_vs, scaling);
      scale(_h_hs, scaling);
      scale(_h_jetPt, scaling);
      scale(_h_sdjetPt, scaling);
      scale(_h_hadroPt, scaling);
      scale(_njets, scaling);
      scale(_h_delta, scaling);
      scale(_h_rg, scaling);
      scale(_h_zg, scaling);
    }
  	private:
      Histo2DPtr _h_2Dbjets, _h_2Dlightjets, _h_2Dlund;
      vector<Histo1DPtr> _h_vs, _h_hs;
      Histo1DPtr _h_jetPt, _h_hadroPt, _h_delta, _h_sdjetPt, _h_rg, _h_zg;
      CounterPtr _njets;

  };

  RIVET_DECLARE_PLUGIN(Lplane);

}
