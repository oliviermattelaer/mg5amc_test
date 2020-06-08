#include <iostream>
#include "Pythia8/Pythia.h"
#include "LHAFortran_aMCatNLO.h"
#include "fstream"

using namespace std;
using namespace Pythia8;

extern "C" { 

  // set up a global instance of pytia8
  Pythia pythia;
  // set up a global instance of LHAup
  //MyLHAupFortran lhareader;
  MyLHAupFortran lhareader(&pythia.settings);
  LHA3FromPythia8 lhawriter(&pythia.event, &pythia.settings, &pythia.info,
    &pythia.particleData);

  PrintFirstEmission printFirstEmission(&lhawriter); 

  // a counter for the number of event
  int iEvent = 0;

  // an initialisation function
  void pythia_init_(char input[500]) {
    string cmdFilePath(input);
    // Remove whitespaces
    cout << "--" << input << "--" << endl;
    cout << "--" << cmdFilePath << "--" << endl;
    while(cmdFilePath.find(" ", 0) != string::npos)
      cmdFilePath.erase(cmdFilePath.begin()+cmdFilePath.find(" ",0));
    if (cmdFilePath!="" && !(fopen(cmdFilePath.c_str(), "r"))) {
      cout<<"Pythia8 input file ' "<<cmdFilePath<<" ' not found."<<endl;
      abort();
    }
    lhareader.setInit();
    // Example of a user hook for storing in the out stream the event after the first emission.
    pythia.setUserHooksPtr(&printFirstEmission);
    if (cmdFilePath!="") {
      cout<<"Initialising Pythia8 from cmd file '"<<cmdFilePath<<"'"<<endl;		
      pythia.readFile(cmdFilePath.c_str());
    } else {
     cout<<"Using default initialization of Pythia8."<<endl;
     pythia.readString("Beams:frameType=5");
     pythia.readString("Check:epTolErr=1.0000000000e-02");
    }
    pythia.setLHAupPtr(& lhareader);
    pythia.init();
    pythia.mergingPtr->setLHAPtr(&lhawriter);
    // Flag that Pythia8 intiialisation has been performed.
    pythia_control_.is_pythia_active = 1;
  }

  // an initialisation function
  void pythia_init_default_(int& idIn1, int& idIn2, int outIDs [10], double masses[26]) {
    lhareader.setInit();
    // Example of a user hook for storing in the out stream the event after the first emission.
    pythia.setUserHooksPtr(&printFirstEmission);

    // Reconstruct the process string.
    string processString = "";
    // Set incoming particles.
    if (idIn1 == 2212) processString += "p";
    if (idIn1 == 11)   processString += "e-";
    if (idIn1 ==-11)   processString += "e+";
    if (idIn2 == 2212) processString += "p";
    if (idIn2 == 11)   processString += "e-";
    if (idIn2 ==-11)   processString += "e+";
    processString += ">";
    // Set outgoing particles.
    bool foundOutgoing = false;
    for (int i=0; i < 10; ++i) {
      if (outIDs[i]==0) continue;
      if (outIDs[i]==2212) {
        processString += "j";
      } else {
        ostringstream proc;
        proc << "{" << pythia.particleData.name(outIDs[i]) << "," << outIDs[i] << "}";
        processString += proc.str();
      }
    }

    // Initialize masses.
    for (int i=1; i <= 25; ++i){
      if (masses[i]<0.) continue;
      stringstream s;
      // Need a non-zero muon mass to get correct Higgs width. Otherwise gets a NAN. Need to be fixed later.
      if (i==13) continue;
      s << i << ":m0 =" << masses[i];
      pythia.readString(s.str());
    }

    cout<<"Using default initialization of Pythia8."<<endl;
    pythia.readString("Beams:frameType=5");
    pythia.readString("Check:epTolErr=1.0000000000e-02");
    // Disallow Pythia to overwrite parts of Les Houches input.
    pythia.readString("LesHouches:setQuarkMass       = 0");
    pythia.readString("LesHouches:setLeptonMass      = 0");
    pythia.readString("LesHouches:mRecalculate       = -1.0");
    pythia.readString("LesHouches:matchInOut         = off");
    // Switch off most of Pythia
    pythia.readString("ProcessLevel:resonanceDecays  = off");
    pythia.readString("BeamRemnants:primordialKT     = off");
    pythia.readString("TimeShower:QEDshowerByQ       = off");
    pythia.readString("TimeShower:QEDshowerByL       = off");
    pythia.readString("TimeShower:QEDshowerByGamma   = off");
    pythia.readString("TimeShower:QEDshowerByOther   = off");
    pythia.readString("SpaceShower:QEDshowerByQ      = off");
    pythia.readString("SpaceShower:QEDshowerByL      = off");
    pythia.readString("PartonLevel:MPI               = off");
    pythia.readString("HadronLevel:all               = off");
    pythia.readString("PartonLevel:Remnants          = off");
    pythia.readString("Check:Event                   = off");
    pythia.readString("PartonLevel:FSRinResonances   = off");
    pythia.readString("PDF:lepton                    = off");
    pythia.readString("Print:quiet                   = on");
    pythia.readString("Beams:setProductionScalesFromLHEF = off");
    pythia.readString("Check:abortIfVeto               = on");
    // Merging settings to be able to use histories only.
    pythia.readString("Merging:nRequested              = 0");
    pythia.readString("Merging:mayRemoveDecayProducts  = on");
    pythia.readString("merging:doptlundmerging         = on");
    pythia.settings.word("Merging:Process", processString);
    pythia.readString("merging:tms                     = 1000000");
    pythia.readString("merging:includeWeightInXSection = off");
    pythia.readString("merging:njetmax                 = 1000");
    pythia.readString("merging:applyveto               = off");


    pythia.setLHAupPtr(& lhareader);
    pythia.init();
    pythia.mergingPtr->setLHAPtr(&lhawriter);
    // Flag that Pythia8 intiialisation has been performed.
    pythia_control_.is_pythia_active = 1;
  }

  // a function to shower and analyse events
  void pythia_setevent_() {
    if (!lhareader.is_initialised()) {
      lhareader.setInit();
      pythia.init();
    }
    //This should set the LHA event using fortran common blocks
    lhareader.setEvent();
  }

  // a function to shower and analyse events
  void pythia_next_() {
    if (!lhareader.is_initialised()) {
      lhareader.setInit();
      pythia.init();
    }
    pythia.next();
    ++iEvent;
  }

  void pythia_get_stopping_info_( double scales [100][100],
    double mass [100][100] ) {
    pythia.mergingPtr->getStoppingInfo(scales, mass);
  }

  void pythia_get_dead_zones_( bool dzone [100][100] ) {
    pythia.mergingPtr->getDeadzones( dzone);
  }

  void pythia_clear_() { pythia.mergingPtr->clear(); }

}
