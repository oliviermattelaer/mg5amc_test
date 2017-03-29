#include <iostream>
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/LHAFortran.h"
#include "Pythia8Plugins/HepMC2.h"
#include "Pythia8Plugins/aMCatNLOHooks.h"
#include "Pythia8Plugins/CombineMatchingInput.h"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"
#include "HepMC/IO_BaseClass.h"
#include "HepMC/IO_HEPEVT.h"
#include "HepMC/HEPEVT_Wrapper.h"
#include "fstream"

using namespace std;
using namespace Pythia8;

// let us add a new class that inherits fomr LHAupFortran
class MyLHAupFortran : public LHAupFortran {
  public:

  MyLHAupFortran(){
    initialised = false;
  }

  //the common blocks should be alredy filled at the fortran level
  //so simply return true
  bool fillHepRup(){
    initialised = true;
    return true;
  }
  bool fillHepEup(){
    return true;
  }
  
  bool is_initialised(){
    return initialised;
  }

  private:
  bool initialised;
};


extern "C" {
  extern struct {
    double EVWGT;
  } cevwgt_;
}
#define cevwgt cevwgt_

extern "C" { 
  void pyabeg_(int&,char(*)[50]);
  void pyaend_(double&);
  void pyanal_(int&,double(*));

  // set up a global instance of pytia8
  Pythia GLOB_pythia;
  // set up a global instance of LHAup
  MyLHAupFortran GLOB_LHAup;
  // a counter for the number of event
  int iEvent = 0;

  // an initialisation function
  void init_pythia_(char input[50]) {

    // look for possible position of the input file
    char py_input[80];
    if (!(fopen(input, "r"))) {
      sprintf(py_input, "../../../MCatNLO/%s", input);
    }
    else {
      sprintf(py_input, "%s", input);
    }

    cout<<"Initialising PYTHIA8 with "<<py_input<<endl;
    GLOB_pythia.readFile(py_input);

    // initialise the analysis
    int cwgtinfo_nn =1;
    char cwgtinfo_weights_info[1024][50];
    sprintf(cwgtinfo_weights_info[0], "%50s", "central value");
    pyabeg_(cwgtinfo_nn,cwgtinfo_weights_info);

    GLOB_pythia.setLHAupPtr(& GLOB_LHAup);
    // strategy is heprup.idwtup

  }

  // a function to shower and analyse events
  void pythia_shower_and_analyse_() {
    int cwgtinfo_nn;
    double cwgt_ww[1024];

    if (!GLOB_LHAup.is_initialised()) {
      GLOB_LHAup.setInit();
      GLOB_pythia.init();
    }
    //This should set the LHA event using fortran common blocks
    GLOB_LHAup.setEvent();

    // pythia does shower and hadronisation here
    GLOB_pythia.next();

    iEvent++;

    // Create the HepMC event object
    HepMC::IO_BaseClass *_hepevtio;
    HepMC::Pythia8ToHepMC ToHepMC;
    HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
    ToHepMC.fill_next_event(GLOB_pythia, hepmcevt);

    hepmcevt->set_event_number(iEvent);

    //define the IO_HEPEVT
    _hepevtio = new HepMC::IO_HEPEVT;
    _hepevtio->write_event(hepmcevt);
    
    //event weight
    cevwgt.EVWGT=hepmcevt->weights()[0];

    cwgt_ww[0]=cevwgt.EVWGT;

    // for the moment, force just one weight
    cwgtinfo_nn = 1;
    pyanal_(cwgtinfo_nn,cwgt_ww);
    delete hepmcevt;
  }

  //a function to close everything
  void end_pythia_(double& norm) {
    pyaend_(norm);
    GLOB_pythia.stat();
  }

}

