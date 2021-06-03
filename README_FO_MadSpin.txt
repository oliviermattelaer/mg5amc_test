
This MG5 version contains the modifications in order for someone to run ttx@NLO followed by LO decays at Fixed Order. The output is an LHE file containing the decayed phase-space points in an event-like format.

- Syntax for LO QCD production:
generate p p > t t~ QCD=2 QED=0 [LOonly=QCD]

- Syntax for NLO QCD production:
generate p p > t t~ QCD=2 QED=0 [QCD]

- In the generated process folder one should choose in Cards/FO_analyse_card.dat
FO_ANALYSIS_FORMAT = LHE

- In the generated process folder one should choose in Cards/madspin_card.dat
# MadSpin FO
set spinmode onshell
set fixed_order True
# specify the decay for the final state particles
decay t > w+ b, w+ > e+ ve
decay t~ > w- b~, w- > mu- vm~


Attention to these points:

- After the production lhe generated file MadSpin generates locally a lot of events and consumes memory, therefore the requested accuracy in the run_card.dat should not be lower than 0.001
- The output lhe file is weighted and one should take it into account while reading the file