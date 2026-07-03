#include <iostream>
#include "TFile.h"
#include "TTree.h"

void read_interactions(const char* filename = "input.root") {

  // Open the ROOT file
  TFile *file = TFile::Open(filename, "READ");
  if (!file || file->IsZombie()) {
    std::cout << "Error opening file!" << std::endl;
    return;
  }

  // Get the TTree
  TTree *tree = (TTree*)file->Get("Tracks");
  if (!tree) {
    std::cout << "Error: TTree 'Tracks' not found!" << std::endl;
    return;
  }

  // ---------------------------
  // Variables to read the tree
  // ---------------------------

  Int_t event;              // Event number
  Int_t n_interaction;      // Number of interactions

  // Max number of interactions: set as 100
  const int MAX = 100;

  Int_t TrackID[MAX];       // Track ID number
  Int_t ParentID[MAX];      // Parent ID number
  Int_t ProcID[MAX];        // Process ID, related to EProcessType within the file include/WCSimEnumerations.hh
  Int_t PID[MAX];           // Particle ID number

  Float_t Time[MAX];        // Interaction time
  Float_t pos[MAX][3];      // Interaction position

  Float_t preP[MAX];        // pre-step momentum
  Float_t preDir[MAX][3];   // pre-step direction

  Float_t postP[MAX];       // post-step momentum
  Float_t postDir[MAX][3];  // post-step direction

  Float_t Ploss[MAX];       // Momentum loss at interaction (preP - postP)

  Int_t n_pip[MAX];         // Number of secondary pi+ at interaction
  Int_t n_pim[MAX];         // Number of secondary pi- at interaction
  Int_t n_muons[MAX];       // Number of secondary mu at interaction
  Int_t n_pi0[MAX];         // Number of secondary pi0 at interaction
  Int_t n_other[MAX];       // Number of secondary charged particles above Cherenkov threshold at interaction

  Int_t isBoundary[MAX];    // Is this position cross a boundary? 0 or 1

  // ---------------------------
  // Set branch addresses
  // ---------------------------

  tree->SetBranchAddress("event", &event);
  tree->SetBranchAddress("n_interaction", &n_interaction);

  tree->SetBranchAddress("TrackID", TrackID);
  tree->SetBranchAddress("ParentID", ParentID);
  tree->SetBranchAddress("ProcID", ProcID);
  tree->SetBranchAddress("PID", PID);

  tree->SetBranchAddress("Time", Time);
  tree->SetBranchAddress("pos", pos);

  tree->SetBranchAddress("preP", preP);
  tree->SetBranchAddress("preDir", preDir);

  tree->SetBranchAddress("postP", postP);
  tree->SetBranchAddress("postDir", postDir);

  tree->SetBranchAddress("Ploss", Ploss);

  tree->SetBranchAddress("n_pip", n_pip);
  tree->SetBranchAddress("n_pim", n_pim);
  tree->SetBranchAddress("n_muons", n_muons);
  tree->SetBranchAddress("n_pi0", n_pi0);
  tree->SetBranchAddress("n_other", n_other);

  tree->SetBranchAddress("isBoundary", isBoundary);

  // ---------------------------
  // Loop over events
  // ---------------------------

  Long64_t nentries = tree->GetEntries();
  std::cout << "Total entries in tree: " << nentries << std::endl;

  for (Long64_t i = 0; i < nentries; i++) {

    if (i >= 10) break;  // Only first 10 events

    tree->GetEntry(i);

    std::cout << "\n==============================" << std::endl;
    std::cout << "Event: " << event << std::endl;
    std::cout << "Number of interactions: " << n_interaction << std::endl;

    for (int j = 0; j < n_interaction; j++) {

      std::cout << "\n  --- <Interaction " << j << " ---" << std::endl;

      // --------------------------------------------------
      // FIRST STEP → Initial particle
      // --------------------------------------------------
      if (j == 0) {
        std::cout << "  [Initial Particle]" << std::endl;

        std::cout << "  PID: " << PID[j] << std::endl;
        std::cout << "  Time: " << Time[j] << " ns" << std::endl;

        std::cout << "  Position: ("
                  << pos[j][0] << ", "
                  << pos[j][1] << ", "
                  << pos[j][2] << ") cm" << std::endl;

        std::cout << "  Momentum: " << preP[j] << " MeV/c" << std::endl;

        std::cout << "  Direction: ("
                  << preDir[j][0] << ", "
                  << preDir[j][1] << ", "
                  << preDir[j][2] << ")" << std::endl;

        continue;
      }

      // --------------------------------------------------
      // BOUNDARY CROSSING
      // --------------------------------------------------
      if (isBoundary[j] == 1) {
        std::cout << "  [Boundary Crossing]" << std::endl;

        std::cout << "  PID: " << PID[j] << std::endl;
        std::cout << "  Time: " << Time[j] << " ns" << std::endl;

        std::cout << "  Position: ("
                  << pos[j][0] << ", "
                  << pos[j][1] << ", "
                  << pos[j][2] << ") cm" << std::endl;

        std::cout << "  Momentum: " << preP[j] << " MeV/c" << std::endl;

        std::cout << "  Direction: ("
                  << preDir[j][0] << ", "
                  << preDir[j][1] << ", "
                  << preDir[j][2] << ")" << std::endl;
      }

      // --------------------------------------------------
      // PHYSICS INTERACTION
      // --------------------------------------------------
      else {
        std::cout << "  [Interaction]" << std::endl;

        /*
        Process ID list can be found in include/WCSimEnumerations.hh
        The function WCSimEnumerations::EnumAsString can be used from src/WCSimEnumerations.cc to have the ID -> name
        Only exception is ID = -1, which is the process "Transportation" (no interaction)
        The process "Scintillation" is used to describe decay at rest while process "Decay" is decay in flight
        Both "Scintillation" and "Decay" have the muon as post-step information (updated)
        The process "pi+Inelastic" can describe absorption, charge exchange, inelastic scattering or even resonnance
        */

        std::cout << "  Process ID: " << ProcID[j] << std::endl;
        std::cout << "  PID: " << PID[j] << std::endl;
        std::cout << "  Time: " << Time[j] << " ns" << std::endl;

        std::cout << "  Position: ("
                  << pos[j][0] << ", "
                  << pos[j][1] << ", "
                  << pos[j][2] << ") cm" << std::endl;

        // Momentum evolution
        std::cout << "  Momentum: " << preP[j]
                  << " -> " << postP[j]
                  << " MeV/c"
                  << " (loss: " << Ploss[j] << " MeV/c)" << std::endl;

        // Direction evolution
        std::cout << "  Direction: ("
                  << preDir[j][0] << ", "
                  << preDir[j][1] << ", "
                  << preDir[j][2] << ") -> ("
                  << postDir[j][0] << ", "
                  << postDir[j][1] << ", "
                  << postDir[j][2] << ")" << std::endl;

        // Secondaries
        std::cout << "  Secondaries: "
                  << "pi+ = " << n_pip[j]
                  << ", pi- = " << n_pim[j]
                  << ", mu = " << n_muons[j]
                  << ", pi0 = " << n_pi0[j]
                  << ", other = " << n_other[j] // Number of secondary charged particles above Cherenkov threshold at interaction
                  << std::endl;
      }
    }
  }

  // Close file
  file->Close();
}