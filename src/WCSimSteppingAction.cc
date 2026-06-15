#include "WCSimSteppingAction.hh"

#include <stdlib.h>
#include <stdio.h>
#include <WCSimRootEvent.hh>
#include <G4SIunits.hh>
#include <G4OpticalPhoton.hh>

#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4VParticleChange.hh"
#include "G4SteppingVerbose.hh"
#include "G4SteppingManager.hh"
#include "G4PVParameterised.hh"
#include "G4PVReplica.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4OpBoundaryProcess.hh"

G4int WCSimSteppingAction::n_photons_through_mPMTLV = 0;
G4int WCSimSteppingAction::n_photons_through_acrylic = 0;
G4int WCSimSteppingAction::n_photons_through_gel = 0;
G4int WCSimSteppingAction::n_photons_on_blacksheet = 0;
G4int WCSimSteppingAction::n_photons_on_smallPMT = 0;

///////////////////////////////////////////////
///// BEGINNING OF WCSIM STEPPING ACTION //////
///////////////////////////////////////////////


WCSimSteppingAction::WCSimSteppingAction(WCSimRunAction *myRun, WCSimDetectorConstruction *myDet) : runAction(myRun), det(myDet) {

}

void WCSimSteppingAction::UserSteppingAction(const G4Step* aStep)
{
    const G4Event *event = G4EventManager::GetEventManager()->GetConstCurrentEvent();
    if(event->IsAborted() || event->GetEventID() < 0)
      return;
  //DISTORTION must be used ONLY if INNERTUBE or INNERTUBEBIG has been defined in BidoneDetectorConstruction.cc
  
  const G4Track* track       = aStep->GetTrack();
  
  // Not used:
  //const G4Event* evt = G4RunManager::GetRunManager()->GetCurrentEvent();
  //G4VPhysicalVolume* volume  = track->GetVolume();
  //G4SDManager* SDman   = G4SDManager::GetSDMpointer();
  //G4HCofThisEvent* HCE = evt->GetHCofThisEvent();

  
  // Debug for photon tracking
  G4StepPoint* thePrePoint = aStep->GetPreStepPoint();
  G4VPhysicalVolume* thePrePV = thePrePoint->GetPhysicalVolume();

  G4StepPoint* thePostPoint = aStep->GetPostStepPoint();
  G4VPhysicalVolume* thePostPV = thePostPoint->GetPhysicalVolume();

  //G4OpBoundaryProcessStatus boundaryStatus=Undefined;
  //static G4ThreadLocal G4OpBoundaryProcess* boundary=NULL;  //doesn't work and needs #include tls.hh from Geant4.9.6 and beyond
  G4OpBoundaryProcess* boundary=NULL;
  
  //find the boundary process only once
  if(!boundary){
    G4ProcessManager* pm
      = aStep->GetTrack()->GetDefinition()->GetProcessManager();
    G4int nprocesses = pm->GetProcessListLength();
    G4ProcessVector* pv = pm->GetProcessList();
    G4int i;
    for( i=0;i<nprocesses;i++){
      if((*pv)[i]->GetProcessName()=="OpBoundary"){
	boundary = (G4OpBoundaryProcess*)(*pv)[i];
	break;
      }
    }
  }

  // Scattering output to fill ---Loris
  bool is_pic = std::abs(track->GetDefinition()->GetPDGEncoding())==211;//checking if pi+/pi- or not

  if(is_pic || std::abs(track->GetDefinition()->GetPDGEncoding())==13){// also accept tracking muon
    auto sctRunAction = const_cast<WCSimRunAction*>(static_cast<const WCSimRunAction*>(G4RunManager::GetRunManager()->GetUserRunAction()));
    bool isNewEvent = sctRunAction && event->GetEventID() != sctRunAction->GetPreviousEvent();

    bool isNewVolume = thePrePoint->GetStepStatus() == fGeomBoundary;
    if (isNewVolume) {
      auto preVol  = thePrePoint->GetPhysicalVolume();
      auto postVol = thePostPoint->GetPhysicalVolume();
      G4String preName = (preVol) ? preVol->GetName() : "OutOfWorld";
      G4String postName = (postVol) ? postVol->GetName() : "OutOfWorld";

      G4cout << "Step crossed a geometry boundary!" << G4endl;
      G4ThreeVector PPposition = thePrePoint->GetPosition();
      G4cout << "=== First step at geometry boundary ===" << G4endl;
      G4cout << "Volume: from " << preName << " to " << postName << G4endl;
      G4cout << "Position: " << PPposition << " mm" << G4endl;
      G4cout << "Time: " << thePrePoint->GetGlobalTime() << " ns" << G4endl;
    }

    const G4VProcess *proc = thePostPoint->GetProcessDefinedStep();
    G4int procID = -2; // default if process not recognized
    G4String processName;
    if(proc!=NULL){
      processName = proc->GetProcessName(); //Cerenkov, Scintilation and hadElastic, etc.
      if (processName == "Transportation") procID = -1; //Transportation in not in the process list
      else procID = WCSimEnumerations::ProcessTypeStringToEnum(proc->GetProcessName());
    }

    if ((processName != "Cerenkov" && processName != "hIoni" && processName != "Transportation" && processName != "muIoni") || isNewEvent || isNewVolume){
      G4int trackID = track->GetTrackID();
      G4int parentID = track->GetParentID();

      const G4ThreeVector position = thePostPoint->GetPosition();
      const G4ThreeVector preposition = thePrePoint->GetPosition();
      // Get momentum vector (has units: MeV/c)
      G4ThreeVector p_before = thePrePoint->GetMomentum();
      G4ThreeVector p_after = thePostPoint->GetMomentum();

      G4double time_ns = thePostPoint->GetGlobalTime();
      G4double pretime_ns = thePrePoint->GetGlobalTime();

      std::cout << "(" << track->GetDefinition()->GetPDGEncoding() << ", ID " << trackID << ") process " << proc->GetProcessType() << " :  " << processName << " , " << p_before.mag() << " MeV/c, " <<position.x()<<" "<<position.y()<<" "<<position.z() << " , " << time_ns << " ns" << std::endl;

      // Look for a secondary pi+ in the current step
      const std::vector<const G4Track*>* secondaries = aStep->GetSecondaryInCurrentStep();
      const G4Track* newPion = nullptr;
      G4double max_p_mag = -1.0;

      int N_pip = 0;
      int N_pim = 0;
      int N_muons = 0;
      int N_pi0 = 0;
      int N_other = 0;

      for (const auto* sec : *secondaries)
      {
        if (std::abs(sec->GetDefinition()->GetPDGEncoding())==211) //(sec->GetDefinition()->GetParticleName() == "pi+")
        {
          G4double p_mag = sec->GetMomentum().mag();
          if (p_mag > max_p_mag)
          {  
            max_p_mag = p_mag;
            newPion = sec;
          }
        }
        
      	G4String name = sec->GetDefinition()->GetParticleName();
      	if (name == "pi0") ++N_pi0;
      	if (sec->GetDefinition()->GetPDGCharge() == 0.) continue; // neutral particles can't emit Cherenkov
      	G4cout << "Secondary: " << sec->GetDefinition()->GetParticleName() << G4endl;
      	G4double beta = sec->GetVelocity() / CLHEP::c_light; // in mm/ns
        G4double n = 1.33; // refractive index for water

        if (name == "pi+") ++N_pip;
        if (name == "pi-") ++N_pim;
        else if (name == "mu+" || name == "mu-") ++N_muons;
        else if (beta <= 1.0 / n) continue; // below Cherenkov threshold, skip
        else ++N_other;
      }

      if (newPion)
      {
        p_after = newPion->GetMomentum();

        G4cout << "\n==== pi+Inelastic with outgoing pi+ (highest momentum) ====" << G4endl;
        G4cout << "Initial pi+ momentum: " << p_before.mag() << " MeV/c" << G4endl;
        G4cout << "Secondary pi+ momentum: " << p_after.mag() << " MeV/c" << G4endl;
        G4cout << "Lost momentum: " << p_before.mag() - p_after.mag() << " MeV/c" << G4endl;
      }

      // Extract direction components (unit vector)
      G4ThreeVector dir, predir;
      predir = (p_before.mag2() > 0) ? p_before.unit() : G4ThreeVector(0, 0, 0);
      dir = (p_after.mag2() > 0) ? p_after.unit() : G4ThreeVector(0, 0, 0);
      G4double dirX = dir.x();
      G4double dirY = dir.y();
      G4double dirZ = dir.z();
      G4double predirX = predir.x();
      G4double predirY = predir.y();
      G4double predirZ = predir.z();

      // Assign and fill
      float momMag, Ploss, premomMag, prePloss;
      premomMag = p_before.mag(); // already in MeV/c
      prePloss = 0.;
      momMag = p_after.mag();
      Ploss = p_before.mag() - p_after.mag();

      if (sctRunAction) {
        if (isNewEvent || isNewVolume) sctRunAction->FillTrackData(
          G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID(), trackID, parentID, procID, track->GetDefinition()->GetPDGEncoding(),
          pretime_ns, preposition.x() / 10., preposition.y() / 10., preposition.z() / 10.,
          premomMag, predirX, predirY, predirZ,
          0, 0, 0, 0,// do not keep post step, will be saved below if process is useful
          prePloss, N_pip, N_pim, N_muons, N_pi0, N_other,
          1//static_cast<int>(isNewEvent || isNewVolume)
        );

        if (processName != "Cerenkov" && processName != "hIoni" && processName != "Transportation" && processName != "muIoni") sctRunAction->FillTrackData(
          G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID(), trackID, parentID, procID, track->GetDefinition()->GetPDGEncoding(),
          time_ns, position.x() / 10., position.y() / 10., position.z() / 10.,
          premomMag, predirX, predirY, predirZ,
          momMag, dirX, dirY, dirZ,
          Ploss, N_pip, N_pim, N_muons, N_pi0, N_other,
          0//static_cast<int>(isNewVolume)
        );
      }
    }
  }

  G4ParticleDefinition *particleType = track->GetDefinition();
  if(particleType == G4OpticalPhoton::OpticalPhotonDefinition() && thePrePV && thePostPV){
    if( (thePrePV->GetName().find("MultiPMT") != std::string::npos) &&
	(thePostPV->GetName().find("vessel") != std::string::npos) )
      n_photons_through_mPMTLV++;
    
    if( (thePostPV->GetName().find("container") != std::string::npos) &&
	(thePrePV->GetName().find("vessel") != std::string::npos))
      n_photons_through_acrylic++;
    
    if( (thePrePV->GetName().find("container") != std::string::npos) )
      n_photons_through_gel++;

    if( (thePrePV->GetName().find("container") != std::string::npos) &&
	(thePostPV->GetName().find("inner") != std::string::npos) )
      n_photons_on_blacksheet++;

    if( (thePrePV->GetName().find("container") != std::string::npos) &&
	(thePostPV->GetName().find("pmt") != std::string::npos) )
      n_photons_on_smallPMT++;
	

    /*
    if( (thePrePV->GetName().find("pmt") != std::string::npos)){
      G4cout << "Photon between " << thePrePV->GetName() <<
	" and " << thePostPV->GetName() << " because " << 
	thePostPoint->GetProcessDefinedStep()->GetProcessName() << 
	" and boundary status: " <<  boundary->GetStatus() << " with track status " << track->GetTrackStatus() << G4endl;
      
	}*/



    if(track->GetTrackStatus() == fStopAndKill){
      if(boundary->GetStatus() == NoRINDEX){
	G4cout << "Optical photon is killed because of missing refractive index in either " << thePrePoint->GetMaterial()->GetName() << " or " << thePostPoint->GetMaterial()->GetName() <<
	  " (transition from " << thePrePV->GetName() << " to " << thePostPV->GetName() << ")" <<
	  " : could also be caused by Overlaps with volumes with logicalBoundaries." << G4endl;
	
      }
      /* Debug :  
      if( (thePrePV->GetName().find("PMT") != std::string::npos) ||
	  (thePrePV->GetName().find("pmt") != std::string::npos)){
	
	if(boundary->GetStatus() != StepTooSmall){
	  //	if(thePostPoint->GetProcessDefinedStep()->GetProcessName() != "Transportation")
	  G4cout << "Killed photon between " << thePrePV->GetName() <<
	    " and " << thePostPV->GetName() << " because " << 
	    thePostPoint->GetProcessDefinedStep()->GetProcessName() << 
	    " and boundary status: " <<  boundary->GetStatus() <<
	    G4endl;
	}
	}	*/
      
    }
  }



}


G4int WCSimSteppingAction::G4ThreeVectorToWireTime(G4ThreeVector *pos3d,
						    G4ThreeVector lArPos,
						    G4ThreeVector start,
						    G4int i)
{
  G4double x0 = start[0]-lArPos[0];
  G4double y0 = start[1]-lArPos[1];
//   G4double y0 = 2121.3;//mm
  G4double z0 = start[2]-lArPos[2];

  G4double dt=0.8;//mm
//   G4double midt = 2651.625;
  G4double pitch = 3;//mm
//   G4double midwir = 1207.10;
  G4double c45 = 0.707106781;
  G4double s45 = 0.707106781;

  G4double w1;
  G4double w2;
  G4double t;

//   G4double xField(0.);
//   G4double yField(0.);
//   G4double zField(0.);

//   if(detector->getElectricFieldDistortion())
//     {
//       Distortion(pos3d->getX(),
// 		 pos3d->getY());
      
//       xField = ret[0];
//       yField = ret[1];
//       zField = pos3d->getZ();
      
//       w1 = (int)(((zField+z0)*c45 + (x0-xField)*s45)/pitch); 
//       w2 = (int)(((zField+z0)*c45 + (x0+xField)*s45)/pitch); 
//       t = (int)(yField+1);

//       //G4cout<<" x orig "<<pos3d->getX()<<" y orig "<<(pos3d->getY()+y0)/dt<<G4endl;
//       //G4cout<<" x new "<<xField<<" y new "<<yField<<G4endl;
//     }
//   else 
//     {
      
  w1 = (int) (((pos3d->getZ()+z0)*c45 + (x0-pos3d->getX())*s45)/pitch); 
  w2 = (int)(((pos3d->getZ()+z0)*c45 + (x0+pos3d->getX())*s45)/pitch); 
  t  = (int)((pos3d->getY()+y0)/dt +1);
//     }

  if (i==0)
    return (int)w1;
  else if (i==1)
    return (int)w2;
  else if (i==2)
    return (int)t;
  else return 0;
} 


void WCSimSteppingAction::Distortion(G4double /*x*/,G4double /*y*/)
{
 
//   G4double theta,steps,yy,y0,EvGx,EvGy,EField,velocity,tSample,dt;
//   y0=2121.3;//mm
//   steps=0;//1 mm steps
//   tSample=0.4; //micros
//   dt=0.8;//mm
//   LiquidArgonMedium medium;  
//   yy=y;
//   while(y<y0 && y>-y0 )
//     {
//       EvGx=FieldLines(x,y,1);
//       EvGy=FieldLines(x,y,2);
//       theta=atan(EvGx/EvGy);
//       if(EvGy>0)
// 	{
// 	  x+=sin(theta);
// 	  y+=cos(theta);
// 	}
//       else
// 	{
// 	  y-=cos(theta);
// 	  x-=sin(theta);
// 	}
//       EField=sqrt(EvGx*EvGx+EvGy*EvGy);//kV/mm
//       velocity=medium.DriftVelocity(EField*10);// mm/microsec
//       steps+=1/(tSample*velocity);

//       //G4cout<<" step "<<steps<<" x "<<x<<" y "<<y<<" theta "<<theta<<" Gx "<<eventaction->Gx->Eval(x,y)<<" Gy "<<eventaction->Gy->Eval(x,y)<<" EField "<<EField<<" velocity "<<velocity<<G4endl;
//     }

//   //numbers
//   //EvGx=FieldLines(0,1000,1);
//   //EvGy=FieldLines(0,1000,2);
//   //EField=sqrt(EvGx*EvGx+EvGy*EvGy);//kV/mm
//   //velocity=medium.DriftVelocity(EField*10);// mm/microsec
//   //G4double quenching;
//   //quenching=medium.QuenchingFactor(2.1,EField*10);
//   //G4cout<<" Gx "<<EvGx<<" Gy "<<EvGy<<" EField "<<EField<<" velocity "<<velocity<<" quenching "<<quenching<<G4endl;


  //ret[0]=5;
//   if(yy>0)
//     ret[1]=2*y0/dt -steps;
//   else
//     ret[1]=steps; 
}


double WCSimSteppingAction::FieldLines(G4double /*x*/,G4double /*y*/,G4int /*coord*/)
{ //0.1 kV/mm = field
  //G4double Radius=302;//mm
//   G4double Radius=602;//mm
//   if(coord==1) //x coordinate
//     return  (0.1*(2*Radius*Radius*x*abs(y)/((x*x+y*y)*(x*x+y*y))));
//   else //y coordinate
//     return 0.1*((abs(y)/y)*(1-Radius*Radius/((x*x+y*y)*(x*x+y*y))) + abs(y)*(2*Radius*Radius*y/((x*x+y*y)*(x*x+y*y))));
  return 0;
}
