/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/** \class LLPModule
 *
 *  Selects candidates from the InputArray according to the efficiency formula.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/LLPModule.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TMath.h"
#include "TString.h"
#include "TFormula.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

#include <algorithm> 
#include <stdexcept>
#include <iostream>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------

LLPModule::LLPModule() :
  fFormula(0), fItInputArray(0)
{
  fFormula = new DelphesFormula;
}

//------------------------------------------------------------------------------

LLPModule::~LLPModule()
{
  if(fFormula) delete fFormula;
}

//------------------------------------------------------------------------------

void LLPModule::Init()
{
  // read parameters


  fFormula->Compile(GetString("EfficiencyFormula", "0.4"));

  ExRootConfParam param = GetParam("PDGCodes");
  Long_t i, size;
  fPdgCodes.clear();

  size = param.GetSize();
  for(i = 0; i < size; ++i)
  {
    fPdgCodes.push_back(param[i].GetInt());
  }
  
  // import input array(s)

  fInputArray = ImportArray(GetString("InputArray", "FastJetFinder/jets"));
  fItInputArray = fInputArray->MakeIterator();

  fMinRadius = GetDouble("MinRadius", 100.0);
  fMaxRadius = GetDouble("MaxRadius", 10000.0);

  // create output array(s)

  fOutputArray = ExportArray(GetString("OutputArray", "jets"));

}

//------------------------------------------------------------------------------

void LLPModule::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void LLPModule::printCandidate(Candidate* currPart, std::string currSpaces) {
    std::cout << currSpaces <<  currPart << "  " << currPart->Status << "  " << currPart->PID << " " << currPart->Momentum.E() << " " << currPart->InitialPosition.T() << "/" << currPart->InitialPosition.X()<< "/" << currPart->InitialPosition.Y()<< "/" << currPart->InitialPosition.Z() << " -> " << currPart->Position.T() << "/" << currPart->Position.X()<< "/" << currPart->Position.Y()<< "/" << currPart->Position.Z() << std::endl;
    if (currPart->D1 > 0) {
        Candidate* currPart1 = static_cast<Candidate*>(fInputArray->At(currPart->D1));	    
	printCandidate(currPart1, currSpaces+"  ");
    }
    if (currPart->D2 > 0 and currPart->D2 != currPart->D1) {
        Candidate* currPart2 = static_cast<Candidate*>(fInputArray->At(currPart->D2));	    
	printCandidate(currPart2, currSpaces+"  ");
    }
}

float decRad(Candidate* c) {
    double x = c->Position.X();
    double y = c->Position.Y();
    double z = c->Position.Z();
    return sqrt(x*x+y*y+z*z);
}



void LLPModule::Process()
{
  Candidate *candidate;
  TLorentzVector candidatePosition, candidateMomentum;

  // loop over all input candidates
  fItInputArray->Reset();
  std::cout << " New Event " << std::endl;
  int counter = 0;
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    candidatePosition = candidate->Position;
    candidateMomentum = candidate->Momentum;

    double d = decRad(candidate);
    if ( (fMinRadius <= d) && (d <= fMaxRadius) ) {	
	std::cout << "Particle " << candidate->PID << " has nontrivial decay vertex! (r=" << d << ")" << std::endl;
	Candidate* mother = static_cast<Candidate*>(fInputArray->At(candidate->M1));
	if (fPdgCodes.size() == 0 || find(fPdgCodes.begin(), fPdgCodes.end(), mother->PID) != fPdgCodes.end()) {
	    std::cout << "Mother " << mother->PID << std::endl;
	    Candidate* dau1 = static_cast<Candidate*>(fInputArray->At(mother->D1));
	    Candidate* dau2 = static_cast<Candidate*>(fInputArray->At(mother->D2));	    
	    std::cout << "Daughter1 " << dau1->PID << " has radius " << decRad(dau1) << std::endl;
	    std::cout << "Daughter2 " << dau2->PID << " has radius " << decRad(dau2) << std::endl;
	}
	std::cout << "  " << std::endl;
    }


	
	/*
	std::cout << "Found Chargino" << std::endl;
	if (candidate->D1 > 0) {
	    std::cout << "  E " << candidate->Momentum.E() << std::endl;
	    std::cout << "  It decays into " << candidate->D1 << "  " << candidate->D2 << std::endl;
	    std::cout << candidate->D0 << "  " << candidate->DZ << std::endl;
	    std::cout << " Its initial position " 

	    Candidate* D1 =  static_cast<Candidate*>(fInputArray->At(candidate->D1));
	    Candidate* D2 =  static_cast<Candidate*>(fInputArray->At(candidate->D2));
	    std::cout << "   D1PID, x " << D1->PID << "  " << D1->Momentum.E() << "  " << D1->InitialPosition.X()<< "/" << D1->InitialPosition.Y()<< "/" << D1->InitialPosition.Z() << " its final position " << D1->Position.X()<< "/" << D1->Position.Y()<< "/" << D1->Position.Z() << std::endl;
	    std::cout << "   D2PID, x " << D2->PID << "  " << D2->Momentum.E() << "  " << D2->InitialPosition.X()<< "/" << D2->InitialPosition.Y()<< "/" << D2->InitialPosition.Z() << " its final position " << D2->Position.X()<< "/" << D2->Position.Y()<< "/" << D2->Position.Z() << std::endl;
	   	
	}
    }
            */

    fOutputArray->Add(candidate);

  }
}

//------------------------------------------------------------------------------
