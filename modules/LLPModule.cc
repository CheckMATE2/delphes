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

  fOutputArrayAll = ExportArray(GetString("OutputArrayAll", "All"));
  fOutputArrayMothers = ExportArray(GetString("OutputArrayMothers", "Mothers"));

}

//------------------------------------------------------------------------------

void LLPModule::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------


float decRad(Candidate* c) {
    double x = c->Position.X();
    double y = c->Position.Y();
    double z = c->Position.Z();
    return sqrt(x*x+y*y+z*z);
}

void LLPModule::storeCandidate(Candidate* candidate) {
    //std::cout << "Start " << fOutputArrayAll->GetEntries() << " and c " << candidate << "(" << candidate->PID << ") decays into " << candidate->D1 << "  " << candidate->D2;
    Candidate* currPart1;
    Candidate* currPart2;
    if (candidate->D1 >= 0) {
	currPart1 = static_cast<Candidate*>(fInputArray->At(candidate->D1));
	currPart2 = static_cast<Candidate*>(fInputArray->At(candidate->D2));
	//std::cout << " = " << currPart1 << " (" << currPart1->PID << ") and " << currPart2 << " (" << currPart2->PID << ")" << std::endl;
    }
    int motherIndex = fOutputArrayAll->GetEntries();
    bool onePartDecay = candidate->D1 == candidate->D2;
    fOutputArrayAll->Add(candidate);
    double y;
    //std::cin >> y;
    if (candidate->D1 >= 0) {
	// Change the mothers index to the new index it gets through the reduced size of the output list
	//std::cout << "Changed " << currPart1 << " mother to " << motherIndex << std::endl;
	currPart1->M1 = motherIndex;
	storeCandidate(currPart1);
	candidate->D1 = motherIndex+1;
	//std::cout << "Changed " << candidate << " D1 to " << motherIndex+1 << std::endl;
	if (onePartDecay) // Indices are equal if it is a 'transisition' P -> P
	    candidate->D2 = motherIndex+1;
    }
    
    // as daughter1 was added next, the new D1 index is motherIndex+1    
    int currIndex = fOutputArrayAll->GetEntries();
    if (!onePartDecay && candidate->D2 >= 0) {	
	currPart2->M1 = motherIndex;
	//std::cout << "Changed " << currPart2 << " mother to " << motherIndex << std::endl;
	storeCandidate(currPart2);
	//daughter2 was added at the  currIndex
	candidate->D2 = currIndex;
	//std::cout << "Changed " << candidate << " D1 to " << currIndex << std::endl;
    }
    //std::cout << "All.size() end " << fOutputArrayAll->GetEntries() << " " << candidate->D1 << "  " << candidate->D2 << std::endl;
    double x;
    //std::cin >> x;
}


void LLPModule::Process()
{
  Candidate *candidate;
  TLorentzVector candidatePosition, candidateMomentum;

  // loop over all input candidates
  fItInputArray->Reset();

  // all mothers appear for several daughters and we need to avoid double counting
  std::deque <Candidate*> storedMothers;
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    candidatePosition = candidate->Position;
    candidateMomentum = candidate->Momentum;

    double d = decRad(candidate);
    if ( (fMinRadius <= d) && (d <= fMaxRadius) ) {	
	Candidate* mother = static_cast<Candidate*>(fInputArray->At(candidate->M1));
	
	if (fPdgCodes.size() == 0 || find(fPdgCodes.begin(), fPdgCodes.end(), mother->PID) != fPdgCodes.end()) {
	    if (find(storedMothers.begin(), storedMothers.end(), mother) == storedMothers.end()) {      
		storedMothers.push_back(mother);
		fOutputArrayMothers->Add(mother);
	    }
	}
    }
  }
    //std::cout << "Found " << storedMothers.size() << " mothers" << std::endl;
    for(int i =0; i < storedMothers.size(); i++) {
	//std::cout << "Mother " << i << ": " << storedMothers[i] << std::endl;
	storeCandidate(storedMothers[i]);
    }
}

//------------------------------------------------------------------------------
