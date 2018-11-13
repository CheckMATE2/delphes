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
    fOutputArrayAll->Add(candidate);
    if (candidate->D1 >= 0) {
	Candidate* currPart1 = static_cast<Candidate*>(fInputArray->At(candidate->D1));
	storeCandidate(currPart1);
    }
    if (candidate->D2 >= 0) {
	Candidate* currPart2 = static_cast<Candidate*>(fInputArray->At(candidate->D2));
	storeCandidate(currPart2);
    }
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
		storeCandidate(mother);
		storedMothers.push_back(mother);
		fOutputArrayMothers->Add(mother);
	    }
	}
    }

  }
}

//------------------------------------------------------------------------------
