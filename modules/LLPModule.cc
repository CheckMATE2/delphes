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

  fIntParam = GetInt("IntParam", 10);

  fDoubleParam = GetDouble("DoubleParam", 1.0);

  fFormula->Compile(GetString("EfficiencyFormula", "0.4"));

  ExRootConfParam param = GetParam("ArrayParam");
  Long_t i, size;
  fArrayParam.clear();

  size = param.GetSize();
  for(i = 0; i < size; ++i)
  {
    fArrayParam.push_back(param[i].GetDouble());
  }
  
  // import input array(s)

  fInputArray = ImportArray(GetString("InputArray", "FastJetFinder/jets"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array(s)

  fOutputArray = ExportArray(GetString("OutputArray", "jets"));

}

//------------------------------------------------------------------------------

void LLPModule::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void LLPModule::Process()
{
  Candidate *candidate;
  TLorentzVector candidatePosition, candidateMomentum;

  // loop over all input candidates
  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    candidatePosition = candidate->Position;
    candidateMomentum = candidate->Momentum;

    if (candidate->PID == 1000024) {
	std::cout << "Found Chargino" << std::endl;
	if (candidate->D1 > 0) {
	    std::cout << "  It decays into " << candidate->D1 << "  " << candidate->D2 << std::endl;
	    GenParticle* D1 =  static_cast<GenParticle *>(fInputArray->At(candidate->D1));
	    GenParticle* D2 =  static_cast<GenParticle *>(fInputArray->At(candidate->D2));
	    std::cout << "   D1PID, x " << D1->PID << "  " << D1->PT << " " << D1->X << "/" << D1->T << std::endl;
	    std::cout << "   D2PID, x " << D2->PID << "  " << D2->PT << " " << D2->X << "/" << D1->T << std::endl;
	}
    }
    fOutputArray->Add(candidate);

  }
}

//------------------------------------------------------------------------------
