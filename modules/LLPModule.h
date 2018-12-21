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

#ifndef LLPModule_h
#define LLPModule_h

/** \class LLPModule
 *
 *  Selects candidates and their decay trees from the InputArray with a decay vertex within a certain region
 *
 *  \author D. Dercks, DESY
 *
 */

#include "classes/DelphesModule.h"
#include "classes/DelphesClasses.h"
#include <deque>

class TObjArray;
class DelphesFormula;

class LLPModule: public DelphesModule
{
public:

  LLPModule();
  ~LLPModule();

  void Init();
  void Process();
  void Finish();
  void storeCandidate(Candidate* candidate); //! Stores candidate and its entire decay tree

private:
  
  Double_t fMinRadius; //! Minimal radius (in mm) for particle->Position to be considered
  Double_t fMaxRadius; //! Maximal radius (in mm) for particle->Position to be considered

  std::deque <Int_t> fPdgCodes;
  
  DelphesFormula *fFormula; //!

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  TObjArray *fOutputArrayAll; //!
  TObjArray *fOutputArrayMothers; //!
  TObjArray *fOutputArrayTracks; //!
  TObjArray *fOutputArrayStableDaughters; //!
  

  ClassDef(LLPModule, 1)
};

#endif
