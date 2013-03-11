/*---------------------------------------------------------------------------------------------*
 * Copyright (C) 2012 Daniel Bolaños - www.bltek.com - Boulder Language Technologies           *
 *                                                                                             *
 * www.bavieca.org is the website of the Bavieca Speech Recognition Toolkit                    *
 *                                                                                             *
 * Licensed under the Apache License, Version 2.0 (the "License");                             *
 * you may not use this file except in compliance with the License.                            *
 * You may obtain a copy of the License at                                                     *
 *                                                                                             *
 *         http://www.apache.org/licenses/LICENSE-2.0                                          *
 *                                                                                             *
 * Unless required by applicable law or agreed to in writing, software                         *
 * distributed under the License is distributed on an "AS IS" BASIS,                           *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.                    *
 * See the License for the specific language governing permissions and                         *
 * limitations under the License.                                                              *
 *---------------------------------------------------------------------------------------------*/


#include "TextAligner.h"

namespace Bavieca {

// constructor
TextAligner::TextAligner(LexiconManager *lexiconManager)
{
	m_lexiconManager = lexiconManager;
}

// destructor
TextAligner::~TextAligner()
{
}

// align two sequences of lexical units
TextAlignment *TextAligner::align(VLexUnit &vLexUnitHyp, VLexUnit &vLexUnitRef, bool bPenalizeTrailingDeletions) {

   // if there are no elements in the reference return error
   if (vLexUnitRef.empty()) {
      return NULL;
   }
 
   TextAlignment *textAlignment = new TextAlignment(m_lexiconManager);
   
   // (3) perform the alignment   
   int iCorrects = 0;
   int iInsertions = 0;
   int iDeletions = 0;
   int iSubstitutions = 0;
   
   int j;
   
   // create the grid
   int iRows = (int)(vLexUnitHyp.size()+1);
   int iColumns = (int)(vLexUnitRef.size()+1);

   int i,pen;
   TAGridElement **grid = new TAGridElement*[iRows];
   grid[0]= new TAGridElement[iRows*iColumns];
   for(i=1;i<iRows;i++) 
      grid[i]=grid[i-1]+iColumns;

   i=j=0;
   pen=0;
   grid[i][j].iScore=pen;
   grid[i][j].prev=NULL;
   i++; j++;
   pen=TEXT_ALIGNMENT_PENALTY_INSERTION;
   for(;i<iRows;i++) {
      grid[i][0].iScore=pen;
      pen+=TEXT_ALIGNMENT_PENALTY_INSERTION;
      grid[i][0].prev=&(grid[i-1][0]);
      grid[i][0].iRef=-1;
      grid[i][0].iTra=i-1;
      grid[i][0].iFlag=TEXT_ALIGNMENT_EVENT_INSERTION;
   }
   pen=TEXT_ALIGNMENT_PENALTY_DELETION;
   for(;j<iColumns;j++) {
      grid[0][j].iScore=pen;
      pen+=TEXT_ALIGNMENT_PENALTY_DELETION;
      grid[0][j].prev=&(grid[0][j-1]);
      grid[0][j].iRef=j-1;
      grid[0][j].iTra=-1;
      grid[0][j].iFlag=TEXT_ALIGNMENT_EVENT_DELETION;
   }

   //Update the grid
   int iW0 = -1;
   int iW = -1; 
   TAGridElement *g= NULL;
   int del,ins,sub,cor;

   // step 1: fill grid    
   for(i=1;i<iRows;i++) {                  /* each row */
      iW0=i-1;
      for(j=1;j<iColumns;j++) {
         iW=j-1;
         grid[i][j].iTra=iW0;
         grid[i][j].iRef=iW;
         // deletions are not considered as such once the partial hypothesis is totally aligned
         if ((i == iRows - 1) && (bPenalizeTrailingDeletions == false)) {  
            del=grid[i][j-1].iScore;
         } 
         // always penalize deletions
         else {
            del=grid[i][j-1].iScore+TEXT_ALIGNMENT_PENALTY_DELETION;
         }
         ins=grid[i-1][j].iScore+TEXT_ALIGNMENT_PENALTY_INSERTION; 
         cor= (vLexUnitHyp[iW0]->iLexUnit == vLexUnitRef[iW]->iLexUnit);
         sub=grid[i-1][j-1].iScore+(cor?0:TEXT_ALIGNMENT_PENALTY_SUBSTITUTION);
         if(sub<del && sub<ins) {
            grid[i][j].iScore=sub;
            grid[i][j].prev=&(grid[i-1][j-1]);
            grid[i][j].iFlag=cor?TEXT_ALIGNMENT_EVENT_CORRECT:TEXT_ALIGNMENT_EVENT_SUBSTITUTION;
         } else if(ins<del) {
            grid[i][j].iScore=ins;
            grid[i][j].prev=&(grid[i-1][j]);
            grid[i][j].iFlag=TEXT_ALIGNMENT_EVENT_INSERTION;
         } else {       // del is the default if ins==del 
            grid[i][j].iScore=del;
            grid[i][j].prev=&(grid[i][j-1]);
            grid[i][j].iFlag=TEXT_ALIGNMENT_EVENT_DELETION;
            grid[i][j].iBookWord=i;
         }   
      }
   }

   // step 2: trace back through grid to get best minimal path    
   g=&(grid[iRows-1][iColumns-1]);
   textAlignment->setScore(g->iScore);
   /*if (vLexUnitHyp.empty() == false) {
   	printf("Alignment score: %d\n",g->iScore);
   } else {
   	printf("Alignment score: 0\n");
   }*/
   g->next=NULL;
   while(g->prev) {     //connect forward links       
      g->prev->next=g;     
      g=g->prev;
   }

   g=grid[0][0].next;
   
   //int iLastAlignedWordInReference = -1;
   
   //printf("hypothesis: %d\n",vLexUnitHyp.size());
   
   while(g) {
      if(g->iFlag&TEXT_ALIGNMENT_EVENT_DELETION) {
         ++iDeletions;
         //printf("deleted word: %s\n",m_lexiconManager->getStrLexUnit(vLexUnitRef[g->iRef]->iLexUnit));
         //m_wordsIncorrect.push_back(reference[g->iRef]); 
         assert((g->iRef >= 0) && (g->iRef < (int)vLexUnitRef.size()));
         textAlignment->addElement(TEXT_ALIGNMENT_EVENT_DELETION,
         										g->iRef,
         										vLexUnitRef[g->iRef]->iLexUnit,
         										-1,
         										m_lexiconManager->m_lexUnitUnknown->iLexUnit);
         
      } else if(g->iFlag&TEXT_ALIGNMENT_EVENT_INSERTION) {
         ++iInsertions;
         //printf("inserted word: %s\n",m_lexiconManager->getStrLexUnit(vLexUnitHyp[g->iTra]->iLexUnit));
         //printf("hyp index = %d\n",g->iTra);
         assert((g->iTra >= 0) && (g->iTra < (int)vLexUnitHyp.size()));
         textAlignment->addElement(TEXT_ALIGNMENT_EVENT_INSERTION,
         										-1,
         										m_lexiconManager->m_lexUnitUnknown->iLexUnit,
         										g->iTra,
         										vLexUnitHyp[g->iTra]->iLexUnit);
         
      } else if(g->iFlag&TEXT_ALIGNMENT_EVENT_SUBSTITUTION) {
         ++iSubstitutions;
         //printf("substituted word: %s -> %s\n",m_lexiconManager->getStrLexUnit(vLexUnitRef[g->iRef]->iLexUnit),m_lexiconManager->getStrLexUnit(vLexUnitHyp[g->iTra]->iLexUnit));
         //m_wordsIncorrect.push_back(reference[g->iRef]);
         //iLastAlignedWordInReference = g->iRef;
         assert((g->iTra >= 0) && (g->iTra < (int)vLexUnitHyp.size()));
         assert((g->iRef >= 0) && (g->iRef < (int)vLexUnitRef.size()));
         textAlignment->addElement(TEXT_ALIGNMENT_EVENT_SUBSTITUTION,
         										g->iRef,
         										vLexUnitRef[g->iRef]->iLexUnit,
         										g->iTra,
         										vLexUnitHyp[g->iTra]->iLexUnit);
      } else if(g->iFlag&TEXT_ALIGNMENT_EVENT_CORRECT) {
         ++iCorrects;
         //printf("correct word: %s -> %s\n",m_lexiconManager->getStrLexUnit(vLexUnitRef[g->iRef]->iLexUnit),m_lexiconManager->getStrLexUnit(vLexUnitHyp[g->iTra]->iLexUnit));
         assert((g->iTra >= 0) && (g->iTra < (int)vLexUnitHyp.size()));
         assert((g->iRef >= 0) && (g->iRef < (int)vLexUnitRef.size()));
         textAlignment->addElement(TEXT_ALIGNMENT_EVENT_CORRECT,
         										g->iRef,
         										vLexUnitRef[g->iRef]->iLexUnit,
         										g->iTra,
         										vLexUnitHyp[g->iTra]->iLexUnit);

      } else {
      	assert(0);
      }
      g=g->next;
   }
 
	// clean-up
	delete [] grid[0];
	delete [] grid;
  
   return textAlignment;
}

};	// end-of-namespace



