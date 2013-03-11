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


#include "ReferenceAligner.h"

namespace Bavieca {

// constructor
ReferenceAligner::ReferenceAligner(LexiconManager *lexiconManager)
{
	m_lexiconManager = lexiconManager;
}

// destructor
ReferenceAligner::~ReferenceAligner()
{
}


// align a recognition path against a sequence of lexical units
// Note: it does not consider as errors deletions after all the words in the hypothesis have been aligned
// - The hypothesis can be empty
// - The reference must contain at least one element
ReferenceAlignment *ReferenceAligner::align(BestPath *bestPath, ReferenceText *referenceText) {

   // (1) convert the best path to a sequence of lexical unit ids (only standard lexical units)
	LBestPathElement *lBestPathElement = bestPath->getBestPathElements();
	VBestPathElement hypothesis;
	for(LBestPathElement::iterator it = lBestPathElement->begin() ; it != lBestPathElement->end() ; ++it) {
		if (m_lexiconManager->isStandard((*it)->lexUnit->iLexUnit)) {
			hypothesis.push_back(*it);	
			//printf("%s\n",m_lexiconManager->getStrLexUnit((*it)->lexUnit->iLexUnit));
		}		
	}
   
   // if there are no elements in the reference return error
   if (referenceText->empty()) {
      return NULL;
   }
 
   ReferenceAlignment *referenceAlignment = new ReferenceAlignment(m_lexiconManager);
   
   // (3) perform the alignment   
   int iCorrects = 0;
   int iInsertions = 0;
   int iDeletions = 0;
   int iSubstitutions = 0;
   
   int j;
   
   // create the grid
   int rows = (int)(hypothesis.size()+1);
   int cols = (int)(referenceText->size()+1);

   int i,pen;
   GridElement **grid=(GridElement**)malloc(rows*sizeof(GridElement*));
   grid[0]=(GridElement*)malloc(rows*cols*sizeof(GridElement));
   for(i=1;i<rows;i++) 
      grid[i]=grid[i-1]+cols;

   i=j=0;
   pen=0;
   grid[i][j].score=pen;
   grid[i][j].prev=NULL;
   i++; j++;
   pen=WA_INS_PEN;
   for(;i<rows;i++) {
      grid[i][0].score=pen;
      pen+=WA_INS_PEN;
      grid[i][0].prev=&(grid[i-1][0]);
      grid[i][0].ref=-1;
      grid[i][0].tra=i-1;//NULL;
      grid[i][0].flag=WA_INS;
   }
   pen=WA_DEL_PEN;
   for(;j<cols;j++) {
      grid[0][j].score=pen;
      pen+=WA_DEL_PEN;
      grid[0][j].prev=&(grid[0][j-1]);
      grid[0][j].ref=j-1;//NULL;
      grid[0][j].tra=-1;
      grid[0][j].flag=WA_DEL;
   }

   //Update the grid
   int w0 = -1;
   int w = -1; 
   GridElement *g= NULL;
   int del,ins,sub,cor;

   // step 1: fill grid    
   for(i=1;i<rows;i++) {                  /* each row */
      w0=i-1;
      for(j=1;j<cols;j++) {
         w=j-1;
         grid[i][j].tra=w0;
         grid[i][j].ref=w;
         if (i == rows - 1) {  //deletions are not such once the partial hypothesis is totally aligned
            del=grid[i][j-1].score;
         } else {
            del=grid[i][j-1].score+WA_DEL_PEN;
         }
         ins=grid[i-1][j].score+WA_INS_PEN; 
         cor= (hypothesis[w0]->lexUnit->iLexUnit == (*referenceText)[w]->iLexUnit);
         sub=grid[i-1][j-1].score+(cor?0:WA_SUB_PEN);
         if(sub<del && sub<ins) {
            grid[i][j].score=sub;
            grid[i][j].prev=&(grid[i-1][j-1]);
            grid[i][j].flag=cor?WA_COR:WA_SUB;
         } else if(ins<del) {
            grid[i][j].score=ins;
            grid[i][j].prev=&(grid[i-1][j]);
            grid[i][j].flag=WA_INS;
         } else {       // del is the default if ins==del 
            grid[i][j].score=del;
            grid[i][j].prev=&(grid[i][j-1]);
            grid[i][j].flag=WA_DEL;
            grid[i][j].bookword=i;
         }   
      }
   }

   // step 2: trace back through grid to get best minimal path    
   g=&(grid[rows-1][cols-1]);
   if (hypothesis.empty() == false) {
   	printf("Alignment score: %d\n",g->score);
   } else {
   	printf("Alignment score: 0\n");
   }
   g->next=NULL;
   while(g->prev) {     //connect forward links       
      g->prev->next=g;     
      g=g->prev;
   }

   g=grid[0][0].next;
   
   cout << "hypothesis: " << hypothesis.size() << endl;
   
   while(g) {
      if(g->flag&WA_DEL) {
         ++iDeletions;
         //m_wordsIncorrect.push_back(reference[g->ref]); 
         assert((g->ref >= 0) && (g->ref < (int)(*referenceText).size()));
         referenceAlignment->addElement(ALIGNMENT_EVENT_DELETION,
         										g->ref,
         										(*referenceText)[g->ref]->iLexUnit,
         										-1,
         										m_lexiconManager->m_lexUnitUnknown->iLexUnit,
         										-1,
         										-1,
         										-1);
         
      } else if(g->flag&WA_INS) {
         ++iInsertions;
         //printf("Inserted word: %s\n",hypothesis[g->tra-1].c_str());
         //printf("hyp index = %d\n",g->tra);
         assert((g->tra >= 0) && (g->tra < (int)hypothesis.size()));
         referenceAlignment->addElement(ALIGNMENT_EVENT_INSERTION,
         										-1,
         										m_lexiconManager->m_lexUnitUnknown->iLexUnit,
         										g->tra,
         										hypothesis[g->tra]->lexUnit->iLexUnit,
         										hypothesis[g->tra]->iFrameStart,
         										hypothesis[g->tra]->iFrameEnd,
         										hypothesis[g->tra]->fScoreConfidence);
         
      } else if(g->flag&WA_SUB) {
         ++iSubstitutions;
         //printf("Substituted word: %s -> %s\n",reference[g->ref]->strWord.c_str(),hypothesis[g->tra].c_str());
         //m_wordsIncorrect.push_back(reference[g->ref]);
         //iLastAlignedWordInReference = g->ref;
         assert((g->tra >= 0) && (g->tra < (int)hypothesis.size()));
         assert((g->ref >= 0) && (g->ref < (int)(*referenceText).size()));
         referenceAlignment->addElement(ALIGNMENT_EVENT_SUBSTITUTION,
         										g->ref,
         										(*referenceText)[g->ref]->iLexUnit,
         										g->tra,
         										hypothesis[g->tra]->lexUnit->iLexUnit,
         										hypothesis[g->tra]->iFrameStart,
         										hypothesis[g->tra]->iFrameEnd,
         										hypothesis[g->tra]->fScoreConfidence);
      } else if(g->flag&WA_COR) {
         ++iCorrects;
         //printf("Correct word: %s -> %s\n",reference[g->ref]->strWord.c_str(),hypothesis[g->tra].c_str());
         assert((g->tra >= 0) && (g->tra < (int)hypothesis.size()));
         assert((g->ref >= 0) && (g->ref < (int)(*referenceText).size()));
         referenceAlignment->addElement(ALIGNMENT_EVENT_CORRECT,
         										g->ref,
         										(*referenceText)[g->ref]->iLexUnit,
         										g->tra,
         										hypothesis[g->tra]->lexUnit->iLexUnit,
         										hypothesis[g->tra]->iFrameStart,
         										hypothesis[g->tra]->iFrameEnd,
         										hypothesis[g->tra]->fScoreConfidence);

      } else {
      	assert(0);
      }
      g=g->next;
   }
 
	// clean-up
   if(grid) {
      if(grid[0]) 
      	free(grid[0]);
      free(grid);
   } 

	
	return referenceAlignment;	
}

};	// end-of-namespace

