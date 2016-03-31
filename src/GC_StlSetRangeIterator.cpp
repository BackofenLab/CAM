/*
 *  Main authors:
 *     Martin Mann http://www.bioinf.uni-freiburg.de/~mmann/
 *
 *  Contributing authors:
 *     Sebastian Will http://www.bioinf.uni-freiburg.de/~will/
 *
 *  Copyright:
 *     Martin Mann, 2007
 *
 *  This file is part of the CPSP-tools package:
 *     http://www.bioinf.uni-freiburg.de/sw/cpsp/
 *
 *  See the file "LICENSE" for information on usage and
 *  redistribution of this file, and for a
 *     DISCLAIMER OF ALL WARRANTIES.
 *
 */

#include "GC_StlSetRangeIterator.h"
  	
  	GC_StlSetRangeIterator::GC_StlSetRangeIterator() :
  		data(NULL), noFurtherRange(true)
  	{
  		getNextRange();
  	}
  		
	GC_StlSetRangeIterator::GC_StlSetRangeIterator(const std::set<int>* data_) :
		data(data_), noFurtherRange(false)
	{
		if (data != NULL) 
			actElem = data->begin();
		getNextRange();
	}
	
	GC_StlSetRangeIterator::~GC_StlSetRangeIterator()
	{
	}
	
	void
	GC_StlSetRangeIterator::getNextRange() {
		if (data==NULL || actElem == data->end()) {
			noFurtherRange = true;
			return;
		}
			// find next range
		nextMin = *actElem;
		nextMax = nextMin;
			// build up new upper bound until end of set reached or gap in 
			// sequence
		while ( (++actElem != data->end()) && (*actElem == (nextMax+1))) {
			nextMax++;
		}
	}
