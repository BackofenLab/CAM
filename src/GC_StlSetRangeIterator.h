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

#ifndef GC_STLSETRANGEITERATOR_HH_
#define GC_STLSETRANGEITERATOR_HH_


#include <set>
#include <iostream>
  	
  		/**
  		 * Provides a constant Gecode RangeIterator of a std::set<int> that 
  		 * calculates the ranges on demand.
  		 */
	class GC_StlSetRangeIterator
	{
	private:
		const std::set<int>* data;
		std::set<int>::const_iterator actElem;
		
		bool noFurtherRange;
		
		int nextMin, nextMax;
		
			//! searchs for the next range and sets the inner members
		void getNextRange();
	public:
		GC_StlSetRangeIterator();
		GC_StlSetRangeIterator(const std::set<int>* data_);
		virtual ~GC_StlSetRangeIterator();
		
		void init(const std::set<int>* const data_) {
			data = data_;
			noFurtherRange = false;
			if (data != NULL) 
				actElem = data->begin();
			getNextRange();
		}
		
		bool operator()(void) const { return !noFurtherRange; }
		void operator++(void) { getNextRange(); }
		
		int min(void) const { return nextMin; }
		int max(void) const { return nextMax; }
		unsigned int width(void) { return nextMax-nextMin+1; }
		
		
	};
	

#endif /*GC_STLSETRANGEITERATOR_HH_*/
