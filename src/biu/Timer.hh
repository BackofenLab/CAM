// $Id: Timer.hh,v 1.1 2013/04/10 21:32:37 naharf Exp $
#ifndef BIU_TIMER_HH_
#define BIU_TIMER_HH_

#include <ctime>

		/**
		 * Timer class to measure runtime in miliseconds.
		 *
		 * @author Martin Mann <mmann@@informatik.uni-freiburg.de>
		 */
	class Timer {
		private:
				//! starting time
			clock_t t0;
		public:
				//! Sets starting time.
			void start(void){
				t0 = clock();
			}
				//! Returns time consumption in miliseconds until now from last
				//! start() call on.
			double stop(void) {
				return (static_cast<double>(clock()-t0) / CLOCKS_PER_SEC) * 1000.0;
			}
	};
	
#endif /*TIMER_HH_*/
