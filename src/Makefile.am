
# additional libs and include paths necessary to compile the C++ programs
AM_CXXFLAGS = 	-I$(top_srcdir)/src \
				-I$(top_builddir)/src/ \
				@AM_CXXFLAGS@

AM_LDFLAGS  = @AM_LDFLAGS@ 
LIBS        = @LIBS@ 

bin_PROGRAMS = cam molKelkule

#############   CYCLIC ITS ATOM MAPPING   ##################################

cam_SOURCES =	\
					AlternatingCycle.h \
					AlternatingCycle.cpp \
					biu/OptionParser.hh \
					biu/OptionParser.cc \
					biu/Timer.hh \
					ChargeChange.h \
					ChargeChange.cpp \
					ConnectedCoverage.h \
					ConnectedCoverage.cpp \
					EdgeDegree.h \
					EdgeDegree.cpp \
					GC_StlSetRangeIterator.h \
					GC_StlSetRangeIterator.cpp \
					ITS.h \
					ITS.cpp \
					ITSAtomLabel.h \
					ITSAtomLabel.cpp \
					ITS_CSP.h \
					ITS_CSP.cpp \
					ITS_Graph_Interface.h \
					ITS_Graph_Interface.cpp \
					ValenceMatrix.h \
					ValenceMatrix.cpp \
					ReactionGraph.h \
					ReactionGraph.cpp \
					MatrixSparse.hh \
					MatrixSparse.icc \
					MinEdgeValence.h \
					MinEdgeValence.cpp \
					NeighborhoodCount.h \
					NeighborhoodCount.cpp \
					ReactionSMILES_Extractor.h \
					ReactionSMILES_Extractor.cpp \
					cyclicAtomMapping.cpp				

#############   AMBIGUOUS RING BOND ENUMERATION   ############################

molKelkule_SOURCES =	\
					biu/OptionParser.hh \
					biu/OptionParser.cc \
					molKekule.cc

####################################################################

valgrind: $(bin_PROGRAMS) 
	valgrind --leak-check=full $< 
