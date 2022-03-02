#	Local configuration for HySGen libraries. HySGen's version of Snap's functions and classes.
LOC_EXGLIB = local_snap/$(GLIB)
LOC_EXSNAP = local_snap/$(SNAP)
LOC_EXGLIBADV = local_snap/$(GLIBADV)
LOC_EXSNAPADV = local_snap/$(SNAPADV)
LOC_EXSNAPEXP = local_snap/$(SNAPEXP)

#	Local g++ include directories (LIBS variable defined inside snap's Makefile.config)
LIBS += -I$(LOC_EXGLIB)  -I$(LOC_EXSNAP) -I$(LOC_EXSNAPADV) 
#-I$(LOC_EXGLIBADV) -I$(LOC_EXSNAPEXP) 

## Main application file
MAIN = hysgen_main
DEPH = $(LOC_EXGLIB)/loc_bd.h $(LOC_EXSNAP)/loc_alg.h $(LOC_EXSNAP)/loc_gbase.h $(LOC_EXSNAP)/loc_graph.h $(LOC_EXSNAP)/loc_subgraph.h $(LOC_EXSNAPADV)/hysgen.h
DEPCPP = $(LOC_EXSNAP)/loc_gbase.cpp $(LOC_EXSNAP)/loc_graph.cpp $(LOC_EXSNAP)/loc_subgraph.cpp $(LOC_EXSNAPADV)/hysgen.cpp
CXXFLAGS += $(CXXOPENMP)
#CXXFLAGS += -g -rdynamic
#CXXFLAGS += -ggdb
#CXXFLAGS += -ggdb3 -rdynamic

