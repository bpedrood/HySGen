/*
#include "Snap.h"
#include "loc_Snap.h"

/////////////////////////////////////////////////
// Graph Base
namespace TSnap {

TStr GetFlagStr(const TGraphFlag& GraphFlag) {
  switch (GraphFlag) {
    case gfUndef : return "Undef";
    case gfDirected : return "Directed";
    case gfMultiGraph : return "Multigraph";
    case gfHyperGraph : return "Hypergraph";
    case gfNodeDat : return "NodeDat";
    case gfEdgeDat : return "EdgeDat";
    case gfSources : return "Sources";
    case gfBipart : return "Bipartite";
    default: FailR("Unknown graph type");
  };
  return TStr();
}

};  // namespace TSnap
*/
