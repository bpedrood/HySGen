#include "Snap.h"
#include "loc_Snap.h"

namespace TSnap {

/////////////////////////////////////////////////
// Graph Algorithms
  
PHGraph GetSubGraph(const TPt<THGraph>& Graph, const TIntV& NIdV, const double& RenumberNodes) {
  PHGraph NewGraphPt = THGraph::New();
  THGraph& NewGraph = *NewGraphPt;
  NewGraph.Reserve(NIdV.Len(), -1);
  for (int n = 0; n < NIdV.Len(); n++) {
    if (Graph->IsNode(NIdV[n])) {
      NewGraph.AddNode(Graph->GetNI(NIdV[n]));
    }
  }
  bool EdgeAddable;
  TIntV NeiNIdV;
  for (THGraph::TEdgeI EI = Graph->BegEI(); EI < Graph->EndEI(); EI++) {
    EdgeAddable = true;
    EI.GetEdge().GetNbrNodes(NeiNIdV);
    for (int e=0; e < NeiNIdV.Len(); e++) {
      if (!NewGraph.IsNode(NeiNIdV[e])) { EdgeAddable = false; break; }
//    for (THGraph::TNodeI NeiNI = EI.BegNI(); NeiNI < EI.EndNI(); NeiNI++) {
//      if (!NewGraph.IsNode(NeiNI.GetId())) { EdgeAddable = false; break; }
    }
    if (EdgeAddable) { NewGraph.AddEdge(EI); }
  }
  NewGraph.Defrag();
  return NewGraphPt;
}

} // namespace TSnap
