#include "Snap.h"
#include "loc_Snap.h"

/////////////////////////////////////////////////
bool THGraph::HasFlag(const TGraphFlag& Flag) const {
  return HasGraphFlag(THGraph::TNet, Flag);
}

// Add a node of ID NId to the graph.
int THGraph::AddNode(int NId, TStr NName) {
  if (NId == -1) {
    NId = MxNId;  MxNId++;
  } else {
    IAssertR(!IsNode(NId), TStr::Fmt("NodeId %d already exists", NId));
    MxNId = TMath::Mx(NId+1, MxNId());
  }
  NodeH.AddDat(NId, TNode(NId, NName));
  return NId;
}

// Add a node of ID NId to the graph.
int THGraph::AddNodeUnchecked(int NId, TStr NName) {
  if (IsNode(NId)) { return -1;}
  MxNId = TMath::Mx(NId+1, MxNId());
  NodeH.AddDat(NId, TNode(NId, NName));
  return NId;
}

// Delete node of ID NId from the graph.
void THGraph::DelNode(const int& NId) {
  AssertR(IsNode(NId), TStr::Fmt("NodeId %d does not exist", NId));
  TNode& Node = GetNode(NId);
  int EId = Node.EIdV.GetVal(0);
  int NumNei;
  TIntSet& ENIdsHS = GetEdge(EId).NeiNIdSH;
  for (int e = 0; e < Node.GetDeg(); e++) {
    // Remove node from neighboring nodes' neighbor list:
    EId = Node.EIdV.GetVal(e);
    ENIdsHS = GetEdge(EId).NeiNIdSH;
    int iKey;
    for (int i=0; i < ENIdsHS.Len(); i++) {
      iKey = ENIdsHS.GetKey(i);
      if (iKey == NId) { continue; }
      GetNode(iKey).DelNeighbor(NId);
    }
    /// Adjust N2Edges
    NumNei = ENIdsHS.Len();
    N2Edges = N2Edges - (NumNei*(NumNei-1))/2 + ((NumNei-1)*(NumNei-2))/2;
    // Remove NId from contributing edges/Delete edges with NId with size of 2:
    if (ENIdsHS.Len() > 2) {
      GetEdge(EId).NeiNIdSH.DelKey(NId);
    } else { DelEdge(EId); } 
    // Remove node from HGraph (list of nodes):
    NodeH.DelKey(NId);
  }
  delete &Node;
}

void THGraph::TNode::UpdEInfo(const int& EId, const TIntSet& ENodesHS){
  if (! EIdV.IsIn(EId)) {
    EIdV.Add(EId);
    for (int j = 0; j < ENodesHS.Len(); j++) {
      if (Id == ENodesHS.GetKey(j)) { continue; }
      if (! NbrNIdENumH.IsKey(ENodesHS.GetKey(j))) {
        NbrNIdENumH.AddDat(ENodesHS.GetKey(j), 1);
      } else {
        NbrNIdENumH.AddDat(ENodesHS.GetKey(j), 1 + NbrNIdENumH.GetDat(ENodesHS.GetKey(j)));
      }
    }
  }
}

inline THGraph::TEdge::TEdge( THGraph* GraphPt, const TIntSet& NodeIdsHS) {
  TInt EIdCandidate = GraphPt->NEdges;
  while (GraphPt->EdgeH.IsKey(EIdCandidate)) { EIdCandidate++; }
  TEdge(EIdCandidate, NodeIdsHS, GraphPt);
}

inline THGraph::TEdge::TEdge( THGraph* GraphPt, const TIntV& NodeIdsV) {
  TInt EIdCandidate = GraphPt->NEdges;
  while (GraphPt->EdgeH.IsKey(EIdCandidate)) { EIdCandidate++; }
  TEdge(EIdCandidate, NodeIdsV, GraphPt);
}

void THGraph::TEdge::UpdNEInfo(const TIntSet& ENodesHS){
  for (int i = 0; i < ENodesHS.Len(); i++) {
    Graph->GetNode(ENodesHS.GetKey(i)).UpdEInfo(Id, ENodesHS);
  }
}

bool THGraph::IsEdge(const TIntSet& NIdH) {
  if (NIdH.Len() < 2) { return false; }
  if (! IsNode(NIdH.GetKey(0))) { return false; }
  int NId = NIdH.GetKey(0);
  TIntSet SharedEIdsH;
  NodeH.GetDat(NId).GetEIDs(SharedEIdsH);
  for (int n = 1; n< NIdH.Len(); n++){
    int N1Id = NIdH.GetKey(n);
    TIntSet NeiEIdH;
    NodeH.GetDat(N1Id).GetEIDs(NeiEIdH);
    TIntersect(SharedEIdsH, NeiEIdH);
    if (SharedEIdsH.Len()==0) { return false; }
  }
  for (THashSetKeyI<TInt> e = SharedEIdsH.BegI(); e < SharedEIdsH.EndI(); e++) {
    if(GetEI(e.GetKey()).Len() == NIdH.Len()) {
      return true;
    }
  }
  return false;
}

/// Add an edge between the nodes in NodeIdsHS set.
int THGraph::AddEdge(const TIntSet& NIdH, int& EId) {
  if (IsEdge(NIdH)) { return -1; }
  EId = TMath::Mx(EId, MxEId());
  MxEId = EId + 1;
  IAssertR(!IsEdgeId(EId), TStr::Fmt("EdgeId %d already exists", EId));
  EdgeH.AddDat(EId, TEdge(EId, NIdH, this));
  EdgeH.GetDat(EId).UpdNEInfo(NIdH);
  NEdges++;
  N2Edges += (NIdH.Len() * (NIdH.Len()-1))/2;
  return EId;
}

int THGraph::AssertNodes(const TIntSet& NodesIS) {
  int NKey;
  for (int N=0; N < NodesIS.Len(); N++) {
    NKey = NodesIS.GetKey(N);
    if (! IsNode(NKey)) {return NKey;}
  }
  return -1;
}

void THGraph::DelEdge(const int& EId) {
  IAssertR(IsEdgeId(EId), TStr::Fmt("EdgeId %d not found", EId));
  TIntSet NodeIdsHS = GetEdge(EId).NeiNIdSH;
  int iKey, jKey;
  for (int i=NodeIdsHS.FFirstKeyId(); NodeIdsHS.FNextKeyId(i); ) {
    iKey = NodeIdsHS.GetKey(i);
    for (int j=NodeIdsHS.FFirstKeyId(); NodeIdsHS.FNextKeyId(j); ) {
      if (i==j) { continue; }
      jKey = NodeIdsHS.GetKey(j);
      GetNode(iKey).DelNeighbor(jKey);
    }
  }
  EdgeH.DelKey(EId);
  int ESize = GetEdge(EId).NeiNIdSH.Len();
  NEdges--;
  N2Edges -= (ESize*(ESize-1))/2;
  delete &GetEdge(EId);
}

// Get a vector IDs of all nodes in the graph.
void THGraph::GetNIdV(TIntV& NIdV) const {
  NIdV.Gen(GetNodes(), 0);
  for (int N=NodeH.FFirstKeyId(); NodeH.FNextKeyId(N); ) {
    NIdV.Add(NodeH.GetKey(N)); }
}

// Defragment the graph.
void THGraph::Defrag(const bool& OnlyNodeLinks) {
  int nKey;
  for (int n = EdgeH.FFirstKeyId(); EdgeH.FNextKeyId(n); ) {
    nKey = EdgeH.GetKey(n);
    if (! EdgeH.GetDat(nKey).NeiNIdSH.IsKeyIdEqKeyN()) {
      EdgeH.GetDat(nKey).NeiNIdSH.Defrag();
    }
  }
  for (int n = NodeH.FFirstKeyId(); NodeH.FNextKeyId(n); ) {
    nKey = NodeH.GetKey(n);
    NodeH.GetDat(nKey).EIdV.Pack();
    if (! NodeH.GetDat(nKey).NbrNIdENumH.IsKeyIdEqKeyN()) {
      NodeH.GetDat(nKey).NbrNIdENumH.Defrag();
    }
  }
  if (! OnlyNodeLinks) {
    if (! NodeH.IsKeyIdEqKeyN()){ NodeH.Defrag(); }
    if (! EdgeH.IsKeyIdEqKeyN()){ EdgeH.Defrag(); }
  }
}

void THGraph::PrintEdge(const int EId) {
  if (! EdgeH.IsKey(EId)) {
    printf("\nEdge Not Found!\n");
  }
  TIntV NV;
  TStr EStr("Edge nodes: ");
  EdgeH.GetDat(EId).GetNodesV(NV);
  for (int i = 0; i < NV.Len(); i++) {
    EStr += (" " + NV[i].GetStr());
  }
  EStr += "\n";
  printf(EStr.GetCStr());
}

// Print the graph in a human readable form to an output stream OutF.
void THGraph::Dump(FILE *OutF) const {
  const int NodePlaces = (int) ceil(log10((double) GetNodes()));
  int NKey;
  for (int N = NodeH.FFirstKeyId(); NodeH.FNextKeyId(N); ) {
    NKey = NodeH.GetKey(N);
    const TNode& Node = NodeH.GetDat(NKey);
    fprintf(OutF, "  %*d [%d] ", NodePlaces, Node.GetId(), Node.GetDeg());
    for (int edge = 0; edge < Node.GetDeg(); edge++) {
      fprintf(OutF, " %*d", NodePlaces, Node.GetNbrNId(edge)); }
    fprintf(OutF, "\n");
  }
  fprintf(OutF, "\n");
}

// Return a small hypergraph on 5 nodes and 5 edges.
PHGraph THGraph::GetSmallGraph() {
  PHGraph Graph = THGraph::New();
  for (int i = 0; i < 5; i++) { Graph->AddNode(i); }
  TIntV EV(3);
  int el[] = {0,1,3};
  for (int i = 0; i < 3; i++)
    EV.Add(el[i]);
  Graph->AddEdge(EV);
  return Graph;
}
