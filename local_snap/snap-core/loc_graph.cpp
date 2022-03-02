#include "Snap.h"
#include "loc_Snap.h"

/////////////////////////////////////////////////
// Bahman: Hyper Graph
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

// // Add a node of ID NId to the graph and create edges to all nodes in vector NbrNIdV. *****/////*****Bahman: Hold for now
// int THGraph::AddNode(const int& NId, const TIntV& NbrNIdV) {
// int NewNId;
// if (NId == -1) {
// NewNId = MxNId;  MxNId++;
// } else {
// IAssertR(! IsNode(NId), TStr::Fmt("NodeId %d already exists", NId));
// NewNId = NId;
// MxNId = TMath::Mx(NewNId+1, MxNId());
// }
// TNode& Node = NodeH.AddDat(NewNId);
// Node.Id = NewNId;
// Node.NIdV = NbrNIdV;
// Node.NIdV.Sort();
// NEdges += Node.GetDeg();
// for (int i = 0; i < NbrNIdV.Len(); i++) {
// GetNode(NbrNIdV[i]).NIdV.AddSorted(NewNId);
// }
// return NewNId;
// }

// // Add a node of ID NId to the graph and create edges to all nodes in the vector NIdVId in the vector pool Pool). *****/////*****Bahman: Hold for now
// int THGraph::AddNode(const int& NId, const TVecPool<TInt>& Pool, const int& NIdVId) {
// int NewNId;
// if (NId == -1) {
// NewNId = MxNId;  MxNId++;
// } else {
// IAssertR(!IsNode(NId), TStr::Fmt("NodeId %d already exists", NId));
// NewNId = NId;
// MxNId = TMath::Mx(NewNId+1, MxNId());
// }
// TNode& Node = NodeH.AddDat(NewNId);
// Node.Id = NewNId;
// Node.NIdV.GenExt(Pool.GetValVPt(NIdVId), Pool.GetVLen(NIdVId));
// Node.NIdV.Sort();
// NEdges += Node.GetDeg();
// return NewNId;
// }


// Delete node of ID NId from the graph. *****/////*****Bahman: Hold for now
void THGraph::DelNode(const int& NId) {
  AssertR(IsNode(NId), TStr::Fmt("NodeId %d does not exist", NId));
  TNode& Node = GetNode(NId);
  int EId = Node.EIdV.GetVal(0);
  int NumNei;
  TIntSet& ENIdsHS = GetEdge(EId).NeiNIdSH;
//  THash<TInt, TNode>& ENIdsHS = GetEdge(EId).NeiNIdSH;
  for (int e = 0; e < Node.GetDeg(); e++) {
    // Remove node from neighboring nodes' neighbor list:
    EId = Node.EIdV.GetVal(e);
    ENIdsHS = GetEdge(EId).NeiNIdSH;
    int iKey;
    for (int i=0; i < ENIdsHS.Len(); i++) {
//    for (int i=ENIdsHS.FFirstKeyId(); ENIdsHS.FNextKeyId(i); ) {
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
    } else { DelEdge(EId); } ////// DelEdge must be implemented. A lot of conflicts must be handled :(
    // Remove node from HGraph (list of nodes):
    NodeH.DelKey(NId);
  }
  delete &Node;
}


//inline THGraph::TEdge::TEdge(const int& EId, const TIntV& NodesIDsV, const THGraph* GraphPt) {
//  Id = EId;
//  for (int i = 0; i < NodesIDsV.Len(); i++){
//    NeiNIdSH.AddDat(NodesIDsV[i], GraphPt->GetNode(NodesIDsV[i]));
//  }
//  Graph = GraphPt;
//}

/*//  void THGraph::UpdENodes(const int& EId, TIntSet ENodesHS){
//    for (int i = 0; i < ENodesHS.Len(); i++) {
//      if (! GetNode(ENodesHS.GetKey(i)).EIdV.IsIn(EId)) {
//        GetNode(ENodesHS.GetKey(i)).EIdV.Add(EId);
//        for (int j = 0; j < ENodesHS.Len(); j++) {
//          if (i == j) { continue; }
//          if (! GetNode(ENodesHS.GetKey(i)).NbrNIdENumH.IsKey(
//          ENodesHS.GetKey(j))) {
//          GetNode(ENodesHS.GetKey(i)).NbrIdV.Add(ENodesHS.GetKey(j));
//          GetNode(ENodesHS.GetKey(i)).NbrNIdENumH.AddDat(ENodesHS.GetKey(j), 1);
//          } else {
//            GetNode(ENodesHS.GetKey(i)).NbrNIdENumH.AddDat
//              (ENodesHS.GetKey(j),
//               1 + ENodesHS.GetKey(i)).NbrNIdENumH.GetDat(ENodesHS.GetKey(j));
//          }
//        }
//      }
//    }
//  }*/

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

//int THGraph::AddEdge(const TIntSet& NIdH, int& EId){
//  THash<TInt, TNode> NodeIdsHS;
//  NodeIdsHS(NIdH.Len());
//  for (int i=0; i < NIdH.Len(); i++) {
//    if (!IsNode(NIdH[i])) { AddNode(NIdH[i]); }
//    NodeIdsHS.AddDat(NIdH[i], GetNode(NIdH[i]));
//  }
//  return AddEdge(NodeIdsHS, EId);
//}


bool THGraph::IsEdge(const TIntSet& NIdH) {
  if (NIdH.Len() < 2) { return false; }
  if (! IsNode(NIdH.GetKey(0))) { return false; }
  int NId = NIdH.GetKey(0);
  TIntSet SharedEIdsH;//(NodeH.GetDat(NId).GetEdges());
  NodeH.GetDat(NId).GetEIDs(SharedEIdsH);
  
//  printf("\n------ Before:\t"); //todo debug
//  for (THashSetKeyI<TInt> e_b = SharedEIdsH.BegI(); e_b < SharedEIdsH.EndI(); e_b++) {printf("%d,\t",e_b.GetKey());}//todo debug
//  printf("\n"); //todo debug
  for (int n = 1; n< NIdH.Len(); n++){
    int N1Id = NIdH.GetKey(n);
    TIntSet NeiEIdH;//(NodeH.GetDat(NId).GetEdges());
    NodeH.GetDat(N1Id).GetEIDs(NeiEIdH);
    TIntersect(SharedEIdsH, NeiEIdH);
    if (SharedEIdsH.Len()==0) { return false; }
  }
//  printf("\n++++++ After:\t"); //todo debug
//  for (THashSetKeyI<TInt> e_a = SharedEIdsH.BegI(); e_a < SharedEIdsH.EndI(); e_a++) {printf("%d,\t",e_a.GetKey());}//todo debug
//  printf("\n"); //todo debug
  
  for (THashSetKeyI<TInt> e = SharedEIdsH.BegI(); e < SharedEIdsH.EndI(); e++) {
    if(GetEI(e.GetKey()).Len() == NIdH.Len()) {
//      {//todo debug
//        TIntV DbgEdgeNIds;
//        GetEI(e.GetKey()).GetNodesV(DbgEdgeNIds);
//        printf("\nEdge ID #%d / %d:\t", e.GetKey(), SharedEIdsH.Len());
//        for (int n_e = 0; n_e < DbgEdgeNIds.Len(); n_e++) {printf("%d,\t",DbgEdgeNIds[n_e]);}
//        printf("\nNew Edge:\t");
//        for (THashSetKeyI<TInt> n_n = NIdH.BegI(); n_n < NIdH.EndI(); n_n++) {printf("%d,\t",n_n.GetKey());}
//        printf("\n");
//      }
      return true;
    }
  }
  return false;
}

/// Add an edge between the nodes in NodeIdsHS set.
int THGraph::AddEdge(const TIntSet& NIdH, int& EId) {
//  /// Verifying edge size
//  if (NIdH.Len() < 2) { AddNode(NIdH.GetKey(0)); return -1; }
//  /// Verifying if it is not a duplicate edge
//  ///TODO: Think about maybe incorporating the repeated edges, maybe in terms of a weight or whatever.
//  ///TODO: generative approach for above:
//  /// - membership generates hyperedge.
//  /// - then again the membership generates k of them. s.t like that to escape
//  /// from falling into considering nonexisting cases of different
//  /// repitition/values. maybe by modeling as some weight.
//  bool FlgNewEdge = false;
//  if (! IsNode(NIdH.GetKey(0))) { FlgNewEdge = true; }
//  else {
//    int NId0 = NIdH.GetKey(0), NIdI;
//    for (int i = 1; i< NIdH.Len(); i++){
//      NIdI = NIdH.GetKey(i);
//      if (! NodeH.GetDat(NId0).HasNeiN(NIdI)) {
//        FlgNewEdge = true;
//        break;
//      }
//    }
//  }
//  if (! FlgNewEdge) {
//    return -1;
//  }
  if (IsEdge(NIdH)) { return -1; }
  /// Adding the new edge
//  for (int i=0; i < NIdH.Len(); i++) {
  EId = TMath::Mx(EId, MxEId());
  MxEId = EId + 1;
  IAssertR(!IsEdgeId(EId), TStr::Fmt("EdgeId %d already exists", EId));
  EdgeH.AddDat(EId, TEdge(EId, NIdH, this));
  EdgeH.GetDat(EId).UpdNEInfo(NIdH);
//  }
  NEdges++;
  N2Edges += (NIdH.Len() * (NIdH.Len()-1))/2;
  return EId;
}

//int THGraph::AddEdge(const TIntV& NodeV) {
//  THash<TInt, TNode> NodeIdsHS;
//  NodeIdsHS(NodeV.Len());
//  for (int i=0; i < NodeV.Len(); i++) {
//    NodeIdsHS.AddDat(NodeV[i], GetNode(NodeV[i]));
//  }
//  return THGraph::AddEdge(NodeIdsHS);
//}

int THGraph::AssertNodes(const TIntSet& NodesIS) {
  int NKey;
  for (int N=0; N < NodesIS.Len(); N++) {
//  for (int N=NodesIS.FFirstKeyId(); NodesIS.FNextKeyId(N); ) {
    NKey = NodesIS.GetKey(N);
    if (! IsNode(NKey)) {return NKey;}
  }
  return -1; // No edge id
}

// // Add an edge between SrcNId and DstNId to the graph.
// int THGraph::AddEdgeUnchecked(const int& SrcNId, const int& DstNId) {
// GetNode(SrcNId).NIdV.Add(DstNId);
// if (SrcNId!=DstNId) { // not a self edge
// GetNode(DstNId).NIdV.Add(SrcNId); }
// NEdges++;
// return -1; // No edge id
// }

// // Add an edge between SrcNId and DstNId to the graph and create the nodes if they don't yet exist.
// int THGraph::AddEdge2(const int& SrcNId, const int& DstNId) {
// if (! IsNode(SrcNId)) { AddNode(SrcNId); }
// if (! IsNode(DstNId)) { AddNode(DstNId); }
// if (GetNode(SrcNId).IsNbrNId(DstNId)) { return -2; } // edge already exists
// GetNode(SrcNId).NIdV.AddSorted(DstNId);
// if (SrcNId!=DstNId) { // not a self edge
// GetNode(DstNId).NIdV.AddSorted(SrcNId); }
// NEdges++;
// return -1; // No edge id
// }

void THGraph::DelEdge(const int& EId) {
  IAssertR(IsEdgeId(EId), TStr::Fmt("EdgeId %d not found", EId));
  TIntSet NodeIdsHS = GetEdge(EId).NeiNIdSH;
//  THash<TInt, TNode> NodeIdsHS = GetEdge(EId).NeiNIdSH;
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

//bool THGraph::IsEdge(const THash<TInt, TNode>& EdgeIdS) const {
//  TIntV EIdsV;
//  int N = EdgeIdS.GetKey(EdgeIdS.FFirstKeyId());
//  int Nnext = EdgeIdS.GetKey(EdgeIdS.FNextKeyId(N));
//  if (GetNode(N).HasNeiN(Nnext)) { //Function to be defined
//    if (EdgeIdS.Len()<3) {return true;}
//    else {
//      // GetEdges must returns the Edge ID of the edges that contain node N.
//      GetNode(N).GetEdges(EIdsV);
//      for (int i=0; i<EIdsV.Len(); i++){ //Iterate over Edge IDs
//        if (EdgeIdS == GetEdge(EIdsV[i]).NeiNodeH) {return true;}
//      }
//    }
//  }
//}

//// Return an iterator referring to edge (SrcNId, DstNId) in the graph.
//THGraph::TEdgeI THGraph::GetEI(const int& EId) const {
//  const int MnNId = TMath::Mn(SrcNId, DstNId);
//  const int MxNId = TMath::Mx(SrcNId, DstNId);
//  const TNodeI SrcNI = GetNI(MnNId);
//  const int NodeN = SrcNI.NodeHI.GetDat().NbrIdV.SearchBin(MxNId);
//  IAssert(NodeN != -1);
//  return TEdgeI(SrcNI, EndNI(), NodeN);
//}

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

//// Check the graph data structure for internal consistency.
//bool THGraph::IsOk(const bool& ThrowExcept) const {
//  bool RetVal = true;
//  for (int N = NodeH.FFirstKeyId(); NodeH.FNextKeyId(N); ) {
//    const TNode& Node = NodeH[N];
//    if (! Node.NIdV.IsSorted()) {
//      const TStr Msg = TStr::Fmt("Neighbor list of node %d is not sorted.", Node.GetId());
//      if (ThrowExcept) { EAssertR(false, Msg); } else { ErrNotify(Msg.CStr()); }
//      RetVal=false;
//    }
//    int prevNId = -1;
//    for (int e = 0; e < Node.GetDeg(); e++) {
//      if (! IsNode(Node.GetNbrNId(e))) {
//        const TStr Msg = TStr::Fmt("Edge %d --> %d: node %d does not exist.",
//                                   Node.GetId(), Node.GetNbrNId(e), Node.GetNbrNId(e));
//        if (ThrowExcept) { EAssertR(false, Msg); } else { ErrNotify(Msg.CStr()); }
//        RetVal=false;
//      }
//      if (e > 0 && prevNId == Node.GetNbrNId(e)) {
//        const TStr Msg = TStr::Fmt("Node %d has duplicate edge %d --> %d.",
//                                   Node.GetId(), Node.GetId(), Node.GetNbrNId(e));
//        if (ThrowExcept) { EAssertR(false, Msg); } else { ErrNotify(Msg.CStr()); }
//        RetVal=false;
//      }
//      prevNId = Node.GetNbrNId(e);
//    }
//  }
//  int EdgeCnt = 0;
//  for (TEdgeI EI = BegEI(); EI < EndEI(); EI++) { EdgeCnt++; }
//  if (EdgeCnt != GetEdges()) {
//    const TStr Msg = TStr::Fmt("Number of edges counter is corrupted: GetEdges():%d, EdgeCount:%d.", GetEdges(), EdgeCnt);
//    if (ThrowExcept) { EAssertR(false, Msg); } else { ErrNotify(Msg.CStr()); }
//    RetVal=false;
//  }
//  return RetVal;
//}

// Print the graph in a human readable form to an output stream OutF.
void THGraph::Dump(FILE *OutF) const {
  const int NodePlaces = (int) ceil(log10((double) GetNodes()));
  fprintf(OutF, "-------------------------------------------------\nUndirected Node Graph: nodes: %d, edges: %d\n", GetNodes(), GetEdges());
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

// Return a small graph on 5 nodes and 5 edges.
PHGraph THGraph::GetSmallGraph() {
  PHGraph Graph = THGraph::New();
  for (int i = 0; i < 5; i++) { Graph->AddNode(i); }
  TIntV EV(3);
//  EV([0,1]);
  int el[] = {0,1,3};
  for (int i = 0; i < 3; i++)
    EV.Add(el[i]);
  Graph->AddEdge(EV);
//  Graph->AddEdge(0,1);  Graph->AddEdge(0,2);
//  Graph->AddEdge(0,3);  Graph->AddEdge(0,4);
//  Graph->AddEdge(1,2);
  return Graph;
}
