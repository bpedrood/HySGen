//todo I THINK THIS PIECE OF CODE IS JUST USELESS. REMOVE IF CODE COULD BE COMPILED SUCCESSFULLY
//namespace TSnap {
//namespace TSnapDetail {
//
//
//template <>
//struct TDelSelfEdges<PHGraph> { // HyperGraph specialization
//  static void Do(const PHGraph& Graph) {
//    TIntV EdgeV;
//    for (typename PHGraph::TObj::TEdgeI EI = Graph->BegEI(); EI < Graph->EndEI(); EI++) {
//      if (EI.Len() <= 1) { // Only contain one [or less!!] nodes means a self-loop
//        IAssert(EI.GetId() >= 0); // real edge id
//        EdgeV.Add(EI.GetId());
//      }
//    }
//    for (int i = 0; i < EdgeV.Len(); i++) { // delete all edges between a pair of nodes
//      Graph->DelEdge(EdgeV[i]);
//    }
//  }
//};
//
//
//} // namespace TSnapDetail
//}; // namespace TSnap
