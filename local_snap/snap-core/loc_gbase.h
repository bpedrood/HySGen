/*
//#//////////////////////////////////////////////
/// Graph Flags, used for quick testing of graph types. ##TGraphFlag
//delete TGraphFlag;
typedef enum TGraphFlag_ {
  gfUndef=0,    ///< default value, no flags
  gfDirected,   ///< directed graph (TNGraph, TNEGraph), else graph is undirected TUNGraph
  gfMultiGraph, ///< have explicit edges (multigraph): TNEGraph, TNodeEdgeNet
  gfNodeDat,    ///< network with data on nodes
  gfEdgeDat,    ///< network with data on edges
  gfSources,    ///< nodes only store out-edges (but not in-edges). See TBigNet
  gfBipart,     ///< bipartite graph
  gfHyperGraph,	///< Bahman: Hypergraph: THGraph
  gfMx          ///< sentinel, last value for iteration
} TGraphFlag;

///// Types for tables, sparse and dense attributes.
//typedef enum TAttrType_ {atInt, atFlt, atStr} TAttrType;

namespace TSnap {

//todo Uncomment upon error in the below object and then figure out a way to tackle the duplicity
///// Tests (at compile time) if the graph is directed.
//template <class TGraph> struct IsDirected   { enum { Val = 0 }; };
///// Tests (at compile time) if the graph is a multigraph with multiple edges between the same nodes.
//template <class TGraph> struct IsMultiGraph { enum { Val = 0 }; };
///// Tests (at compile time) if the graph is a network with data on nodes.
//template <class TGraph> struct IsNodeDat    { enum { Val = 0 }; };
///// Tests (at compile time) if the graph is a network with data on edges.
//template <class TGraph> struct IsEdgeDat    { enum { Val = 0 }; };
///// Tests (at compile time) if the nodes store only out-edges, but not in-edges.
//template <class TGraph> struct IsSources    { enum { Val = 0 }; };
///// Tests (at compile time) if the graph is a bipartite graph type.
//template <class TGraph> struct IsBipart     { enum { Val = 0 }; };

/// Tests (at compile time) if the graph is a HYPERGRAPH. E.g., is used in subgraph.h:113
  template <class TGraph> struct IsHyperGraph { enum { Val = 0 }; }; //Is set to 1 in graph.h:1227. How really that line is get called? No idea yet!

/// For quick testing of the properties of the graph/network object (see TGraphFlag).
#ifdef HasGraphFlag
#undef HasGraphFlag
#define HasGraphFlag(TGraph, Flag) \
  ((Flag)==gfDirected ? TSnap::IsDirected<TGraph::TNet>::Val : \
  (Flag)==gfMultiGraph ? TSnap::IsMultiGraph<TGraph::TNet>::Val : \
  (Flag)==gfHyperGraph ? TSnap::IsHyperGraph<TGraph::TNet>::Val : \
  (Flag)==gfNodeDat ? TSnap::IsNodeDat<TGraph::TNet>::Val : \
  (Flag)==gfEdgeDat ? TSnap::IsEdgeDat<TGraph::TNet>::Val : \
  (Flag)==gfSources ? TSnap::IsSources<TGraph::TNet>::Val : \
  (Flag)==gfBipart ? TSnap::IsBipart<TGraph::TNet>::Val : 0)
//#define HasGraphFlag(TGraph, Flag) \
//  (\
//    (Flag)==gfHyperGraph ? TSnap::IsHyperGraph<TGraph::TNet>::Val : 0\
//  )
#endif
//#define HasGraphFlag(TGraph, Flag) \
//  (\
//    (Flag)==gfDirected ? TSnap::IsDirected<TGraph::TNet>::Val : \
//    (Flag)==gfMultiGraph ? TSnap::IsMultiGraph<TGraph::TNet>::Val : \
//    (Flag)==gfNodeDat ? TSnap::IsNodeDat<TGraph::TNet>::Val : \
//    (Flag)==gfEdgeDat ? TSnap::IsEdgeDat<TGraph::TNet>::Val : \
//    (Flag)==gfSources ? TSnap::IsSources<TGraph::TNet>::Val : \
//    (Flag)==gfBipart ? TSnap::IsBipart<TGraph::TNet>::Val : 0 \
//    (Flag)==gfHyperGraph ? TSnap::IsHyperGraph<TGraph::TNet>::Val : \
//  )


}  // namespace TSnap
*/
