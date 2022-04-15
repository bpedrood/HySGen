//#//////////////////////////////////////////////
/// Undirected graphs

/// [Undirected] Hyper-graphs
class THGraph;
/// Pointer to a hypergraph graph (THGraph)
typedef TPt<THGraph> PHGraph;

//#//////////////////////////////////////////////
/// Hypergraph.
class THGraph {
public:
  typedef THGraph TNet;
  typedef TPt<THGraph> PNet;
public:
  class TNode {
  private:
    TInt Id;
    TStr Name;
    TIntV EIdV; // Vector of Edge IDs
    THash<TInt, TInt> NbrNIdENumH; 
  public:
    TNode() : Id(-1), Name(""), EIdV(), NbrNIdENumH() { }
    TNode(const int& NId) : Id(NId), Name(""), EIdV(), NbrNIdENumH() { }
    TNode(const int& NId, const TStr& NName) : Id(NId), Name(NName), EIdV(), NbrNIdENumH() { }
    TNode(const TNode& Node) : Id(Node.Id), Name(Node.Name), EIdV(Node.EIdV), NbrNIdENumH(Node.NbrNIdENumH) { }
    TNode(TSIn& SIn) : Id(SIn), Name(SIn), EIdV(SIn), NbrNIdENumH(SIn) { }
    void Save(TSOut& SOut) const { Id.Save(SOut); Name.Save(SOut); EIdV.Save(SOut); NbrNIdENumH.Save(SOut);}
    int GetId() const { return Id; }
    TStr GetName() const { return Name; }
    int GetNbrNodes() const { return NbrNIdENumH.Len(); }
    void GetNbrNodes(TIntV& VNbrV) const { NbrNIdENumH.GetKeyV(VNbrV); }
    void GetNbrEdges(TIntV& ENbrV) const { ENbrV = EIdV; }
    int GetDeg() const { return EIdV.Len(); }
    int GetInDeg() const { return GetDeg(); }
    int GetOutDeg() const { return GetDeg(); }
    int GetNbrEId(const int& NE) const { return EIdV[NE]; } // Returns N-th neighboring edge.
    int GetNbrNId(const int& NN) const { return NbrNIdENumH.GetKey(NN); } // Returns N-th neighboring node.
    bool IsNbrNId(const int& NId) const { return HasNeiN(NId); }
    bool IsInNId(const int& NId) const { return HasNeiN(NId); }
    bool IsOutNId(const int& NId) const { return HasNeiN(NId); }
    void SortNIdV() { NbrNIdENumH.SortByKey();}
    bool HasNeiN(const TNode &Node) const { return NbrNIdENumH.IsKey(Node.GetId()); } //Checks for a neighboring node Node.
    bool HasNeiN(const int& NId) const { return NbrNIdENumH.IsKey(NId); } //Checks for a neighboring node Node.
    void UpdEInfo(const int& EId, const TIntSet& ENodesHS);
    void GetEIDs(TIntV& NeiEIdV) const { NeiEIdV=EIdV; } // Returns the edges containing the Node.
    void GetEIDs(TIntSet& NeiEIdH) const {TIntSet EH(EIdV); NeiEIdH=EH; } // Returns the edges containing the Node.
    int GetEdges() const { return GetDeg(); } // Returns the number of edges containing the Node.
    int Get2Edges() const { return (NbrNIdENumH.Len() * (NbrNIdENumH.Len()-1))/2; } // Returns the number of edges containing the Node.
    void AddNeighbor(const int& NId) {
      if (NbrNIdENumH.IsKey(NId)) {
        NbrNIdENumH.AddDat(NId, NbrNIdENumH.GetDat(NId)+1);
      }else {
        NbrNIdENumH.AddDat(NId,1);
      }
    }
    int DelNeighbor(const int& NId) {
      if (NbrNIdENumH.IsKey(NId)) {
        if (NbrNIdENumH.GetDat(NId)>1) {
          NbrNIdENumH.AddDat(NId, NbrNIdENumH.GetDat(NId)-1);
        } else {
          NbrNIdENumH.DelKey(NId);
        }
      }
      return NbrNIdENumH.GetDat(NId);
    }
    bool operator == (const TNode& Node) const {
      TIntV* ThisNNeiV;
      TIntV* InpNNeiV;
      GetNbrNodes(*ThisNNeiV); Node.GetNbrNodes(*InpNNeiV);
      if (*ThisNNeiV == *InpNNeiV && EIdV == Node.EIdV) {
        IAssertR(Id != Node.Id, "All node's neighbors are the same, but IDs don't match.");
        IAssertR(Name.EqI(Node.Name), "All node's neighbors are the same, but Names don't match.");
        return true;
      }
      return false;
    }
//    ~TNode() {}
    friend class THGraph;
  };
  /// Node Iterator
  class TNodeI {
  private:
    typedef THash<TInt, TNode>::TIter THashIter;
    THashIter NodeHI;
  public:
    TNodeI() : NodeHI() { }
    TNodeI(const THashIter& NodeHIter) : NodeHI(NodeHIter) { }
    TNodeI(const TNodeI& NodeI) : NodeHI(NodeI.NodeHI) { }
    TNodeI& operator = (const TNodeI& NodeI) { NodeHI = NodeI.NodeHI; return *this; }
    TNodeI& operator++ (int) { NodeHI++; return *this; }
    TNodeI& operator-- (int) { NodeHI--; return *this; }
    bool operator < (const TNodeI& NodeI) const { return NodeHI < NodeI.NodeHI; }
    bool operator == (const TNodeI& NodeI) const { return NodeHI == NodeI.NodeHI; }
    /// Returns ID of the current node.
    int GetId() const { return NodeHI.GetDat().GetId(); }
    TStr GetName() const { return NodeHI.GetDat().GetName(); }
    /// Returns number of neighboring nodes of the current node.
    int GetNbrNodes() const { return NodeHI.GetDat().GetNbrNodes(); }
    void GetNbrNodes(TIntV& VNbrV) const { NodeHI.GetDat().GetNbrNodes(VNbrV); }
    /// Returns the ID on neighboring edges.
    void GetEIDs(TIntV& NeiEIdV) { NodeHI.GetDat().GetEIDs(NeiEIdV); }
    void GetEIDs(TIntSet& NeiEIdH) { NodeHI.GetDat().GetEIDs(NeiEIdH); }
    /// Returns degree of the current node.
    int GetDeg() const { return NodeHI.GetDat().GetDeg(); }
    /// Returns in-degree of the current node (returns same as value GetDeg() since the graph is undirected).
    int GetInDeg() const { return NodeHI.GetDat().GetInDeg(); }
    /// Returns out-degree of the current node (returns same as value GetDeg() since the graph is undirected).
    int GetOutDeg() const { return NodeHI.GetDat().GetOutDeg(); }
    /// Returns ID of NodeN-th neighboring node.
    int GetNbrNId(const int& NodeN) const { return NodeHI.GetDat().GetNbrNId(NodeN); }
    /// Tests whether node with ID NId points to the current node.
    bool HasNeiInN(const int& NId) const { return NodeHI.GetDat().IsInNId(NId); }
    /// Tests whether the current node points to node with ID NId.
    bool HasNeiOutN(const int& NId) const { return NodeHI.GetDat().IsOutNId(NId); }
    /// Tests whether node with ID NId is a neighbor of the current node.
    bool HasNeiN(const int& NId) const { return NodeHI.GetDat().IsNbrNId(NId); }
    /// Returns number of [2]edges if being flatten to a regular graph.
    int Get2Edges() const { return NodeHI.GetDat().Get2Edges(); }
    /// Query a neighboring edge by Edge Id.
    bool HasEdge(const int& EId) {
      TIntV EV = NodeHI.GetDat().EIdV;
      for (int i = 0; i < EV.Len(); i++) {
        if (EId == EV[i]) { return true; }
        return false;
      }
    }
    friend class THGraph;
  };
  class TEdge {
  private:
    TInt Id;
    TIntSet NeiNIdSH;
    THGraph* Graph;
  public:
    TEdge() { TEdge(NULL); }
    TEdge(THGraph* GraphPt) : Id(-1), NeiNIdSH(), Graph(GraphPt) { }
    TEdge(const int& EId, const TIntSet& NodesIDsH, THGraph* GraphPt) : Id(EId), NeiNIdSH(NodesIDsH), Graph(GraphPt) { }
    TEdge(const int& EId, const TIntV& NodesIDsV, THGraph* GraphPt) : Id(EId), NeiNIdSH(NodesIDsV), Graph(GraphPt) { }
    TEdge(const TEdge& Edge) : Id(Edge.Id), NeiNIdSH(Edge.NeiNIdSH), Graph(Edge.Graph) { }
    TEdge(THGraph* GraphPt, const TIntSet& NodeIdsHS);
    TEdge(THGraph* GraphPt, const TIntV& NodeIdsV);
    TEdge(TSIn& SIn) : Id(SIn), NeiNIdSH(SIn), Graph(NULL) { }
    TEdge(TSIn& SIn, THGraph* GraphPt) : Id(SIn), NeiNIdSH(SIn), Graph(GraphPt) { }
    TEdge(const THashSet<TInt>& NIdH, const int& EId) {
      TIntV NIdV;
      NIdH.GetKeyV(NIdV);
      TEdge(EId, NIdV, Graph);
    }
    void Save(TSOut& SOut) const { Id.Save(SOut); NeiNIdSH.Save(SOut); }
    /// Gets edge ID.
    int GetId() const { return Id; }
    /// Gets the number of including nodes.
    int Len() const { return NeiNIdSH.Len(); }
    /// Gets the nodes included in the edge.
    void GetNodesV(TIntV& NeiV) const { NeiNIdSH.GetKeyV(NeiV); }
    /// Gets a vector of including nodes' Ids.
    void GetNbrNodes(TIntV& NIdsV) { NeiNIdSH.GetKeyV(NIdsV);}
    bool HasNode(const TNode& Node) const { return NeiNIdSH.IsKey(Node.GetId()); } // Checks the edge for node Node.
    void UpdNEInfo(const TIntSet& ENodesHS);
//    ~TEdge() {}
    friend class THGraph;
  };
  /// Edge iterator. Only forward iteration (operator++) is supported.
  class TEdgeI {
  private:
    typedef THash<TInt, TEdge>::TIter THashIter;
    THashIter EdgeHI;
  public:
    TEdgeI() : EdgeHI() { }
    TEdgeI(const THashIter& EdgeHIter) : EdgeHI(EdgeHIter) { }
    TEdgeI(const TEdgeI& EdgeI) : EdgeHI(EdgeI.EdgeHI) { }
    TEdgeI& operator = (const TEdgeI& EdgeI) { if (this!=&EdgeI) { EdgeHI=EdgeI.EdgeHI; }  return *this; }
    TEdgeI& operator ++ (int) { EdgeHI++; return *this; }
    bool operator < (const TEdgeI& EdgeI) const { return EdgeHI < EdgeI.EdgeHI; }
    bool operator == (const TEdgeI& EdgeI) const { return EdgeHI == EdgeI.EdgeHI; }
    /// Gets edge ID.
    int GetId() const { return EdgeHI.GetDat().GetId(); }
    TEdge GetEdge() {return EdgeHI.GetDat();}
    /// Gets the nodes included in the edge.
    void GetNodesV(TIntV& NeiV) const { EdgeHI.GetDat().GetNodesV(NeiV); }
    /// Get ID of Neighboring Nodes of the edge with id EId
    void GetNbrNodes(TIntV& NIDV) { EdgeHI.GetDat().GetNbrNodes(NIDV); }
    /// Gets the number of including nodes.
    int Len() const { return EdgeHI.GetDat().Len(); }
    int GetDeg() const { return Len(); }
    friend class THGraph;
  };
private:
  TCRef CRef;
  TInt MxNId, MxEId, NEdges, N2Edges;
  THash<TInt, TNode> NodeH; 
  THash<TInt, TEdge> EdgeH; 
private:
  TNode& GetNode(const int& NId) { return NodeH.GetDat(NId); }
  const TNode& GetNode(const int& NId) const { return NodeH.GetDat(NId); }
  TEdge& GetEdge(const int& EId) { return EdgeH.GetDat(EId); }
  const TEdge& GetEdge(const int& EId) const { return EdgeH.GetDat(EId); }
  int AssertNodes(const TIntSet& NodesIS);
  /// To be utilized for preventing addition of duplicate hyperedges.
  static void TIntersect(TIntSet& S, const TIntSet& S2) {
    /// S <- Intersect(S, S2)
    TIntV MarkDel(S.Len(),0); //todo maybe remove if inneffective (comparing to before)
    for (THashSetKeyI<TInt> i = S.BegI(); i < S.EndI(); i++) {
      if (!S2.IsKey(i.GetKey())) { MarkDel.Add(i.GetKey()); }
    }
    if (MarkDel.Len() > 0) { for (int j = 0; j < MarkDel.Len(); j++) S.DelKey(MarkDel[j]); }
  }
public:
  THGraph() : CRef(), MxNId(0), MxEId(0), NEdges(0), N2Edges(0), NodeH(), EdgeH() { }
  /// Constructor that reserves enough memory for a graph of Nodes nodes and Edges edges.
  explicit THGraph(const int& Nodes, const int& Edges) : MxNId(0), NEdges(0), N2Edges(0) { Reserve(Nodes, Edges); }
  THGraph(const THGraph& Graph) : MxNId(Graph.MxNId), MxEId(Graph.MxEId), NEdges(Graph.NEdges), N2Edges(Graph.N2Edges), NodeH(Graph.NodeH), EdgeH(Graph.EdgeH) { }
  /// Constructor that loads the graph from a (binary) stream SIn.
  THGraph(TSIn& SIn) : MxNId(SIn), MxEId(SIn), NEdges(SIn), N2Edges(SIn), NodeH(SIn), EdgeH(SIn) { }
  /// Saves the graph to a (binary) stream SOut.
  void Save(TSOut& SOut) const { MxNId.Save(SOut); MxEId.Save(SOut); NEdges.Save(SOut); N2Edges.Save(SOut); NodeH.Save(SOut); EdgeH.Save(SOut); }
  /// Static constructor that returns a pointer to the graph. Call: PUNGraph Graph = THGraph::New().
  static PHGraph New() { return new THGraph(); }
  /// Static constructor that returns a pointer to the graph and reserves enough memory for Nodes nodes and Edges edges. 
  static PHGraph New(const int& Nodes, const int& Edges) { return new THGraph(Nodes, Edges); }
  /// Static constructor that loads the graph from a stream SIn and returns a pointer to it.
  static PHGraph Load(TSIn& SIn) { return PHGraph(new THGraph(SIn)); }
  /// Allows for run-time checking the type of the graph (see the TGraphFlag for flags).
  bool HasFlag(const TGraphFlag& Flag) const;
  THGraph& operator = (const THGraph& Graph) {
    if (this!=&Graph) { MxNId=Graph.MxNId; MxEId=Graph.MxEId; NEdges=Graph.NEdges; N2Edges=Graph.N2Edges; NodeH=Graph.NodeH; EdgeH=Graph.EdgeH; } return *this; }
  /// Returns the number of nodes in the graph.
  int GetNodes() const { return NodeH.Len(); }
  /// Adds a node of ID NId to the network, noop if the node already exists. 
  int AddNodeUnchecked(int NId = -1, TStr NName = "");
  /// Adds a node of ID NId and Name NName to the graph. 
  int AddNode(int NId = -1, TStr NName = "");
  /// Adds a node of ID NodeI.GetId() to the graph.
  int AddNode(const TNodeI& NodeI) { return AddNode(NodeI.GetId(), NodeI.GetName()); }
  /// Deletes node of ID NId from the graph. 
  void DelNode(const int& NId);
  /// Deletes node of ID NodeI.GetId() from the graph.
  void DelNode(const TNode& NodeI) { DelNode(NodeI.GetId()); }
  /// Tests whether ID NId is a node.
  bool IsNode(const int& NId) const { return NodeH.IsKey(NId); }
  /// Returns an iterator referring to the first node in the graph.
  TNodeI BegNI() const { return TNodeI(NodeH.BegI()); }
  /// Returns an iterator referring to the past-the-end node in the graph.
  TNodeI EndNI() const { return TNodeI(NodeH.EndI()); }
  /// Returns an iterator referring to the node of ID NId in the graph.
  TNodeI GetNI(const int& NId) const { return TNodeI(NodeH.GetI(NId)); }
  /// Returns an ID that is larger than any node ID in the graph.
  int GetMxNId() const { return MxNId; }
  int GetVNbrNodes(int& NId) const { return GetNode(NId).GetNbrNodes(); }
  void GetVNbrEdges(int& NId, TIntV& EIDV) const { return GetNode(NId).GetNbrEdges(EIDV); }
  /// Returns the number of edges in the graph.
  int GetEdges() const { return EdgeH.Len(); }
  /// Returns the number of equivalent 2-edges in the graph.
  int Get2Edges() const { return N2Edges; }
  /// Get ID of Neighboring Nodes of the edge with id EId
  void GetENbrNodes(int EId, TIntV& NIDV) { GetEdge(EId).GetNbrNodes(NIDV); }
  /// Adds an edge between the nodes in NodeIdsHS set
  int AddEdge(const THash<TInt, TNode>& NodeIdsHS){
    TIntV NIdV;
    NodeIdsHS.GetKeyV(NIdV);
    return AddEdge(NIdV);
  }
  int AddEdge(const TIntSet& NIdH, int& EId);
  int AddEdge(const TIntSet& NIdH){ int _EId = -1; return AddEdge(NIdH, _EId); }
  int AddEdge(const TIntV& NodeV) { TIntSet NodeIdsHS(NodeV); return AddEdge(NodeIdsHS); }
  int AddEdge(const TEdgeI& EI){
    int _EId = EI.GetId();
    return AddEdge(EI.EdgeHI.GetDat().NeiNIdSH, _EId);
  }
  int AddEdge(const TEdgeI& EI, int& EId){ return AddEdge(EI.EdgeHI.GetDat().NeiNIdSH, EId); }
  /// Deletes the edge with id=EId.
  void DelEdge(const int& EId);
  void DelEdge(const TEdge& Edge) { DelEdge(Edge.Id); }
  /// Tests whether an edge with edge ID EId exists in the graph.
  bool IsEdgeId(const int& EId) const { return EdgeH.IsKey(EId); }
  /// Tests whether an edge exactly with the nodes specified in the input.
  bool IsEdge(const TIntSet& NIdH);
  /// Returns an iterator referring to the first edge in the graph.
  TEdgeI BegEI() const { return TEdgeI(EdgeH.BegI()); }
  /// Returns an iterator referring to the past-the-end edge in the graph.
  TEdgeI EndEI() const { return TEdgeI(EdgeH.EndI()); }
  /// Returns an iterator referring to edge with edge ID EId.
  TEdgeI GetEI(const int& EId) const { return TEdgeI(EdgeH.GetI(EId)); }
  /// Returns an ID of a random node in the graph.
  int GetRndNId(TRnd& Rnd=TInt::Rnd) { return NodeH.GetKey(NodeH.GetRndKeyId(Rnd, 0.8)); }
  /// Returns an interator referring to a random node in the graph.
  TNodeI GetRndNI(TRnd& Rnd=TInt::Rnd) { return GetNI(GetRndNId(Rnd)); }
  /// Returns an ID of a random edge in the graph.
  int GetRndEId(TRnd& Rnd=TInt::Rnd) { return EdgeH.GetKey(EdgeH.GetRndKeyId(Rnd, 0.8)); }
  /// Returns an interator referring to a random edge in the graph.
  TEdgeI GetRndEI(TRnd& Rnd=TInt::Rnd) { return GetEI(GetRndEId(Rnd)); }
  /// Returns a random edge in the graph.
  TEdge GetRndEdge(TRnd& Rnd=TInt::Rnd) { return GetEdge(GetRndEId(Rnd)); }
  /// Gets a vector IDs of all nodes in the graph.
  void GetNIdV(TIntV& NIdV) const;
  // Get a vector IDs of all edges in the graph.
  void GetEIdV(TIntV& EIdV) const { EdgeH.GetKeyV(EIdV); }
  /// Tests whether the graph is empty (has zero nodes).
  bool Empty() const { return GetNodes()==0; }
  /// Deletes all nodes and edges from the graph.
  void Clr() { MxNId=0; MxEId=0; NEdges=0; N2Edges=0; NodeH.Clr(); EdgeH.Clr(); }
  /// Reserves memory for a graph of Nodes nodes and Edges edges.
  void Reserve(const int& Nodes, const int& Edges) {
    if (Nodes>0) { NodeH.Gen(Nodes/2); } if (Edges>0) { EdgeH.Gen(Edges/2); }
  }
  /// Defragments the graph. 
  void Defrag(const bool& OnlyNodeLinks=false);
  /// Checks the graph data structure for internal consistency. 
  bool IsOk(const bool& ThrowExcept=true) const;
  /// Print the graph in a human readable form to an output stream OutF.
  void Dump(FILE *OutF=stdout) const;
  /// Returns a small graph on 5 nodes and 5 edges. 
  static PHGraph GetSmallGraph();
  void PrintEdge(const int EId);
  friend class TPt<THGraph>;
};

/*      SHOULD BE UNCOMMENTED IF THE CODE IS ADDED TO THE SNAP'S ORIGINAL LIBRARY
// set flags
namespace TSnap {
  template <> struct IsHyperGraph<THGraph> { enum { Val = 1 }; };
}
*/
