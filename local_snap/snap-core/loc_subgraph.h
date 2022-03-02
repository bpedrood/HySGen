/*! \file subgraph.h
    \brief Functions and templates to generate subgraphs.
*/

/// Main namespace for all the Snap global entities.
namespace TSnap {

PHGraph GetSubGraph(const TPt<THGraph>& Graph, const TIntV& NIdV, const double& RenumberNodes=1.0);

} // namespace TSnap
