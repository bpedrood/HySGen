#ifndef Pedrood_CmntyHyper_h
#define Pedrood_CmntyHyper_h
#include <execinfo.h>
#include "Snap.h"
#include "loc_Snap.h"
#include <sstream> 
#include <stdlib.h> 


class THysgenUtil {
public:
  static void DumpCmtyVV(const TStr OutFNm, TVec<TIntV>& CmtyVV, TIntStrH& NIDNmH);
  static void DumpCmtyVH(const TStr OutFNm, TVec<TIntFltH>& CmtyVH, TIntStrH& NIDNmH);
  static PHGraph LoadEdgeList(const TStr& InFNm, TStrHash<TInt>& NodeNameH,
                              const TSsFmt SsFmt = ssfTabSep);
  
  template<class PHGraph>
  static double GetConductance(const PHGraph& Graph, const TIntSet& CmtyS, const int N2Edges);
  
  template<class PHGraph>
  static void GetNbhCom(const PHGraph& Graph, const int NID, TIntSet& NBCmtyS);
  
  template<class PHGraph>
  static void GetPhiNIdPrV(const PHGraph &G, TFltIntPrV &PhiNIdPrV, const int MinComSiz);
  
  static void GetBinLocs(const int& DecNum, TIntV& LocsV, const TIntV& NodMapV);
  
  static double Min(double	N1, double N2) {
    return N1 < N2 ? N1 : N2;
  }
};

class THysgen {
private:
  PHGraph G; ///< The graph to fit
  TVec<TIntFltH> S; ///< Membership of node \c NId to each of the communities
  ///< Product of the membership of all the nodes in edge \c EId to community \c CId
  ///< value saved for efficient calculation.
  THash<TInt,TIntFltH> ProbEdgCommHH;
  ///< Prob that no community generates Edge \c EId. Used for efficient calculation.
  THash<TInt,TFlt> ProbNotEdgH;
  ///< Prob that SOME community generates Edge \c EId.
  ///< Used to compensate when numerical precision of the "double" datatype
  ///< is not sufficient for error-free computation of ProbNotEdgH.
  THash<TInt,TFlt> ProbEdgH;
  ///< Auxiliary Matrix of (C+1, C+1) to help in the proposed DP to compute the Edge probability
  TVec<TFltV> AuxDPEdgVV;
  TFltV ProbENoiseV; ///< Prob that an edge \c EId is generated by noise (leak)
  TRnd Rnd; ///< Random number generator
  TIntSet NIdToIdx;
  TInt NumComs; ///< The total number of communities
  TIntV NumCIdNV; ///< Number of nodes (members) in each community
  TIntSet HOEIDS; //Set of edges to hold out for cross validation
  TVec<TIntSet> HONEIdsV; ///< For each node, this keeps the ID of neighboring hold out Edges
  TVec<TIntSet> HONNIdsV; ///< For each node, this keeps the ID of Nodes in the neighboring hold out Edges
  TVec<TIntSet> HOKIDSV; ///< set of attribute index (k) to hold out
  TInt MinTayN; ///< Minimum Number of terms in Taylor expansion that we like to expand.
  TInt TayN; ///< MAXIMUM Number of terms in Taylor expansion that we like to expand.
  TFlt TayThresh; ///< Minimum value that a term (n) in Taylor expansion must have to be stored.
  TVec<TFltV> SumPrPsblEdgesPow_nVV;
  TFlt SNoise;
  TFlt ThreshLearnRate;
  TInt DebugPrcCmpCnt;
  TInt DebugPrcApproxCmpCnt;
  TIntIntH DebugPrcCmpCntH;
public:
  TFlt InitVal;   ///< default value of S for community initialization
  TFlt InitNullS; ///< Default membership value for all the nodes to the noise community
  TFlt MinVal;    ///< minimum value of S
  TFlt MaxVal;    ///< maximum value of S (for numerical reasons)
  TFlt NegWgt;
  TFlt RegCoef;   ///< L1 regularization coefficient
  TFlt PrNoCom;
private:
  void GetUpdatedNodP(TIntFltH &SNew, const int &UID, const TIntFltH GradUH,
                      double &StepSize);
  bool AcceptStepSA(const int &UID, const TIntFltH &SNew, const int &Iter,
                    const int &MaxIter, const double &SAParamK);
public:
  THysgen(const PHGraph& GraphPt, const int& InitComs,
          const int RndSeed = 0, const double _InitVal = 0.1, const double InitNulS = 0.03,
          const double NoiseConstS = 0.01):
          Rnd(RndSeed), RegCoef(1),MinVal(0.0), MaxVal(1.0), NegWgt(1.0), InitVal(_InitVal),
          TayN(50), MinTayN(10), TayThresh(0.00001),
          InitNullS(InitNulS), SNoise(NoiseConstS),
          DebugPrcCmpCnt(0), DebugPrcApproxCmpCnt(0), DebugPrcCmpCntH(5) {
    SNoise = (InitNulS>0.0) ? THysgenUtil::Min(NoiseConstS, InitNulS / 2) : NoiseConstS;
    ThreshLearnRate = TayThresh / GraphPt->GetNodes();
    SetGraph(GraphPt); ComInit(InitComs);
  }
  void Save(TSOut& SOut) {
    G->Save(SOut);
    S.Save(SOut);
    ProbEdgCommHH.Save(SOut);
    ProbNotEdgH.Save(SOut);
    ProbEdgH.Save(SOut);
    AuxDPEdgVV.Save(SOut);
    ProbENoiseV.Save(SOut);
    NIdToIdx.Save(SOut);
    RegCoef.Save(SOut);
    NumComs.Save(SOut);
    NumCIdNV.Save(SOut);
    HONEIdsV.Save(SOut);
    HONNIdsV.Save(SOut);
    HOKIDSV.Save(SOut);
    TayN.Save(SOut);
    MinTayN.Save(SOut);
    TayThresh.Save(SOut);
    SumPrPsblEdgesPow_nVV.Save(SOut);
    InitVal.Save(SOut);
    MinVal.Save(SOut);
    MaxVal.Save(SOut);
    NegWgt.Save(SOut);
    PrNoCom.Save(SOut);
    InitNullS.Save(SOut);
  }
  void Load(TSIn& SIn, const int& RndSeed = 0) {
    G->Load(SIn);
    S.Load(SIn);
    ProbEdgCommHH.Load(SIn);
    ProbNotEdgH.Load(SIn);
    ProbEdgH.Load(SIn);
    AuxDPEdgVV.Load(SIn);
    ProbENoiseV.Load(SIn);
    NIdToIdx.Load(SIn);
    RegCoef.Load(SIn);
    NumComs.Load(SIn);
    NumCIdNV.Load(SIn);
    HONEIdsV.Load(SIn);
    HONNIdsV.Load(SIn);
    HOKIDSV.Load(SIn);
    TayN.Load(SIn);
    MinTayN.Load(SIn);
    TayThresh.Load(SIn);
    SumPrPsblEdgesPow_nVV.Load(SIn);
    InitVal.Load(SIn);
    MinVal.Load(SIn);
    MaxVal.Load(SIn);
    NegWgt.Load(SIn);
    PrNoCom.Load(SIn);
    InitNullS.Load(SIn);
  }

  void SetGraph(const PHGraph& GraphPt);
  void SetRegCoef(const double _RegCoef) { RegCoef = _RegCoef; }
  double GetRegCoef() { return RegCoef; }
  void LoadComInit(const TStr& InFNm, TSsFmt SsFmt=ssfTabSep);
  void ComInit(const int InitComs, const int MinComSiz=5, const double PerturbDens=0.0);
  void UniformComInit(const int InitComs);
  void RandomComPerturb(double Density = 1.0);
  
  /// Initialize with the neighborhood communities (Gleich et.al. KDD'12)
  void NeighborComInit(const int MinComSiz, const bool& IsInit = false);
  void NeighborComInit(TFltIntPrV& PhiNIdPrV, const bool& IsInit = false);
  
  TInt GetNumComs() { return NumComs; }
  void SetCmtyVV(const TVec<TIntV>& CmtyVV);
  double Likelihood();
  double Likelihood(const int UID, const TIntFltH& SU);
  double LikelihoodForRow(const int UID);
  double LikelihoodForRow(const int UID, const TIntFltH& SU, const bool CmprDirct_vs_Taylor=false);
  void GradientForRow(const int UId, TIntFltH& GradNod, const TIntSet& CIDSet, const bool CmprDirct_vs_Taylor=false);
  
  /// Dump community affiliation into a text file with node names
  void GetCmtyVV(TVec<TIntFltH>& CmtyVH, TVec<TIntV>& CmtyVV, TVec<TFltV>& WckVV,
                 const double Thres, const int MinSz = 3);
  void GetCmtyVV(TVec<TIntFltH>& CmtyVH, TVec<TIntV>& CmtyVV, TVec<TFltV>& WckVV,
                 const int MinSz = 3) {
    printf("\n+-+-+-+- Threshold = %f -+-+-+-+\n",
           sqrt(2.0 * (double) G->GetEdges() / G->GetNodes() / G->GetNodes()));
    GetCmtyVV(CmtyVH, CmtyVV, WckVV,
              sqrt(2.0 * (double) G->GetEdges() / G->GetNodes() / G->GetNodes()), MinSz);
  }
  void GetCmtyVV(TVec<TIntFltH>& CmtyVH, TVec<TIntV>& CmtyVV,
                 const double Thres, const int MinSz = 3) {
    TVec<TFltV> TmpV;
    GetCmtyVV(CmtyVH, CmtyVV, TmpV, Thres, MinSz);
  }
  void GetCmtyVV(TVec<TIntFltH>& CmtyVH, TVec<TIntV>& CmtyVV) {
    TVec<TFltV> TmpVV;
    GetCmtyVV(CmtyVH, CmtyVV, TmpVV, sqrt(2.0 * (double) G->GetEdges() / G->GetNodes() / G->GetNodes()), 3);
  }
  void GetCmtyVVUnSorted(TVec<TIntV>& CmtyVV);
  void GetCmtyVVUnSorted(TVec<TIntV>& CmtyVV, const double Thres, const int MinSz = 3);
  
  double GetStepSizeByLineSearch(const int UID, const TIntFltH &DeltaH,
                                 TIntFltH &SearchVecH, const double &stepSize,
                                 const double &CtrlParam, const double &ReductionRatio,
                                 const int MaxIter);
  
  int MLEGradAscent(const double& Thres, const int& MaxIter, const TStr PlotNm,
                    const double StepSize = 0.3, const double StepCtrlParam = 0.3, const double StepReductionRatio = 0.3);
  
  void inline GetNCom(TIntFltH& NIdH, const int& NID) {
    NIdH = S[NID];
  }
  double inline GetNCom(const int& NID, const int& CID) {
    if (S[NID].IsKey(CID)) {
      return S[NID].GetDat(CID);
    } else {
      return 0.0;
    }
  }
  
  void inline DelNCom(const int &UID, const int &CID) {
    if (S[UID].IsKey(CID)) {
      UpdatePrAllEdgesS(UID, CID, 0.0);
      UpdateUEdgesProb(UID, CID, 0.0);
      S[UID].DelKey(CID);
      NumCIdNV[CID] --;
    }
  }
  
  void inline AddNCom(const TIntV& NIdV, const int& CID, const double& Val, const bool& IsInit = false) {
    for (int ui = 0; ui < NIdV.Len(); ui++){
      AddNCom(NIdV[ui], CID, Val, IsInit);
    }
  }
  void inline AddNCom(const int& UId, const int& CID, const double& Val, const bool& IsInit = false) {
    if (Val == 0.0) {
      DelNCom(UId, CID);
      return;
    }
    if (!IsInit){ 
      UpdatePrAllEdgesS(UId, CID, Val);
      UpdateUEdgesProb(UId, CID, Val);
    }
    if (S[UId].IsKey(CID)) {
      NumCIdNV[CID] --;
    }
    S[UId].AddDat(CID, Val);
    NumCIdNV[CID] ++;
  }
  
  double inline GetPrE(const int &EId) {
    if (ProbNotEdgH.GetDat(EId)() == -1.0) {
      return ProbEdgH.GetDat(EId)();
    }
    return 1.0 - ProbNotEdgH.GetDat(EId);
  }
  /// Computes Prob(EId) when SU is different from S[UId].
  /// @param PrEOutCH is an output vector of \prod_{v \in e} S_vc for each c
  double GetPrE(const int &EId, const int &UId, TIntFltH &PrEOutCH,
                       const TIntFltH &SU);
  
  double GetPrEPrecisionApprox(const TIntFltH& ECH, TVec<TFltV>& DPMatVV,
                                const double PrENoise);
  
  double GetPrEPrecision(const TIntFltH& ECH, TVec<TFltV>& DPMatVV, const double PrENoise, TInt LineNo);
  
  double inline GetENoiseProb(const int Size) {
    if (Size > ProbENoiseV.Len()) { AddENoiseProb(Size); }
    return ProbENoiseV[Size];
  }
  
  void inline AddENoiseProb(const int Size) {
    for (int i = ProbENoiseV.Len(); i <= Size; i++) {
      ProbENoiseV[i] = ProbENoiseV[i-1] * SNoise;
    }
  }
  
  template <class TInt1, class TInt2>
  double inline GetECom(const TInt1 &EId, const TInt2 &CId) {
    if (ProbEdgCommHH.IsKey(EId) && ProbEdgCommHH.GetDat(EId).IsKey(CId)) {
      return ProbEdgCommHH.GetDat(EId).GetDat(CId);
    }
    else {
      return 0.0;
    }
  }
  
  void inline AddECom(const int& EId, const TIntFltH& ProdH) {
    double PrENoise = GetENoiseProb(G->GetEI(EId).Len());
    double PrNotE = 1 - PrENoise;
    if (ProdH.Len() != 0) {
      ProbEdgCommHH.AddDat(EId, ProdH);
      for (TIntFltH::TIter HI = ProdH.BegI(); HI < ProdH.EndI(); HI++) {
        PrNotE *= 1.0 - HI.GetDat();
      }
    }
    if (PrNotE >= 1.0 && SNoise > 0) {
      ProbEdgH.AddDat(EId, GetPrEPrecision(ProdH, AuxDPEdgVV, PrENoise, 730));
      ProbNotEdgH.AddDat(EId,-1);
      return;
    }
    ProbNotEdgH.AddDat(EId,PrNotE);
  }
  void inline AddECom(const int& EId, const int& CId, const double& PrECNew) {
    double PrECOld = GetECom(EId, CId);
    if (! ProbEdgCommHH.IsKey(EId)) {
      TIntFltH EmptyEProdH(NumComs);
      ProbEdgCommHH.AddDat(EId, EmptyEProdH);
    }
    ProbEdgCommHH.GetDat(EId).AddDat(CId, PrECNew);
    UpdateProbNotEdgH(EId, PrECNew, PrECOld);
  }
  
  void inline DelECom(const int& EId, const int& CId) {
    double PrECOld = GetECom(EId, CId);
    if (PrECOld > 0) {
      ProbEdgCommHH.GetDat(EId).DelKey(CId);
      UpdateProbNotEdgH(EId, 0.0, PrECOld);
    }
  }
  
  void InitEdgeProb();
  void UpdateUEdgesProb(const int& UId, const int& CId, const double& SUNew);
  void InitPrAllEdgesS(const double& DefVal, const bool& IsEqualComms=false);
  void UpdatePrAllEdgesS(const int &UID, const int &CID, const double& SNNew);
  void UpdatePrAllEdgesS(TFltV &PsiV, const int &UID, const int &CID, const double& SNNew, const bool IsApplyChange);
  void UpdatePrAllEdgesS(const int &UID, const int &CID, const TFltV& SNodNewV){
    UpdatePrAllEdgesS(UID, CID, SNodNewV[CID]);
  }
  double PredictAllCEdgesS_direct(const int &UID, const int &CID,
                                  const bool IsLikelihood = false, const bool Verbose= false);
  double PredictAllCEdgesS(const int &UID, const int &CID,
                           const bool IsLikelihood = false, const bool Verbose=false);
  double PredictAllCEdgesS(const int &UID, const int &CID, const double& SNNew,
                           const bool IsLikelihood = true, const bool Verbose=false);
  void inline UpdateProbNotEdgH(const int &EId, const double &PrECNew,
                                const double &PrECOld) {
    double PrENoise = GetENoiseProb(G->GetEI(EId).Len());    
    if ((PrECNew>0.0 && 1-PrECNew>=1.0) || (ProbNotEdgH.IsKey(EId) && ProbNotEdgH.GetDat(EId)()==-1.0))  {      
      TIntFltH& ProdH = ProbEdgCommHH.GetDat(EId); 
      ProbEdgH.AddDat(EId, GetPrEPrecisionApprox(ProdH, AuxDPEdgVV, PrENoise));
      ProbNotEdgH.AddDat(EId,-1);
      return;
    }
    double PrNotE;
    if (PrECOld < 1.0){
      PrNotE = ProbNotEdgH.GetDat(EId) * (1.0 - PrECNew) / (1.0 - PrECOld);
    } else {
      PrNotE = 1.0 - PrENoise;
      for (TIntFltH::TIter HI = ProbEdgCommHH.GetDat(EId).BegI();
           HI < ProbEdgCommHH.GetDat(EId).EndI(); HI++) {
        PrNotE *= 1.0 - HI.GetDat(); 
      }
    }
    ProbNotEdgH.AddDat(EId, PrNotE);
    if (PrNotE >= 1) {
      TIntFltH& ProdH = ProbEdgCommHH.GetDat(EId); 
      ProbEdgH.AddDat(EId, GetPrEPrecision(ProdH, AuxDPEdgVV, PrENoise,798));
      ProbNotEdgH.AddDat(EId,-1);
    }
  }
  
  double inline DotProduct(const TIntFltH& UV, const TIntFltH& VV) {
    double DP = 0;
    if (UV.Len() > VV.Len()) {
      for (TIntFltH::TIter HI = UV.BegI(); HI < UV.EndI(); HI++) {
        if (VV.IsKey(HI.GetKey())) {
          DP += VV.GetDat(HI.GetKey()) * HI.GetDat();
        }
      }
    } else {
      for (TIntFltH::TIter HI = VV.BegI(); HI < VV.EndI(); HI++) {
        if (UV.IsKey(HI.GetKey())) {
          DP += UV.GetDat(HI.GetKey()) * HI.GetDat();
        }
      }
    }
    return DP;
  }
  
  double inline DotProduct(const int& UID, const int& VID) {
    return DotProduct(S[UID], S[VID]);
  }
  void inline Normalize(const TIntFltH& UV, TIntFltH& UVNmd) {
    double Nrm = Norm2(UV);
    for (TIntFltH::TIter HI = UV.BegI(); HI < UV.EndI(); HI++) {
      UVNmd.AddDat(HI.GetKey(), HI.GetDat()/Nrm);
    }
  }
  void inline Normalize(TIntFltH& UV) {
    double Nrm = Norm2(UV);
    for (TIntFltH::TIter HI = UV.BegI(); HI < UV.EndI(); HI++) {
      UV.AddDat(HI.GetKey(), HI.GetDat()/Nrm);
    }
  }
  void inline NormalizeIfLarge(const TIntFltH& UV, TIntFltH& UVNmd) {
    double Nrm = Norm2(UV);
    if (Nrm >1) { Normalize(UV, UVNmd); }
    else { 
      for (TIntFltH::TIter HI = UV.BegI(); HI < UV.EndI(); HI++) {
        UVNmd.AddDat(HI.GetKey(), HI.GetDat());
      }
    }
  }
  double inline Sum(const TIntFltH& UV) {
    double N = 0.0;
    for (TIntFltH::TIter HI = UV.BegI(); HI < UV.EndI(); HI++) {
      N += HI.GetDat();
    }
    return N;
  }
  double inline Norm2(const TIntFltH& UV) {
    double N = 0.0;
    for (TIntFltH::TIter HI = UV.BegI(); HI < UV.EndI(); HI++) {
      N += HI.GetDat() * HI.GetDat();
    }
    return sqrt(N);
  }
  double inline Sigmoid(const double X) {
    return 1.0 / ( 1.0 + exp(-X));
  }
  
  void PrintComms() {
    printf("List of communities:\n");
    for (int c = 0; c < GetNumComs(); c++) {
      printf("comm %d has %d members:\t\t", c, NumCIdNV[c].Val);
      for (THGraph::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
        if (GetNCom(NI.GetId(),c) > InitNullS) {
          printf("%s,%0.2f\t",NI.GetName().CStr(), GetNCom(NI.GetId(),c));
        }
      }
      printf("\n");
    }
  }
  template <class TObj2d> void PrintObj(const TObj2d& O2D) {
    printf("Obj = [");
    for (typename TObj2d::TIter i = O2D.BegI(); i<O2D.EndI(); i++) {
      printf("(%d,%f), ", i.GetKey(), i.GetDat());
    }
    printf("]\n");
  }
  template <class TObj2d> void PrintObj(const TObj2d& O2D, double coef) {
    printf("Obj = [");
    for (typename TObj2d::TIter i = O2D.BegI(); i<O2D.EndI(); i++) {
      printf("(%d,%f), ", i.GetKey(), i.GetDat()*coef);
    }
    printf("]\n");
  }
  
};


#endif
