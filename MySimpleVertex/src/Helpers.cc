#include "MyVertexFinder.h"

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMatrixTSym.h"

// ACTS
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include <Acts/Vertexing/ImpactPointEstimator.hpp>
#include "Acts/Vertexing/FullBilloirVertexFitter.hpp"

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

using namespace Acts;
using namespace Acts::UnitLiterals;

/*  Helpers and Utilities */
const Acts::MagneticFieldContext &MyVertexFinder::magneticFieldContext() const
{
    return m_magneticFieldContext;
}

const Acts::GeometryContext &MyVertexFinder::geometryContext() const
{
    return m_geometryContext;
}

std::shared_ptr<Acts::MagneticFieldProvider> MyVertexFinder::magneticField() const
{
    return m_magneticField;
}

double MyVertexFinder::projSV_PV(const TVector3 &SV, const TVector3 &PV, const TLorentzVector &Direction) const
{
    TVector3 SV_PV(SV.x() - PV.x(), SV.y() - PV.y(), SV.z() - PV.z());
    return Direction.Vect().Unit() * SV_PV.Unit();
}

double MyVertexFinder::MomProjDist(const TVector3 &SecVrt, const TVector3 &primVrt, const TLorentzVector &Mom) const
{
    TVector3 vv = SecVrt - primVrt;
    return (vv.x() * Mom.X() + vv.y() * Mom.Y() + vv.z() * Mom.Z()) / Mom.P();
}

void MyVertexFinder::printWrkSet(const std::vector<WrkVrt> *WrkVrtSet, const std::string &name) const
{
    int nGoodV = 0;
    for (int iv = 0; iv < (int)WrkVrtSet->size(); iv++)
    {
        streamlog_out(DEBUG0) << name
                              << "= " << (*WrkVrtSet)[iv].vertex.position()[0]
                              << ", " << (*WrkVrtSet)[iv].vertex.position()[1]
                              << ", " << (*WrkVrtSet)[iv].vertex.position()[2]
                              << " NTrk=" << (*WrkVrtSet)[iv].selTrk.size()
                              << " is good=" << std::boolalpha << (*WrkVrtSet)[iv].Good << std::noboolalpha
                              << "  Chi2=" << (*WrkVrtSet)[iv].vertex.fitQuality().first
                              << "  Mass=" << (*WrkVrtSet)[iv].vertexMom.M()
                              << "  detached=" << (*WrkVrtSet)[iv].detachedTrack
                              << "  proj.dist=" << (*WrkVrtSet)[iv].projectedVrt
                              << " trk=";
        for (int kk = 0; kk < (int)(*WrkVrtSet)[iv].selTrk.size(); kk++)
        {
            streamlog_out(DEBUG0) << ", " << (*WrkVrtSet)[iv].selTrk[kk];
        }
        /*
        for (int kk = 0; kk < (int)(*WrkVrtSet)[iv].selTrk.size(); kk++)
        {
            streamlog_out(DEBUG0) << ", " << momAtVrt((*WrkVrtSet)[iv].trkAtVrt[kk]).Perp();
        }
        */
        streamlog_out(DEBUG0) << std::endl;
        if ((*WrkVrtSet)[iv].Good)
            nGoodV++;
    }
    streamlog_out(DEBUG0) << name << " N=" << nGoodV << std::endl;
}

int MyVertexFinder::nTrkCommon(std::vector<WrkVrt> *wrkVrtSet, int V1, int V2)
    const
{
    int nTrk_V1 = (*wrkVrtSet).at(V1).selTrk.size();
    if (nTrk_V1 < 2)
        return 0; /* Bad vertex */
    int nTrk_V2 = (*wrkVrtSet).at(V2).selTrk.size();
    if (nTrk_V2 < 2)
        return 0; /* Bad vertex */
    int nTrkCom = 0;
    if (nTrk_V1 < nTrk_V2)
    {
        for (int i = 0; i < nTrk_V1; i++)
        {
            int trk = (*wrkVrtSet)[V1].selTrk[i];
            for (int j = 0; j < nTrk_V2; j++)
            {
                if (trk == (*wrkVrtSet)[V2].selTrk[j])
                {
                    nTrkCom++;
                    break;
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < nTrk_V2; i++)
        {
            int trk = (*wrkVrtSet)[V2].selTrk[i];
            for (int j = 0; j < nTrk_V1; j++)
            {
                if (trk == (*wrkVrtSet)[V1].selTrk[j])
                {
                    nTrkCom++;
                    break;
                }
            }
        }
    }
    return nTrkCom;
}

double MyVertexFinder::mergeAndRefitVertices(WrkVrt &v1, WrkVrt &v2, WrkVrt &newvrt,
                                             ROOT::VecOps::RVec<Acts::BoundTrackParameters> &AllTrackList)
    const
{
    if (!v1.Good)
        return -1.; // bad vertex
    if (!v2.Good)
        return -1.; // bad vertex

    newvrt.Good = true;
    int nTrk_V1 = v1.selTrk.size();
    int nTrk_V2 = v2.selTrk.size();
    streamlog_out(DEBUG0) << " v1 " << nTrk_V1 << " v2 " << nTrk_V2 << std::endl;

    newvrt.selTrk.resize(nTrk_V1 + nTrk_V2);
    std::copy(v1.selTrk.begin(), v1.selTrk.end(), newvrt.selTrk.begin());
    std::copy(v2.selTrk.begin(), v2.selTrk.end(), newvrt.selTrk.begin() + nTrk_V1);

    std::deque<long int>::iterator TransfEnd;
    sort(newvrt.selTrk.begin(), newvrt.selTrk.end());
    TransfEnd = unique(newvrt.selTrk.begin(), newvrt.selTrk.end());
    newvrt.selTrk.erase(TransfEnd, newvrt.selTrk.end());

    using Propagator = Acts::Propagator<Acts::EigenStepper<>>;
    using Estimator = Acts::ImpactPointEstimator<Acts::BoundTrackParameters, Propagator>;

    // Set up EigenStepper
    Acts::EigenStepper<> vstepper(magneticField());

    using Linearizer = Acts::HelicalTrackLinearizer<Propagator>;
    Linearizer::Config ltConfig(magneticField(), std::make_shared<Propagator>(std::move(vstepper)));
    Linearizer linearizer(ltConfig);

    // Set up Billoir Vertex Fitter
    using VertexFitter = Acts::FullBilloirVertexFitter<Acts::BoundTrackParameters, Linearizer>;
    VertexFitter::Config vertexFitterCfg;
    VertexFitter billoirFitter(vertexFitterCfg);
    VertexFitter::State state(magneticField()->makeCache(magneticFieldContext()));
    Acts::VertexingOptions<Acts::BoundTrackParameters> vfOptions(geometryContext(), magneticFieldContext());

    std::vector<const Acts::BoundTrackParameters *> fitTrackList;
    TLorentzVector sum_p4;
    streamlog_out(DEBUG0) << " vmerged " << newvrt.selTrk.size() << std::endl;

    for (int i = 0; i < (int)newvrt.selTrk.size(); i++)
    {
        fitTrackList.push_back(&AllTrackList[newvrt.selTrk[i]]);
        TVector3 p3_i = {AllTrackList[newvrt.selTrk[i]].momentum()[0], AllTrackList[newvrt.selTrk[i]].momentum()[1], AllTrackList[newvrt.selTrk[i]].momentum()[2]};
        TLorentzVector p4_i;
        p4_i.SetPtEtaPhiM(AllTrackList[newvrt.selTrk[i]].transverseMomentum(), p3_i.PseudoRapidity(), p3_i.Phi(), m_massPi);
        sum_p4 += p4_i;
    }

    auto result = billoirFitter.fit(fitTrackList, linearizer, vfOptions, state);
    if (not result.ok())
        return -1.;

    Acts::Vertex<Acts::BoundTrackParameters> fittedVertex = result.value();

    if (fittedVertex.fitQuality().first > 500.)
        return -1.;
    newvrt.vertex = fittedVertex;
    newvrt.Good = true;
    newvrt.vertexMom = sum_p4;
    newvrt.chi2 = fittedVertex.fitQuality().first;

    return TMath::Prob(newvrt.chi2, fittedVertex.fitQuality().second);
}

double MyVertexFinder::refitVertex(WrkVrt &vrt,
                                   ROOT::VecOps::RVec<Acts::BoundTrackParameters> &selectedTracks) const
{
    int i, j;
    int nth = vrt.selTrk.size();

    if (nth < 2)
        return -1.;

    std::vector<const Acts::BoundTrackParameters *> fitTrackList;
    TLorentzVector sum_p4;
    for (int i = 0; i < nth; i++)
    {
        fitTrackList.push_back(&selectedTracks[vrt.selTrk[i]]);
        TVector3 p3_i = {selectedTracks[vrt.selTrk[i]].momentum()[0], selectedTracks[vrt.selTrk[i]].momentum()[1], selectedTracks[vrt.selTrk[i]].momentum()[2]};
        TLorentzVector p4_i;
        p4_i.SetPtEtaPhiM(selectedTracks[vrt.selTrk[i]].transverseMomentum(), p3_i.PseudoRapidity(), p3_i.Phi(), m_massPi);
        sum_p4 += p4_i;
    }

    vrt.Good = false;
    vrt.chi2PerTrk.resize(nth);
    // for (i = 0; i < nth; i++)
    //     vrt.chi2PerTrk[i] = 100000. + i; // VK safety

    using Propagator = Acts::Propagator<Acts::EigenStepper<>>;
    using Estimator = Acts::ImpactPointEstimator<Acts::BoundTrackParameters, Propagator>;

    // Set up EigenStepper
    Acts::EigenStepper<> vstepper(magneticField());

    using Linearizer = Acts::HelicalTrackLinearizer<Propagator>;
    Linearizer::Config ltConfig(magneticField(), std::make_shared<Propagator>(std::move(vstepper)));
    Linearizer linearizer(ltConfig);

    // Set up Billoir Vertex Fitter
    using VertexFitter = Acts::FullBilloirVertexFitter<Acts::BoundTrackParameters, Linearizer>;
    VertexFitter::Config vertexFitterCfg;
    VertexFitter billoirFitter(vertexFitterCfg);
    VertexFitter::State state(magneticField()->makeCache(magneticFieldContext()));
    Acts::VertexingOptions<Acts::BoundTrackParameters> vfOptions(geometryContext(), magneticFieldContext());

    auto result = billoirFitter.fit(fitTrackList, linearizer, vfOptions, state);
    if (not result.ok())
    {
        vrt.Good = false;
        return -1.;
    }

    Acts::Vertex<Acts::BoundTrackParameters> fittedVertex = result.value();
    vrt.vertex = fittedVertex;
    vrt.Good = true;
    vrt.vertexMom = sum_p4;
    vrt.chi2 = fittedVertex.fitQuality().first;

    return TMath::Prob(vrt.chi2, fittedVertex.fitQuality().second);
}

double MyVertexFinder::refineVerticesWithCommonTracks(WrkVrt &v1, WrkVrt &v2,
                                                      ROOT::VecOps::RVec<Acts::BoundTrackParameters> &allTrackList)
    const
{
    if (!v1.Good || !v2.Good)
        return -1.; // bad vertex check for safety

    int ntv1 = v1.selTrk.size();
    int ntv2 = v2.selTrk.size();
    if (ntv1 < 2 || ntv2 < 2)
        return -1.; // for safety
    double prb_v1 = TMath::Prob(v1.chi2, 2 * ntv1 - 3);
    double prb_v2 = TMath::Prob(v2.chi2, 2 * ntv2 - 3);
    WrkVrt *badV = &v1;
    WrkVrt *goodV = &v2;
    bool swap = false;
    //---- Good vertex selection
    if (prb_v1 > 0.01 && prb_v2 > 0.01)
    { // Both vertices are good Prob>1%
        if (ntv1 == ntv2)
        {
            if (prb_v1 > prb_v2)
                swap = true;
        } // If multiplicities are equal- select better Chi2
        else
        {
            if (ntv1 > ntv2)
                swap = true;
        }
    }
    if (prb_v1 < 0.01 && prb_v2 < 0.01)
    { // Both vertices are bad Prob<1%
        if (prb_v1 > prb_v2)
            swap = true; // select better Chi2
    }
    if (prb_v1 > 0.01 && prb_v2 < 0.01)
    { // Second vertex is bad Prob<1%
        if (prb_v1 > prb_v2)
            swap = true;
    }
    if (swap)
    {
        badV = &v2;
        goodV = &v1;
    }
    int badVNtrk = (*badV).selTrk.size();
    //-----------------
    unsigned int it = 0;
    while (it < (*badV).selTrk.size())
    {
        int trk = (*badV).selTrk[it];
        if (std::find((*goodV).selTrk.begin(), (*goodV).selTrk.end(), trk) != (*goodV).selTrk.end())
        {
            (*badV).selTrk.erase((*badV).selTrk.begin() + it);
            (*badV).detachedTrack = trk;
        }
        else
            it++;
    }
    if ((*badV).selTrk.size() < 2)
    {
        (*badV).Good = false;
        if ((*badV).selTrk.size() == 1 && m_multiWithOneTrkVrt)
        { // Special case if 1-track vertices are allowed
            (*badV).vertexCharge = allTrackList[(*badV).selTrk.at(0)].charge();
            TVector3 p3_i = {allTrackList[(*badV).selTrk.at(0)].momentum()[0], allTrackList[(*badV).selTrk.at(0)].momentum()[1], allTrackList[(*badV).selTrk.at(0)].momentum()[2]};
            TLorentzVector p4_i;
            p4_i.SetPtEtaPhiM(allTrackList[(*badV).selTrk.at(0)].transverseMomentum(), p3_i.PseudoRapidity(), p3_i.Phi(), m_massPi);
            (*badV).vertexMom = p4_i;
            if (badVNtrk >= 2)
                (*badV).Good = true; // For 1-track vertices
        }
        return -1.;
    }
    return refitVertex((*badV), allTrackList);
}

double MyVertexFinder::vrtVrtDist(const Acts::Vertex<Acts::BoundTrackParameters> &vrt1, const Acts::Vertex<Acts::BoundTrackParameters> &vrt2) const
{
    double distx = vrt1.position()[0] - vrt2.position()[0];
    double disty = vrt1.position()[1] - vrt2.position()[1];
    double distz = vrt1.position()[2] - vrt2.position()[2];

    auto vrtErr1 = vrt1.covariance();
    auto vrtErr2 = vrt2.covariance();

    TMatrixTSym<double> primCovMtx(3); // Create
    primCovMtx(0, 0) = vrtErr1(0, 0) + vrtErr2(0, 0);
    primCovMtx(0, 1) = primCovMtx(1, 0) = vrtErr1(1, 0) + vrtErr2(1, 0);
    primCovMtx(1, 1) = vrtErr1(1, 1) + vrtErr2(1, 1);
    primCovMtx(0, 2) = primCovMtx(2, 0) = vrtErr1(2, 0) + vrtErr2(2, 0);
    primCovMtx(1, 2) = primCovMtx(2, 1) = vrtErr1(2, 1) + vrtErr2(2, 1);
    primCovMtx(2, 2) = vrtErr1(2, 2) + vrtErr2(2, 2);

    TMatrixTSym<double> wgtMtx = primCovMtx.Invert();

    double signif =
        distx * wgtMtx(0, 0) * distx + disty * wgtMtx(1, 1) * disty + distz * wgtMtx(2, 2) * distz + 2. * distx * wgtMtx(0, 1) * disty + 2. * distx * wgtMtx(0, 2) * distz + 2. * disty * wgtMtx(1, 2) * distz;
    signif = std::sqrt(std::abs(signif));
    if (signif != signif)
        signif = 0.;
    return signif;
}

double MyVertexFinder::minVrtVrtDist(std::vector<WrkVrt> *wrkVrtSet, int &V1, int &V2, std::vector<double> &checked)
    const
{
    V1 = V2 = -1;
    double foundMinVrtDst = 1000000.;

    for (int iv = 0; iv < (int)wrkVrtSet->size() - 1; iv++)
    {
        if ((*wrkVrtSet)[iv].selTrk.size() < 2)
            continue; /* Bad vertices */
        if (!(*wrkVrtSet)[iv].Good)
            continue; /* Bad vertices */
        for (int jv = iv + 1; jv < (int)wrkVrtSet->size(); jv++)
        {
            if ((*wrkVrtSet)[jv].selTrk.size() < 2)
                continue; /* Bad vertices */
            if (!(*wrkVrtSet)[jv].Good)
                continue; /* Bad vertices */
            double tmp = std::abs((*wrkVrtSet)[iv].vertex.position()[0] - (*wrkVrtSet)[jv].vertex.position()[0]) + std::abs((*wrkVrtSet)[iv].vertex.position()[1] - (*wrkVrtSet)[jv].vertex.position()[1]) + std::abs((*wrkVrtSet)[iv].vertex.position()[2] - (*wrkVrtSet)[jv].vertex.position()[2]);
            if (tmp > 20.)
                continue;
            double tmpDst = vrtVrtDist((*wrkVrtSet)[iv].vertex, (*wrkVrtSet)[jv].vertex);
            if (std::find(checked.begin(), checked.end(), tmpDst) != checked.end())
                continue; // Already tried
            if (tmpDst < foundMinVrtDst)
            {
                foundMinVrtDst = tmpDst;
                V1 = iv;
                V2 = jv;
            }
        }
    }
    return foundMinVrtDst;
}
