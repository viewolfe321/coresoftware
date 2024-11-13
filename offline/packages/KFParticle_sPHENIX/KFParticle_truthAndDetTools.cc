#include "KFParticle_truthAndDetTools.h"

#include <fun4all/Fun4AllReturnCodes.h> //for Fun4All to stop yelling at me

#include "KFParticle_Tools.h"  // for KFParticle_Tools

#include <g4eval/SvtxEvalStack.h>   // for SvtxEvalStack
#include <g4eval/SvtxTrackEval.h>   // for SvtxTrackEval
#include <g4eval/SvtxTruthEval.h>   // for SvtxTruthEval
#include <g4eval/SvtxVertexEval.h>  // for SvtxVertexEval

#include <globalvertex/SvtxVertex.h>         // for SvtxVertex
#include <globalvertex/SvtxVertexMap.h>      // for SvtxVertexMap, SvtxVer...
#include <trackbase/InttDefs.h>              // for getLadderPhiId, getLad...
#include <trackbase/MvtxDefs.h>              // for getChipId, getStaveId
#include <trackbase/TpcDefs.h>               // for getSectorId, getSide
#include <trackbase/TrkrCluster.h>           // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>  // for TrkrClusterContainer
#include <trackbase/TrkrDefs.h>              // for getLayer, getTrkrId
#include <trackbase_historic/SvtxPHG4ParticleMap_v1.h>
#include <trackbase_historic/SvtxTrack.h>     // for SvtxTrack, SvtxTrack::...
#include <trackbase_historic/SvtxTrackMap.h>  // for SvtxTrackMap, SvtxTrac...

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

#include <g4main/PHG4Particle.h>            // for PHG4Particle
#include <g4main/PHG4TruthInfoContainer.h>  // for PHG4TruthInfoContainer
#include <g4main/PHG4VtxPoint.h>            // for PHG4VtxPoint

#include <phhepmc/PHHepMCGenEvent.h>     // for PHHepMCGenEvent
#include <phhepmc/PHHepMCGenEventMap.h>  // for PHHepMCGenEventMap

#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/getClass.h>        // for getClass

#include <calobase/RawTowerGeomContainer.h> // for calo matching
#include <calobase/RawCluster.h> // for calo matching
#include <calobase/RawClusterUtility.h> // for calo matching
#include <calobase/RawTowerDefs.h> // for calo matching
#include <calobase/RawTowerGeom.h> // for calo matching
#include <calobase/RawTowerv2.h> // for calo matching
#include <calobase/RawTowerContainer.h> // for calo matching
#include <calobase/TowerInfoContainer.h> // for calo matching
#include <calobase/TowerInfo.h> // for calo matching
#include <calobase/TowerInfoDefs.h> // for calo matching

#include <KFParticle.h>  // for KFParticle
#include <TString.h>     // for TString, operator+
#include <TTree.h>       // for TTree

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <HepMC/GenEvent.h>   // for GenEvent::particle_con...
#include <HepMC/GenVertex.h>  // for GenVertex::particle_it...
#pragma GCC diagnostic pop

#include <HepMC/GenParticle.h>    // for GenParticle
#include <HepMC/IteratorRange.h>  // for parents
#include <HepMC/SimpleVector.h>   // for FourVector

#include <algorithm>  // for max, find
#include <cmath>      // for pow, sqrt
#include <cstdlib>    // for NULL, abs
#include <iostream>   // for operator<<, endl, basi...
#include <iterator>   // for end, begin
#include <map>        // for _Rb_tree_iterator, map
#include <memory>     // for allocator_traits<>::va...
#include <utility>    // for pair
#include <vector>     //for vector

class PHNode;

std::map<std::string, int> Use =
    {
        {"MVTX", 1},
        {"INTT", 1},
        {"TPC", 1},
        {"TPOT", 1},
        {"EMCAL", 0},
        {"OHCAL", 0},
        {"IHCAL", 0}};

KFParticle_truthAndDetTools::KFParticle_truthAndDetTools()
  : m_svtx_evalstack(nullptr)
{
}  // Constructor

KFParticle_truthAndDetTools::~KFParticle_truthAndDetTools() = default;  // Destructor

SvtxTrack *KFParticle_truthAndDetTools::getTrack(unsigned int track_id, SvtxTrackMap *trackmap)
{
  SvtxTrack *matched_track = nullptr;

  for (auto &iter : *trackmap)
  {
    if (iter.first == track_id)
    {
      matched_track = iter.second;
    }
  }

  return matched_track;
}

GlobalVertex *KFParticle_truthAndDetTools::getVertex(unsigned int vertex_id, GlobalVertexMap *vertexmap)
{
  GlobalVertex *matched_vertex = vertexmap->get(vertex_id);

  return matched_vertex;
}

PHG4Particle *KFParticle_truthAndDetTools::getTruthTrack(SvtxTrack *thisTrack, PHCompositeNode *topNode)
{
  /*
   * There are two methods for getting the truth rack from the reco track
   * 1. (recommended) Use the reco -> truth tables (requires SvtxPHG4ParticleMap). Introduced Summer of 2022
   * 2. Get truth track via nClusters. Older method and will work with older DSTs
   */

  PHG4Particle *particle = nullptr;

  PHNodeIterator nodeIter(topNode);
  PHNode *findNode = dynamic_cast<PHNode *>(nodeIter.findFirst("SvtxPHG4ParticleMap"));
  if (findNode)
  {
    findNode = dynamic_cast<PHNode *>(nodeIter.findFirst("G4TruthInfo"));
    if (findNode)
    {
      m_truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    }
    else
    {
      std::cout << "KFParticle truth matching: G4TruthInfo does not exist" << std::endl;
    }

    SvtxPHG4ParticleMap_v1 *dst_reco_truth_map = findNode::getClass<SvtxPHG4ParticleMap_v1>(topNode, "SvtxPHG4ParticleMap");

    std::map<float, std::set<int>> truth_set = dst_reco_truth_map->get(thisTrack->get_id());
    const auto &best_weight = truth_set.rbegin();
    int best_truth_id = *best_weight->second.rbegin();
    particle = m_truthinfo->GetParticle(best_truth_id);
  }
  else
  {
    std::cout << __FILE__ << ": SvtxPHG4ParticleMap not found, reverting to max_truth_particle_by_nclusters()" << std::endl;

    if (!m_svtx_evalstack)
    {
      m_svtx_evalstack = new SvtxEvalStack(topNode);
      // clustereval = m_svtx_evalstack->get_cluster_eval();
      // hiteval = m_svtx_evalstack->get_hit_eval();
      trackeval = m_svtx_evalstack->get_track_eval();
      trutheval = m_svtx_evalstack->get_truth_eval();
      vertexeval = m_svtx_evalstack->get_vertex_eval();
    }

    m_svtx_evalstack->next_event(topNode);

    particle = trackeval->max_truth_particle_by_nclusters(thisTrack);
  }
  return particle;
}

void KFParticle_truthAndDetTools::initializeTruthBranches(TTree *m_tree, int daughter_id, const std::string &daughter_number, bool m_constrain_to_vertex_truthMatch)
{
  m_tree->Branch(TString(daughter_number) + "_true_ID", &m_true_daughter_id[daughter_id], TString(daughter_number) + "_true_ID/I");
  if (m_constrain_to_vertex_truthMatch)
  {
    m_tree->Branch(TString(daughter_number) + "_true_IP", &m_true_daughter_ip[daughter_id], TString(daughter_number) + "_true_IP/F");
    m_tree->Branch(TString(daughter_number) + "_true_IP_xy", &m_true_daughter_ip_xy[daughter_id], TString(daughter_number) + "_true_IP_xy/F");
  }
  m_tree->Branch(TString(daughter_number) + "_true_px", &m_true_daughter_px[daughter_id], TString(daughter_number) + "_true_px/F");
  m_tree->Branch(TString(daughter_number) + "_true_py", &m_true_daughter_py[daughter_id], TString(daughter_number) + "_true_py/F");
  m_tree->Branch(TString(daughter_number) + "_true_pz", &m_true_daughter_pz[daughter_id], TString(daughter_number) + "_true_pz/F");
  m_tree->Branch(TString(daughter_number) + "_true_p", &m_true_daughter_p[daughter_id], TString(daughter_number) + "_true_p/F");
  m_tree->Branch(TString(daughter_number) + "_true_pT", &m_true_daughter_pt[daughter_id], TString(daughter_number) + "_true_pT/F");
  m_tree->Branch(TString(daughter_number) + "_true_EV_x", &m_true_daughter_vertex_x[daughter_id], TString(daughter_number) + "_true_EV_x/F");
  m_tree->Branch(TString(daughter_number) + "_true_EV_y", &m_true_daughter_vertex_y[daughter_id], TString(daughter_number) + "_true_EV_y/F");
  m_tree->Branch(TString(daughter_number) + "_true_EV_z", &m_true_daughter_vertex_z[daughter_id], TString(daughter_number) + "_true_EV_z/F");
  if (m_constrain_to_vertex_truthMatch)
  {
    m_tree->Branch(TString(daughter_number) + "_true_PV_x", &m_true_daughter_pv_x[daughter_id], TString(daughter_number) + "_true_PV_x/F");
    m_tree->Branch(TString(daughter_number) + "_true_PV_y", &m_true_daughter_pv_y[daughter_id], TString(daughter_number) + "_true_PV_y/F");
    m_tree->Branch(TString(daughter_number) + "_true_PV_z", &m_true_daughter_pv_z[daughter_id], TString(daughter_number) + "_true_PV_z/F");
  }
  m_tree->Branch(TString(daughter_number) + "_true_track_history_PDG_ID", &m_true_daughter_track_history_PDG_ID[daughter_id]);
  m_tree->Branch(TString(daughter_number) + "_true_track_history_PDG_mass", &m_true_daughter_track_history_PDG_mass[daughter_id]);
  m_tree->Branch(TString(daughter_number) + "_true_track_history_px", &m_true_daughter_track_history_px[daughter_id]);
  m_tree->Branch(TString(daughter_number) + "_true_track_history_py", &m_true_daughter_track_history_py[daughter_id]);
  m_tree->Branch(TString(daughter_number) + "_true_track_history_pz", &m_true_daughter_track_history_pz[daughter_id]);
  m_tree->Branch(TString(daughter_number) + "_true_track_history_pE", &m_true_daughter_track_history_pE[daughter_id]);
  m_tree->Branch(TString(daughter_number) + "_true_track_history_pT", &m_true_daughter_track_history_pT[daughter_id]);
}

void KFParticle_truthAndDetTools::fillTruthBranch(PHCompositeNode *topNode, TTree * /*m_tree*/, const KFParticle &daughter, int daughter_id, const KFParticle &kfvertex, bool m_constrain_to_vertex_truthMatch)
{
  float true_px, true_py, true_pz, true_p, true_pt;

  PHNodeIterator nodeIter(topNode);
  PHNode *findNode = dynamic_cast<PHNode *>(nodeIter.findFirst(m_trk_map_node_name_nTuple));
  if (findNode)
  {
    dst_trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trk_map_node_name_nTuple);
  }
  else
  {
    std::cout << "KFParticle truth matching: " << m_trk_map_node_name_nTuple << " does not exist" << std::endl;
  }
  findNode = dynamic_cast<PHNode *>(nodeIter.findFirst(m_vtx_map_node_name_nTuple));
  if (findNode)
  {
    dst_vertexmap = findNode::getClass<SvtxVertexMap>(topNode, m_vtx_map_node_name_nTuple);
  }
  else
  {
    std::cout << "KFParticle truth matching: " << m_vtx_map_node_name_nTuple << " does not exist" << std::endl;
  }
  auto globalvertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!globalvertexmap)
  {
    std::cout << "KFParticle truth matching: GlobalVertexMap does not exist" << std::endl;
  }
  track = getTrack(daughter.Id(), dst_trackmap);
  g4particle = getTruthTrack(track, topNode);

  bool isParticleValid = g4particle == nullptr ? false : true;

  true_px = isParticleValid ? (Float_t) g4particle->get_px() : 0.;
  true_py = isParticleValid ? (Float_t) g4particle->get_py() : 0.;
  true_pz = isParticleValid ? (Float_t) g4particle->get_pz() : 0.;
  true_p = sqrt(pow(true_px, 2) + pow(true_py, 2) + pow(true_pz, 2));
  true_pt = sqrt(pow(true_px, 2) + pow(true_py, 2));

  m_true_daughter_px[daughter_id] = true_px;
  m_true_daughter_py[daughter_id] = true_py;
  m_true_daughter_pz[daughter_id] = true_pz;
  m_true_daughter_p[daughter_id] = true_p;
  m_true_daughter_pt[daughter_id] = true_pt;
  m_true_daughter_id[daughter_id] = isParticleValid ? g4particle->get_pid() : 0;

  if (!m_svtx_evalstack)
  {
    m_svtx_evalstack = new SvtxEvalStack(topNode);
    // clustereval = m_svtx_evalstack->get_cluster_eval();
    // hiteval = m_svtx_evalstack->get_hit_eval();
    trackeval = m_svtx_evalstack->get_track_eval();
    trutheval = m_svtx_evalstack->get_truth_eval();
    vertexeval = m_svtx_evalstack->get_vertex_eval();
  }

  if (isParticleValid)
  {
    g4vertex_point = trutheval->get_vertex(g4particle);
  }

  m_true_daughter_vertex_x[daughter_id] = isParticleValid ? g4vertex_point->get_x() : 0.;
  m_true_daughter_vertex_y[daughter_id] = isParticleValid ? g4vertex_point->get_y() : 0.;
  m_true_daughter_vertex_z[daughter_id] = isParticleValid ? g4vertex_point->get_z() : 0.;

  if (m_constrain_to_vertex_truthMatch)
  {
    // Calculate true DCA
    GlobalVertex *recoVertex = getVertex(kfvertex.Id(), globalvertexmap);
    auto svtxviter = recoVertex->find_vertexes(GlobalVertex::SVTX);
    // check that it contains a track vertex
    if (svtxviter == recoVertex->end_vertexes())
    {
      std::cout << "Have a global vertex with no track vertex... shouldn't happen in KFParticle_truthAndDetTools::fillTruthBranch..." << std::endl;
    }

    auto svtxvertexvector = svtxviter->second;
    SvtxVertex *svtxvertex = nullptr;

    for (auto &vertex_iter : svtxvertexvector)
    {
      svtxvertex = dst_vertexmap->find(vertex_iter->get_id())->second;
    }

    PHG4VtxPoint *truePoint = vertexeval->max_truth_point_by_ntracks(svtxvertex);

    if (truePoint == nullptr && isParticleValid)
    {
      PHG4Particle *g4mother = m_truthinfo->GetPrimaryParticle(g4particle->get_primary_id());
      truePoint = m_truthinfo->GetVtx(g4mother->get_vtx_id());  // Note, this may not be the PV for a decay with tertiaries
    }

    KFParticle trueKFParticleVertex;

    float f_vertexParameters[6] = {0};

    if (truePoint == nullptr)
    {
      std::cout << "KFParticle truth matching: This event has no PHG4VtxPoint information!\n";
      std::cout << "Your truth track DCA will be measured wrt a reconstructed vertex!" << std::endl;

      f_vertexParameters[0] = recoVertex->get_x();
      f_vertexParameters[1] = recoVertex->get_y();
      f_vertexParameters[2] = recoVertex->get_z();
    }
    else
    {
      f_vertexParameters[0] = truePoint->get_x();
      f_vertexParameters[1] = truePoint->get_y();
      f_vertexParameters[2] = truePoint->get_z();
    }

    float f_vertexCovariance[21] = {0};

    trueKFParticleVertex.Create(f_vertexParameters, f_vertexCovariance, 0, -1);

    KFParticle trueKFParticle;

    float f_trackParameters[6] = {m_true_daughter_vertex_x[daughter_id],
                                  m_true_daughter_vertex_y[daughter_id],
                                  m_true_daughter_vertex_z[daughter_id],
                                  true_px,
                                  true_py,
                                  true_pz};

    float f_trackCovariance[21] = {0};

    trueKFParticle.Create(f_trackParameters, f_trackCovariance, 1, -1);

    m_true_daughter_ip[daughter_id] = trueKFParticle.GetDistanceFromVertex(trueKFParticleVertex);
    m_true_daughter_ip_xy[daughter_id] = trueKFParticle.GetDistanceFromVertexXY(trueKFParticleVertex);

    m_true_daughter_pv_x[daughter_id] = truePoint == nullptr ? -99. : truePoint->get_x();
    m_true_daughter_pv_y[daughter_id] = truePoint == nullptr ? -99. : truePoint->get_y();
    m_true_daughter_pv_z[daughter_id] = truePoint == nullptr ? -99. : truePoint->get_z();
  }
}

void KFParticle_truthAndDetTools::fillGeant4Branch(PHG4Particle *particle, int daughter_id)
{
  Float_t pT = sqrt(pow(particle->get_px(), 2) + pow(particle->get_py(), 2));

  m_true_daughter_track_history_PDG_ID[daughter_id].push_back(particle->get_pid());
  m_true_daughter_track_history_PDG_mass[daughter_id].push_back(0);
  m_true_daughter_track_history_px[daughter_id].push_back((Float_t) particle->get_px());
  m_true_daughter_track_history_py[daughter_id].push_back((Float_t) particle->get_py());
  m_true_daughter_track_history_pz[daughter_id].push_back((Float_t) particle->get_pz());
  m_true_daughter_track_history_pE[daughter_id].push_back((Float_t) particle->get_e());
  m_true_daughter_track_history_pT[daughter_id].push_back((Float_t) pT);
}

void KFParticle_truthAndDetTools::fillHepMCBranch(HepMC::GenParticle *particle, int daughter_id)
{
  const HepMC::FourVector &myFourVector = particle->momentum();

  m_true_daughter_track_history_PDG_ID[daughter_id].push_back(particle->pdg_id());
  m_true_daughter_track_history_PDG_mass[daughter_id].push_back((Float_t) particle->generatedMass());
  m_true_daughter_track_history_px[daughter_id].push_back((Float_t) myFourVector.px());
  m_true_daughter_track_history_py[daughter_id].push_back((Float_t) myFourVector.py());
  m_true_daughter_track_history_pz[daughter_id].push_back((Float_t) myFourVector.pz());
  m_true_daughter_track_history_pE[daughter_id].push_back((Float_t) myFourVector.e());
  m_true_daughter_track_history_pT[daughter_id].push_back((Float_t) myFourVector.perp());
}

int KFParticle_truthAndDetTools::getHepMCInfo(PHCompositeNode *topNode, TTree * /*m_tree*/, const KFParticle &daughter, int daughter_id)
{
  // Make dummy particle for null pointers and missing nodes
  HepMC::GenParticle *dummyParticle = new HepMC::GenParticle();
  HepMC::FourVector dummyFourVector(0, 0, 0, 0);
  dummyParticle->set_momentum(dummyFourVector);
  dummyParticle->set_pdg_id(0);
  dummyParticle->set_generated_mass(0.);

  PHNodeIterator nodeIter(topNode);
  PHNode *findNode = dynamic_cast<PHNode *>(nodeIter.findFirst(m_trk_map_node_name_nTuple));
  if (findNode)
  {
    dst_trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trk_map_node_name_nTuple);
  }
  else
  {
    std::cout << "KFParticle truth matching: " << m_trk_map_node_name_nTuple << " does not exist" << std::endl;
  }

  track = getTrack(daughter.Id(), dst_trackmap);
  g4particle = getTruthTrack(track, topNode);

  bool isParticleValid = g4particle == nullptr ? false : true;

  if (!isParticleValid)
  {
    std::cout << "KFParticle truth matching: this track is a ghost" << std::endl;
    fillHepMCBranch(dummyParticle, daughter_id);
    return 0;
  }

  m_geneventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
  if (!m_geneventmap)
  {
    std::cout << "KFParticle truth matching: Missing node PHHepMCGenEventMap" << std::endl;
    std::cout << "You will have no mother information" << std::endl;
    fillHepMCBranch(dummyParticle, daughter_id);
    return 0;
  }

  m_genevt = m_geneventmap->get(1);
  if (!m_genevt)
  {
    std::cout << "KFParticle truth matching: Missing node PHHepMCGenEvent" << std::endl;
    std::cout << "You will have no mother information" << std::endl;
    fillHepMCBranch(dummyParticle, daughter_id);
    return 0;
  }

  // Start by looking for our particle in the Geant record
  // Any decay that Geant4 handles will not be in the HepMC record
  // This can happen if you limit the decay volume in the generator
  if (g4particle->get_parent_id() != 0)
  {
    findNode = dynamic_cast<PHNode *>(nodeIter.findFirst("G4TruthInfo"));
    if (findNode)
    {
      m_truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    }
    else
    {
      std::cout << "KFParticle truth matching: G4TruthInfo does not exist" << std::endl;
    }
    while (g4particle->get_parent_id() != 0)
    {
      g4particle = m_truthinfo->GetParticle(g4particle->get_parent_id());
      fillGeant4Branch(g4particle, daughter_id);
    }
  }

  HepMC::GenEvent *theEvent = m_genevt->getEvent();
  HepMC::GenParticle *prevParticle = nullptr;

  // int forbiddenPDGIDs[] = {21, 22};  //Stop tracing history when we reach quarks, gluons and photons
  int forbiddenPDGIDs[] = {0};  // 20230921 - Request made to have gluon information to see about gluon-splitting

  for (HepMC::GenEvent::particle_const_iterator p = theEvent->particles_begin(); p != theEvent->particles_end(); ++p)
  {
    if (((*p)->barcode() == g4particle->get_barcode()))
    {
      prevParticle = *p;
      while (!prevParticle->is_beam())
      {
        bool breakOut = false;
        for (HepMC::GenVertex::particle_iterator mother = prevParticle->production_vertex()->particles_begin(HepMC::parents);
             mother != prevParticle->production_vertex()->particles_end(HepMC::parents); ++mother)
        {
          if (std::find(std::begin(forbiddenPDGIDs), std::end(forbiddenPDGIDs),
                        abs((*mother)->pdg_id())) != std::end(forbiddenPDGIDs))
          {
            breakOut = true;
            break;
          }

          fillHepMCBranch((*mother), daughter_id);
          prevParticle = *mother;
        }
        if (breakOut)
        {
          break;
        }
      }
    }
  }

  return 0;
}  // End of function

void KFParticle_truthAndDetTools::initializeCaloBranches(TTree *m_tree, int daughter_id, const std::string &daughter_number)
{
  m_tree->Branch(TString(daughter_number) + "_EMCAL_DeltaPhi", &detector_emcal_deltaphi[daughter_id], TString(daughter_number) + "_EMCAL_DeltaPhi/F");
  m_tree->Branch(TString(daughter_number) + "_EMCAL_DeltaEta", &detector_emcal_deltaeta[daughter_id], TString(daughter_number) + "_EMCAL_DeltaEta/F");
  m_tree->Branch(TString(daughter_number) + "_EMCAL_energy_3x3", &detector_emcal_energy_3x3[daughter_id], TString(daughter_number) + "_EMCAL_energy_3x3/F");
  m_tree->Branch(TString(daughter_number) + "_EMCAL_energy_5x5", &detector_emcal_energy_5x5[daughter_id], TString(daughter_number) + "_EMCAL_energy_5x5/F");
  m_tree->Branch(TString(daughter_number) + "_EMCAL_energy_cluster", &detector_emcal_cluster_energy[daughter_id], TString(daughter_number) + "_EMCAL_energy_cluster/F");
  m_tree->Branch(TString(daughter_number) + "_IHCAL_DeltaPhi", &detector_ihcal_deltaphi[daughter_id], TString(daughter_number) + "_IHCAL_DeltaPhi/F");
  m_tree->Branch(TString(daughter_number) + "_IHCAL_DeltaEta", &detector_ihcal_deltaeta[daughter_id], TString(daughter_number) + "_IHCAL_DeltaEta/F");
  m_tree->Branch(TString(daughter_number) + "_IHCAL_energy_3x3", &detector_ihcal_energy_3x3[daughter_id], TString(daughter_number) + "_IHCAL_energy_3x3/F");
  m_tree->Branch(TString(daughter_number) + "_IHCAL_energy_5x5", &detector_ihcal_energy_5x5[daughter_id], TString(daughter_number) + "_IHCAL_energy_5x5/F");
  m_tree->Branch(TString(daughter_number) + "_IHCAL_energy_cluster", &detector_ihcal_cluster_energy[daughter_id], TString(daughter_number) + "_IHCAL_energy_cluster/F");
  m_tree->Branch(TString(daughter_number) + "_OHCAL_DeltaPhi", &detector_ohcal_deltaphi[daughter_id], TString(daughter_number) + "_OHCAL_DeltaEta/F");
  m_tree->Branch(TString(daughter_number) + "_OHCAL_DeltaEta", &detector_ohcal_deltaeta[daughter_id], TString(daughter_number) + "_OHCAL_DeltaEta/F");
  m_tree->Branch(TString(daughter_number) + "_OHCAL_energy_3x3", &detector_ohcal_energy_3x3[daughter_id], TString(daughter_number) + "_OHCAL_energy_3x3/F");
  m_tree->Branch(TString(daughter_number) + "_OHCAL_energy_5x5", &detector_ohcal_energy_5x5[daughter_id], TString(daughter_number) + "_OHCAL_energy_5x5/F");
  m_tree->Branch(TString(daughter_number) + "_OHCAL_energy_cluster", &detector_ohcal_cluster_energy[daughter_id], TString(daughter_number) + "_OHCAL_energy_cluster/F");


}

// double KFParticle_truthAndDetTools::Get_CAL_e3x3(
//   SvtxTrack* track, 
//   RawTowerContainer* _towersRawOH, 
//   RawTowerGeomContainer* _geomOH, 
//   int what, 
//   double Zvtx, 
//   double &dphi, 
//   double &deta) {

// double e3x3 = 0.;
// double pathlength = 999.;
// // vector<double> proj;
// // for (SvtxTrack::StateIter stateiter = track->begin_states(); stateiter != track->end_states(); ++stateiter)
// // {
// //   SvtxTrackState *trackstate = stateiter->second;
// //   if(trackstate) { proj.push_back(trackstate->get_pathlength()); }
// // }
// // if(what==0)      { pathlength = proj[proj.size()-3]; } // CEMC
// // else if(what==1) { pathlength = proj[proj.size()-2]; } // HCALIN
// // else if(what==2) { pathlength = proj[proj.size()-1]; } // HCALOUT
// // else { dphi = 9999.; deta = 9999.; return e3x3;}

//   double track_eta = 999.;
//   double track_phi = 999.;
//   SvtxTrackState* trackstate = track->get_state(pathlength);
//   if(trackstate) {
//     double track_x = trackstate->get_x();
//     double track_y = trackstate->get_y();
//     double track_z = trackstate->get_z() - Zvtx;
//     double track_r = sqrt(track_x*track_x+track_y*track_y);
//     track_eta = asinh( track_z / track_r );
//     track_phi = atan2( track_y, track_x );
//   } else { cout << "track state not found!" << endl; dphi = 9999.; deta = 9999.; return e3x3; }

//   double dist = 9999.;
//   RawTower* thetower = nullptr;
//   RawTowerContainer::ConstRange begin_end_rawOH = _towersRawOH->getTowers();
//   for (RawTowerContainer::ConstIterator rtiter = begin_end_rawOH.first; rtiter != begin_end_rawOH.second; ++rtiter) {
//     RawTower *tower = rtiter->second;
//     RawTowerGeom *tower_geom = _geomOH->get_tower_geometry(tower->get_key());
//     //double tower_phi  = tower_geom->get_phi();
//     double tower_x = tower_geom->get_center_x();
//     double tower_y = tower_geom->get_center_y();
//     double tower_z = tower_geom->get_center_z() - Zvtx; // correct for event vertex
//     double tower_r = sqrt(pow(tower_x,2)+pow(tower_y,2));
//     double tower_eta = asinh( tower_z / tower_r );
//     double tower_phi = atan2( tower_y , tower_x );
//     double tmpdist = sqrt(pow(track_eta-tower_eta,2)+pow(track_phi-tower_phi,2));
//     if(tmpdist<dist) { dist = tmpdist; thetower = tower; deta = fabs(track_eta-tower_eta); dphi = fabs(track_phi-tower_phi); }
//   }
//   //cout << "dist: " << dist << " " << deta << " " << dphi << endl;

//   if(!thetower) { dphi = 9999.; deta = 9999.; return e3x3; }
//   RawTowerGeom *thetower_geom = _geomOH->get_tower_geometry(thetower->get_key());
//   unsigned int ieta = thetower_geom->get_bineta();
//   unsigned int jphi = thetower_geom->get_binphi();

//   unsigned int maxbinphi = 63; if(what==0) maxbinphi = 255;
//   unsigned int maxbineta = 23; if(what==0) maxbineta = 93;
//     // calculate e3x3
//     for(unsigned int i=0; i<=2; i++) {
//       for(unsigned int j=0; j<=2; j++) {
//         unsigned int itmp = ieta-1+i;
//         unsigned int jtmp = 0;
//         if(jphi==0 && j==0) { jtmp = maxbinphi; }      // wrap around
//         else if(jphi==maxbinphi && j==2) { jtmp = 0; } // wrap around
//         else { jtmp = jphi-1+j; }
//         if(itmp>=0 && itmp<=maxbineta) { 
//           RawTower* tmptower = _towersRawOH->getTower(itmp,jtmp);
//           if(tmptower) { e3x3 += tmptower->get_energy(); }
//         }
//       }
//     }

//   return e3x3;
// }




float KFParticle_truthAndDetTools::get_e3x3(RawCluster *cluster, RawTowerContainer *Towers, int layer){
  
  std::vector<float> _emcal_tower_phi;
  std::vector<float> _emcal_tower_eta;
  std::vector<float> _emcal_tower_e;
  std::vector<int> _emcal_tower_phi_bin;
  std::vector<int> _emcal_tower_eta_bin;
  float energy;
  float tmp = 0;
  int index = -1;
  int ijk = 0;
  float e3x3 = 0;

  RawCluster::TowerConstRange towers = cluster->get_towers();
  RawCluster::TowerConstIterator toweriter;

  // TowerInfo *towerInfo = nullptr;

  for (toweriter = towers.first; toweriter != towers.second; ++toweriter)
  {

    energy = toweriter->second;

    // _emcal_tower_cluster_id.push_back(clusIter_EMC->first);
    RawTowerGeom *tower_geom = EMCalGeo->get_tower_geometry(toweriter->first);
    _emcal_tower_phi.push_back(tower_geom->get_phi());
    _emcal_tower_eta.push_back(tower_geom->get_eta());
    _emcal_tower_e.push_back(energy);

    _emcal_tower_phi_bin.push_back(tower_geom->get_bineta());
    _emcal_tower_eta_bin.push_back(tower_geom->get_binphi());

    // unsigned int key = TowerInfoDefs::encode_emcal(tower_geom->get_bineta(), tower_geom->get_binphi());
    // towerInfo = EMCAL_Container->get_tower_at_key(key);
    // _emcal_tower_status.push_back(towerInfo->get_status());

    if(energy > tmp){
      index = ijk;
      tmp = energy;
    }

    ijk++;

  }

  int ieta = _emcal_tower_eta_bin[index];
  int jphi = _emcal_tower_phi_bin[index];

  int maxbinphi = 63; if(layer==0) maxbinphi = 255;
  int maxbineta = 23; if(layer==0) maxbineta = 93;
    // calculate e3x3
    for(int i=0; i<=2; i++) {
      for(int j=0; j<=2; j++) {
        int itmp = ieta-1+i;
        int jtmp = 0;
        if(jphi==0 && j==0) { jtmp = maxbinphi; }      // wrap around
        else if(jphi==maxbinphi && j==2) { jtmp = 0; } // wrap around
        else { jtmp = jphi-1+j; }
        if(itmp>=0 && itmp<=maxbineta) { 
          RawTower* tmptower = Towers->getTower(itmp,jtmp);
          if(tmptower) { e3x3 += tmptower->get_energy(); }
        }
      }
    }

  return e3x3;

}


void KFParticle_truthAndDetTools::fillCaloBranch(PHCompositeNode *topNode,
                                                 TTree * /*m_tree*/, const KFParticle &daughter, int daughter_id)
{
  PHNodeIterator nodeIter(topNode);
  PHNode *findNode = dynamic_cast<PHNode *>(nodeIter.findFirst(m_trk_map_node_name_nTuple));
  if (findNode)
  {
    dst_trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trk_map_node_name_nTuple);
  }
  else
  {
    std::cout << "KFParticle truth matching: " << m_trk_map_node_name_nTuple << " does not exist" << std::endl;
  }

  track = getTrack(daughter.Id(), dst_trackmap);

  // ⭐ Ｓｔａｒｔ Ｖａｌｅｒｉｅ＇ｓ ＷＩＰ ⭐
  if ( !clustersEM ) {
    clustersEM = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_CEMC");
    if (!clustersEM)
    {
      // std::cout << "TrackCaloMatch::process_event : FATAL ERROR, cannot find cluster container " << "CLUSTER_CEMC" << std::endl;
      std::cout << __FILE__ << "::" << __func__ << " : FATAL ERROR, cannot find cluster container " << "CLUSTER_CEMC" << std::endl;
      //return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  if(!EMCalGeo)
  {
    EMCalGeo = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
    if(!EMCalGeo)
    {
      std::cout << __FILE__ << "::" << __func__ << " : FATAL ERROR, cannot find cluster container " << "TOWERGEOM_CEMC" << std::endl;
      //return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  if(!_towersEM) 
  {
    _towersEM = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_CEMC");
    if(!_towersEM) 
    {
      std::cout << __FILE__ << "::" << __func__ << " : FATAL ERROR, cannot find cluster container " << "TOWER_CALIB_CEMC" << std::endl;
      //return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  if ( !clustersIH ) {
    clustersIH = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_HCALIN");
    if (!clustersIH)
    {
      std::cout << __FILE__ << "::" << __func__ << " : FATAL ERROR, cannot find cluster container " << "CLUSTER_HCALIN" << std::endl;
      //return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  if(!IHCalGeo)
  {
    IHCalGeo = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
    if(!IHCalGeo)
    {
      std::cout << __FILE__ << "::" << __func__ << " : FATAL ERROR, cannot find cluster container " << "TOWERGEOM_HCALIN" << std::endl;
      //return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  if(!_towersIH) 
  {
    _towersIH = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALIN");
    if(!_towersIH) 
    {
      std::cout << __FILE__ << "::" << __func__ << " : FATAL ERROR, cannot find cluster container " << "TOWER_CALIB_HCALIN" << std::endl;
      //return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
   if ( !clustersOH ) {
    clustersOH = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_HCALOUT");
    if (!clustersOH)
    {
      std::cout << __FILE__ << "::" << __func__ << " : FATAL ERROR, cannot find cluster container " << "CLUSTER_HCALOUT" << std::endl;
      //return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  if(!OHCalGeo)
  {
    OHCalGeo = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
    if(!OHCalGeo)
    {
      std::cout << __FILE__ << "::" << __func__ << " : FATAL ERROR, cannot find cluster container " << "TOWERGEOM_HCALOUT" << std::endl;
      //return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  if(!_towersOH) 
  {
    _towersOH = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALOUT");
    if(!_towersOH) 
    {
      std::cout << __FILE__ << "::" << __func__ << " : FATAL ERROR, cannot find cluster container " << "TOWER_CALIB_HCALOUT" << std::endl;
      //return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  // Radii for track projections
  double caloRadiusEMCal;
  double caloRadiusIHCal;
  double caloRadiusOHCal;
  caloRadiusEMCal = EMCalGeo->get_radius();
  caloRadiusOHCal = OHCalGeo->get_radius();
  caloRadiusIHCal = IHCalGeo->get_radius();

  // Create variables and containers etc.
  bool is_match;
  int index;
  float radius_scale;
  RawCluster *cluster;


//EMCAL******************************************************
  SvtxTrackState *thisState = nullptr;
  thisState = track->get_state(caloRadiusEMCal);
  float _track_phi_emc = NAN;
  float _track_eta_emc = NAN;
  float _track_x_emc = NAN;
  float _track_y_emc = NAN;
  float _track_z_emc = NAN;

  _track_phi_emc = atan2(thisState->get_y(), thisState->get_x());
  _track_eta_emc = asinh(thisState->get_z()/sqrt(thisState->get_x()*thisState->get_x() + thisState->get_y()*thisState->get_y()));
  _track_x_emc = thisState->get_x();
  _track_y_emc = thisState->get_y();
  _track_z_emc = thisState->get_z();

  // EMCal variables and vectors
  float _emcal_phi = NAN;
  float _emcal_eta = NAN;
  float _emcal_x = NAN;
  float _emcal_y = NAN;
  float _emcal_z = NAN;
  radius_scale = NAN;
  float _emcal_3x3 = NAN;
  //float _emcal_5x5 = NAN;
  float _emcal_clusE = NAN;
  std::vector<float> v_emcal_phi;
  std::vector<float> v_emcal_eta;
  std::vector<float> v_emcal_x;
  std::vector<float> v_emcal_y;
  std::vector<float> v_emcal_z;
  std::vector<float> v_emcal_dphi;
  std::vector<float> v_emcal_deta;
  std::vector<float> v_emcal_dr;
  std::vector<float> v_emcal_3x3;
  //std::vector<float> v_emcal_5x5;
  std::vector<float> v_emcal_clusE;
  
  // Set variables for matching
  is_match = false;
  index = -1;
  
  // Create objects, containers, iterators for clusters
  cluster = nullptr;
  RawClusterContainer::Range begin_end_EMC = clustersEM->getClusters();
  RawClusterContainer::Iterator clusIter_EMC;

  // Loop over the EMCal clusters
  for (clusIter_EMC = begin_end_EMC.first; clusIter_EMC != begin_end_EMC.second; ++clusIter_EMC)
  {
    
    //Minimum energy cut
    cluster = clusIter_EMC->second;
    if(cluster->get_energy() < m_emcal_e_low_cut)
    {
      continue;
    }

    // Get cluster information
    _emcal_phi = atan2(cluster->get_y(), cluster->get_x());
    _emcal_eta = asinh(cluster->get_z()/sqrt(cluster->get_x()*cluster->get_x() + cluster->get_y()*cluster->get_y()));
    _emcal_x = cluster->get_x();
    _emcal_y = cluster->get_y();
    radius_scale = m_emcal_radius_user / sqrt(_emcal_x*_emcal_x+_emcal_y*_emcal_y);
    _emcal_z = radius_scale*cluster->get_z();
    _emcal_3x3 = get_e3x3(cluster, _towersEM, 0); //0 for emcal
    _emcal_clusE = cluster->get_energy();

    //Variables to determine potential matches
    float dphi = PiRange(_track_phi_emc - _emcal_phi);
    float dz = _track_z_emc - _emcal_z;
    float deta = abs(_emcal_eta - _track_eta_emc);
    float dr = sqrt((dphi*dphi + deta*deta));

    // Requirements for a possible match
    if(fabs(dphi)<m_dphi_cut && fabs(dz)<m_dz_cut)
    {
      // Add potential match's information to vectors
      v_emcal_phi.push_back(_emcal_phi);
      v_emcal_eta.push_back(_emcal_eta);
      v_emcal_x.push_back(_emcal_x);
      v_emcal_y.push_back(_emcal_y);
      v_emcal_z.push_back(_emcal_z);
      v_emcal_dphi.push_back(dphi);
      v_emcal_deta.push_back(deta);
      v_emcal_dr.push_back(dr);
      v_emcal_3x3.push_back(_emcal_3x3);
      v_emcal_clusE.push_back(_emcal_clusE);
      
      is_match = true;
    }
  }

 // Find the closest match from all potential matches
  if (is_match == true)
  {
  float tmp = 99999;
  for(long unsigned int i = 0; i < v_emcal_dr.size(); i++){
    if(v_emcal_dr[i] < tmp){
      index = i;
      tmp = v_emcal_dr[i];
    }
  }
  }

  // Print out statements
  if(index != -1){
      std::cout<<"matched tracks!!!"<<std::endl;
      std::cout<<"emcal x = "<<v_emcal_x[index]<<" , y = "<<v_emcal_y[index]<<" , z = "<<v_emcal_z[index]<<" , phi = "<<v_emcal_phi[index]<<" , eta = "<<v_emcal_eta[index]<<std::endl;
      std::cout<<"track projected x = "<<_track_x_emc<<" , y = "<<_track_y_emc<<" , z = "<<_track_z_emc<<" , phi = "<<_track_phi_emc<<" , eta = "<<_track_eta_emc<<std::endl;
      std::cout<<"track px = "<<track->get_px()<<" , py = "<<track->get_py()<<" , pz = "<<track->get_pz()<<" , pt = "<<track->get_pt()<<" , p = "<<track->get_p()<<" , charge = "<<track->get_charge()<<std::endl;

  }

  // Save values to the branches!
  if(index == -1){
  detector_emcal_deltaphi[daughter_id] = NAN;
  detector_emcal_deltaeta[daughter_id] = NAN;
  detector_emcal_energy_3x3[daughter_id] = NAN;
  //detector_emcal_energy_5x5[daughter_id] = NAN;
  detector_emcal_cluster_energy[daughter_id] = NAN;
  }
  else{
  detector_emcal_deltaphi[daughter_id] = v_emcal_phi[index];
  detector_emcal_deltaeta[daughter_id] = v_emcal_eta[index];
  detector_emcal_energy_3x3[daughter_id] = v_emcal_3x3[index];
  //detector_emcal_energy_5x5[daughter_id] = _emcal_5x5;
  detector_emcal_cluster_energy[daughter_id] = v_emcal_clusE[index];
  }

//HCAL*******************************************************

//INNER
  // Track projection
  thisState = nullptr;
  thisState = track->get_state(caloRadiusIHCal);
  float _track_phi_ihc = NAN;
  float _track_eta_ihc = NAN;
  float _track_x_ihc = NAN;
  float _track_y_ihc = NAN;
  float _track_z_ihc = NAN;

  _track_phi_ihc = atan2(thisState->get_y(), thisState->get_x());
  _track_eta_ihc = asinh(thisState->get_z()/sqrt(thisState->get_x()*thisState->get_x() + thisState->get_y()*thisState->get_y()));
  _track_x_ihc = thisState->get_x();
  _track_y_ihc = thisState->get_y();
  _track_z_ihc = thisState->get_z();

  // IHCal variables and vectors
  float _ihcal_phi = NAN;
  float _ihcal_eta = NAN;
  float _ihcal_x = NAN;
  float _ihcal_y = NAN;
  float _ihcal_z = NAN;
  float _ihcal_3x3 = NAN;
  //float _ihcal_5x5 = NAN;
  float _ihcal_clusE = NAN;
  radius_scale = NAN;
  std::vector<float> v_ihcal_phi;
  std::vector<float> v_ihcal_eta;
  std::vector<float> v_ihcal_x;
  std::vector<float> v_ihcal_y;
  std::vector<float> v_ihcal_z;
  std::vector<float> v_ihcal_dphi;
  std::vector<float> v_ihcal_deta;
  std::vector<float> v_ihcal_dr;
  std::vector<float> v_ihcal_3x3;
  //std::vector<float> v_ihcal_5x5;
  std::vector<float> v_ihcal_clusE;
  
  // Reset variables for matching
  is_match = false;
  index = -1;
  
  // Create objects, containers, iterators for clusters
  cluster = nullptr;
  RawClusterContainer::Range begin_end_IHC = clustersIH->getClusters();
  RawClusterContainer::Iterator clusIter_IHC;

  // Loop over the IHCal clusters
  for (clusIter_IHC = begin_end_IHC.first; clusIter_IHC != begin_end_IHC.second; ++clusIter_IHC)
  {
    
    // Minimum energy cut
    cluster = clusIter_IHC->second;
    if(cluster->get_energy() < m_ihcal_e_low_cut)
    {
      continue;
    }

    // Get cluster information
    _ihcal_phi = atan2(cluster->get_y(), cluster->get_x());
    _ihcal_eta = asinh(cluster->get_z()/sqrt(cluster->get_x()*cluster->get_x() + cluster->get_y()*cluster->get_y()));
    _ihcal_x = cluster->get_x();
    _ihcal_y = cluster->get_y();
    radius_scale = m_ihcal_radius_user / sqrt(_ihcal_x*_ihcal_x+_ihcal_y*_ihcal_y);
    _ihcal_z = radius_scale*cluster->get_z();
    _ihcal_3x3 = get_e3x3(cluster, _towersIH, 1); //1 for ihcal
    _ihcal_clusE = cluster->get_energy();

    // Variables to determine potential matches
    float dphi = PiRange(_track_phi_ihc - _ihcal_phi);
    float dz = _track_z_ihc - _ihcal_z;
    float deta = abs(_ihcal_eta - _track_eta_ihc);
    float dr = sqrt((dphi*dphi + deta*deta));

    // Requirements for a possible match
    if(fabs(dphi)<m_dphi_cut && fabs(dz)<m_dz_cut)
    {
      //Add potential match's information to vectors
      v_ihcal_phi.push_back(_ihcal_phi);
      v_ihcal_eta.push_back(_ihcal_eta);
      v_ihcal_x.push_back(_ihcal_x);
      v_ihcal_y.push_back(_ihcal_y);
      v_ihcal_z.push_back(_ihcal_z);
      v_ihcal_dphi.push_back(dphi);
      v_ihcal_deta.push_back(deta);
      v_ihcal_dr.push_back(dr);
      v_ihcal_3x3.push_back(_ihcal_3x3);
      v_ihcal_clusE.push_back(_ihcal_clusE);
      
      is_match = true;
    }
  }

 // Find the closest match from all potential matches
  if (is_match == true)
  {
  float tmp = 99999;
  for(long unsigned int i = 0; i < v_ihcal_dr.size(); i++){
    if(v_ihcal_dr[i] < tmp){
      index = i;
      tmp = v_ihcal_dr[i];
    }
  }
  }

  // Print out statihents
  if(index != -1){
      std::cout<<"matched tracks!!!"<<std::endl;
      std::cout<<"ihcal x = "<<v_ihcal_x[index]<<" , y = "<<v_ihcal_y[index]<<" , z = "<<v_ihcal_z[index]<<" , phi = "<<v_ihcal_phi[index]<<" , eta = "<<v_ihcal_eta[index]<<std::endl;
      std::cout<<"track projected x = "<<_track_x_ihc<<" , y = "<<_track_y_ihc<<" , z = "<<_track_z_ihc<<" , phi = "<<_track_phi_ihc<<" , eta = "<<_track_eta_ihc<<std::endl;
      std::cout<<"track px = "<<track->get_px()<<" , py = "<<track->get_py()<<" , pz = "<<track->get_pz()<<" , pt = "<<track->get_pt()<<" , p = "<<track->get_p()<<" , charge = "<<track->get_charge()<<std::endl;

  }
  
  // Save values to the branches!
  if(index == -1){
  detector_ihcal_deltaphi[daughter_id] = NAN;
  detector_ihcal_deltaeta[daughter_id] = NAN;
  detector_ihcal_energy_3x3[daughter_id] = NAN;
  // detector_ihcal_energy_5x5[daughter_id] = NAN;
  detector_ihcal_cluster_energy[daughter_id] = NAN;
  }
  else{
  detector_ihcal_deltaphi[daughter_id] = v_ihcal_phi[index];
  detector_ihcal_deltaeta[daughter_id] = v_ihcal_eta[index];
  detector_ihcal_energy_3x3[daughter_id] = v_ihcal_3x3[index];
  // detector_ihcal_energy_5x5[daughter_id] = _ihcal_5x5;
  detector_ihcal_cluster_energy[daughter_id] = v_ihcal_clusE[index];
  }





//OUTER
  // Track projection
  thisState = nullptr;
  thisState = track->get_state(caloRadiusOHCal);
  float _track_phi_ohc = NAN;
  float _track_eta_ohc = NAN;
  float _track_x_ohc = NAN;
  float _track_y_ohc = NAN;
  float _track_z_ohc = NAN;

  _track_phi_ohc = atan2(thisState->get_y(), thisState->get_x());
  _track_eta_ohc = asinh(thisState->get_z()/sqrt(thisState->get_x()*thisState->get_x() + thisState->get_y()*thisState->get_y()));
  _track_x_ohc = thisState->get_x();
  _track_y_ohc = thisState->get_y();
  _track_z_ohc = thisState->get_z();

  // OHCal variables and vectors
  float _ohcal_phi = NAN;
  float _ohcal_eta = NAN;
  float _ohcal_x = NAN;
  float _ohcal_y = NAN;
  float _ohcal_z = NAN;
  float _ohcal_3x3 = NAN;
  //float _ohcal_5x5 = NAN;
  float _ohcal_clusE = NAN;
  radius_scale = NAN;
  std::vector<float> v_ohcal_phi;
  std::vector<float> v_ohcal_eta;
  std::vector<float> v_ohcal_x;
  std::vector<float> v_ohcal_y;
  std::vector<float> v_ohcal_z;
  std::vector<float> v_ohcal_dphi;
  std::vector<float> v_ohcal_deta;
  std::vector<float> v_ohcal_dr;
  std::vector<float> v_ohcal_3x3;
  //std::vector<float> v_ohcal_5x5;
  std::vector<float> v_ohcal_clusE;
  
  // Reset variables for matching
  is_match = false;
  index = -1;
  
  // Create objects, containers, iterators for clusters
  cluster = nullptr;
  RawClusterContainer::Range begin_end_OHC = clustersOH->getClusters();
  RawClusterContainer::Iterator clusIter_OHC;

  // Loop over the OHCal clusters
  for (clusIter_OHC = begin_end_OHC.first; clusIter_OHC != begin_end_OHC.second; ++clusIter_OHC)
  {
    
    // Minimum energy cut
    cluster = clusIter_OHC->second;
    if(cluster->get_energy() < m_ohcal_e_low_cut)
    {
      continue;
    }

    // Get cluster information
    _ohcal_phi = atan2(cluster->get_y(), cluster->get_x());
    _ohcal_eta = asinh(cluster->get_z()/sqrt(cluster->get_x()*cluster->get_x() + cluster->get_y()*cluster->get_y()));
    _ohcal_x = cluster->get_x();
    _ohcal_y = cluster->get_y();
    radius_scale = m_ohcal_radius_user / sqrt(_ohcal_x*_ohcal_x+_ohcal_y*_ohcal_y);
    _ohcal_z = radius_scale*cluster->get_z();
    _ohcal_3x3 = get_e3x3(cluster, _towersOH, 2); //2 for ohcal
    _ohcal_clusE = cluster->get_energy();

    // Variables to determine potential matches
    float dphi = PiRange(_track_phi_ohc - _ohcal_phi);
    float dz = _track_z_ohc - _ohcal_z;
    float deta = abs(_ohcal_eta - _track_eta_ohc);
    float dr = sqrt((dphi*dphi + deta*deta));

    // Requirohents for a possible match
    if(fabs(dphi)<m_dphi_cut && fabs(dz)<m_dz_cut)
    {
      //Add potential match's information to vectors
      v_ohcal_phi.push_back(_ohcal_phi);
      v_ohcal_eta.push_back(_ohcal_eta);
      v_ohcal_x.push_back(_ohcal_x);
      v_ohcal_y.push_back(_ohcal_y);
      v_ohcal_z.push_back(_ohcal_z);
      v_ohcal_dphi.push_back(dphi);
      v_ohcal_deta.push_back(deta);
      v_ohcal_dr.push_back(dr);
      v_ohcal_3x3.push_back(_ohcal_3x3);
      v_ohcal_clusE.push_back(_ohcal_clusE);
      
      is_match = true;
    }
  }

 // Find the closest match from all potential matches
  if (is_match == true)
  {
  float tmp = 99999;
  for(long unsigned int i = 0; i < v_ohcal_dr.size(); i++){
    if(v_ohcal_dr[i] < tmp){
      index = i;
      tmp = v_ohcal_dr[i];
    }
  }
  }

  // Print out statohents
  if(index != -1){
      std::cout<<"matched tracks!!!"<<std::endl;
      std::cout<<"ohcal x = "<<v_ohcal_x[index]<<" , y = "<<v_ohcal_y[index]<<" , z = "<<v_ohcal_z[index]<<" , phi = "<<v_ohcal_phi[index]<<" , eta = "<<v_ohcal_eta[index]<<std::endl;
      std::cout<<"track projected x = "<<_track_x_ohc<<" , y = "<<_track_y_ohc<<" , z = "<<_track_z_ohc<<" , phi = "<<_track_phi_ohc<<" , eta = "<<_track_eta_ohc<<std::endl;
      std::cout<<"track px = "<<track->get_px()<<" , py = "<<track->get_py()<<" , pz = "<<track->get_pz()<<" , pt = "<<track->get_pt()<<" , p = "<<track->get_p()<<" , charge = "<<track->get_charge()<<std::endl;

  }

  // Save values to the branches!
  if(index == -1){
  detector_ohcal_deltaphi[daughter_id] = NAN;
  detector_ohcal_deltaeta[daughter_id] = NAN;
  detector_ohcal_energy_3x3[daughter_id] = NAN;
  //detector_ohcal_energy_5x5[daughter_id] = NAN;
  detector_ohcal_cluster_energy[daughter_id] = NAN;
  }
  else{
  detector_ohcal_deltaphi[daughter_id] = v_ohcal_phi[index];
  detector_ohcal_deltaeta[daughter_id] = v_ohcal_eta[index];
  detector_ohcal_energy_3x3[daughter_id] = v_ohcal_3x3[index];
  //detector_ohcal_energy_5x5[daughter_id] = _ohcal_5x5;
  detector_ohcal_cluster_energy[daughter_id] = v_ohcal_clusE[index];
  }

  // ⭐ Ｅｎｄ Ｖａｌｅｒｉｅ＇ｓ ＷＩＰ ⭐

}

void KFParticle_truthAndDetTools::initializeDetectorBranches(TTree *m_tree, int daughter_id, const std::string &daughter_number)
{
  m_tree->Branch(TString(daughter_number) + "_residual_x", &residual_x[daughter_id]);
  m_tree->Branch(TString(daughter_number) + "_residual_y", &residual_y[daughter_id]);
  m_tree->Branch(TString(daughter_number) + "_residual_z", &residual_z[daughter_id]);
  m_tree->Branch(TString(daughter_number) + "_layer", &detector_layer[daughter_id]);

  for (auto const &subdetector : Use)
  {
    if (subdetector.second)
    {
      initializeSubDetectorBranches(m_tree, subdetector.first, daughter_id, daughter_number);
    }
  }
}

void KFParticle_truthAndDetTools::initializeSubDetectorBranches(TTree *m_tree, const std::string &detectorName, int daughter_id, const std::string &daughter_number)
{
  if (detectorName == "MVTX")
  {
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_staveID", &mvtx_staveID[daughter_id]);
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_chipID", &mvtx_chipID[daughter_id]);
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_nHits", &detector_nHits_MVTX[daughter_id]);
  }
  if (detectorName == "INTT")
  {
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_ladderZID", &intt_ladderZID[daughter_id]);
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_ladderPhiID", &intt_ladderPhiID[daughter_id]);
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_nHits", &detector_nHits_INTT[daughter_id]);
  }
  if (detectorName == "TPC")
  {
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_sectorID", &tpc_sectorID[daughter_id]);
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_side", &tpc_side[daughter_id]);
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_nHits", &detector_nHits_TPC[daughter_id]);
  }
  if (detectorName == "TPOT")
  {
    m_tree->Branch(TString(daughter_number) + "_" + TString(detectorName) + "_nHits", &detector_nHits_TPOT[daughter_id]);
  }
}

void KFParticle_truthAndDetTools::fillDetectorBranch(PHCompositeNode *topNode,
                                                     TTree * /*m_tree*/, const KFParticle &daughter, int daughter_id)
{
  PHNodeIterator nodeIter(topNode);

  PHNode *findNode = dynamic_cast<PHNode *>(nodeIter.findFirst(m_trk_map_node_name_nTuple));
  if (findNode)
  {
    dst_trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trk_map_node_name_nTuple);
  }
  else
  {
    std::cout << "KFParticle truth matching: " << m_trk_map_node_name_nTuple << " does not exist" << std::endl;
  }

  std::string geoName = "ActsGeometry";
  findNode = dynamic_cast<PHNode *>(nodeIter.findFirst(geoName));
  if (findNode)
  {
    geometry = findNode::getClass<ActsGeometry>(topNode, geoName);
  }
  else
  {
    std::cout << "KFParticle detector info: " << geoName << " does not exist" << std::endl;
  }

  findNode = dynamic_cast<PHNode *>(nodeIter.findFirst("TRKR_CLUSTER"));
  if (findNode)
  {
    dst_clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  }
  else
  {
    std::cout << "KFParticle detector info: TRKR_CLUSTER does not exist" << std::endl;
  }

  track = getTrack(daughter.Id(), dst_trackmap);
  detector_nHits_MVTX[daughter_id] = 0;
  detector_nHits_INTT[daughter_id] = 0;
  detector_nHits_TPC[daughter_id] = 0;
  detector_nHits_TPOT[daughter_id] = 0;

  TrackSeed *silseed = track->get_silicon_seed();
  TrackSeed *tpcseed = track->get_tpc_seed();

  if (silseed)
  {
    for (auto cluster_iter = silseed->begin_cluster_keys(); cluster_iter != silseed->end_cluster_keys(); ++cluster_iter)
    {
      const auto &cluster_key = *cluster_iter;
      const auto trackerID = TrkrDefs::getTrkrId(cluster_key);
  
      detector_layer[daughter_id].push_back(TrkrDefs::getLayer(cluster_key));
  
      unsigned int staveId, chipId, ladderZId, ladderPhiId, sectorId, side;
      staveId = chipId = ladderZId = ladderPhiId = sectorId = side = std::numeric_limits<unsigned int>::quiet_NaN();
  
      if (Use["MVTX"] && trackerID == TrkrDefs::mvtxId)
      {
        staveId = MvtxDefs::getStaveId(cluster_key);
        chipId = MvtxDefs::getChipId(cluster_key);
        ++detector_nHits_MVTX[daughter_id];
      }
      else if (Use["INTT"] && trackerID == TrkrDefs::inttId)
      {
        ladderZId = InttDefs::getLadderZId(cluster_key);
        ladderPhiId = InttDefs::getLadderPhiId(cluster_key);
        ++detector_nHits_INTT[daughter_id];
      }
  
      mvtx_staveID[daughter_id].push_back(staveId);
      mvtx_chipID[daughter_id].push_back(chipId);
      intt_ladderZID[daughter_id].push_back(ladderZId);
      intt_ladderPhiID[daughter_id].push_back(ladderPhiId);
      tpc_sectorID[daughter_id].push_back(sectorId);
      tpc_side[daughter_id].push_back(side);
    }
  }

  if (tpcseed)
  {
    for (auto cluster_iter = tpcseed->begin_cluster_keys(); cluster_iter != tpcseed->end_cluster_keys(); ++cluster_iter)
    {
      const auto &cluster_key = *cluster_iter;
      const auto trackerID = TrkrDefs::getTrkrId(cluster_key);

      detector_layer[daughter_id].push_back(TrkrDefs::getLayer(cluster_key));
  
      unsigned int staveId, chipId, ladderZId, ladderPhiId, sectorId, side;
      staveId = chipId = ladderZId = ladderPhiId = sectorId = side = std::numeric_limits<unsigned int>::quiet_NaN();
  
      if (Use["TPC"] && trackerID == TrkrDefs::tpcId)
      {
        sectorId = TpcDefs::getSectorId(cluster_key);
        side = TpcDefs::getSide(cluster_key);
        ++detector_nHits_TPC[daughter_id];
      }
      else if (Use["TPOT"] && trackerID == TrkrDefs::micromegasId)
      {
        ++detector_nHits_TPOT[daughter_id];
      }
  
      mvtx_staveID[daughter_id].push_back(staveId);
      mvtx_chipID[daughter_id].push_back(chipId);
      intt_ladderZID[daughter_id].push_back(ladderZId);
      intt_ladderPhiID[daughter_id].push_back(ladderPhiId);
      tpc_sectorID[daughter_id].push_back(sectorId);
      tpc_side[daughter_id].push_back(side);
    }
  }

  for (auto state_iter = track->begin_states();
       state_iter != track->end_states();
       ++state_iter)
  {
    SvtxTrackState* tstate = state_iter->second;
    if (tstate->get_pathlength() != 0) //The first track state is an extrapolation so has no cluster
    {
      auto stateckey = tstate->get_cluskey();
      TrkrCluster *cluster = dst_clustermap->findCluster(stateckey);
      auto global = geometry->getGlobalPosition(stateckey, cluster);
  
      residual_x[daughter_id].push_back(global.x() - tstate->get_x());
      residual_y[daughter_id].push_back(global.y() - tstate->get_y());
      residual_z[daughter_id].push_back(global.z() - tstate->get_z()); 
    }
  }
}

void KFParticle_truthAndDetTools::allPVInfo(PHCompositeNode *topNode,
                                            TTree * /*m_tree*/,
                                            const KFParticle &motherParticle,
                                            std::vector<KFParticle> daughters,
                                            std::vector<KFParticle> intermediates)
{
  KFParticle_Tools kfpTupleTools;
  std::vector<KFParticle> primaryVertices = kfpTupleTools.makeAllPrimaryVertices(topNode, m_vtx_map_node_name_nTuple);

  for (auto &primaryVertice : primaryVertices)
  {
    allPV_x.push_back(primaryVertice.GetX());
    allPV_y.push_back(primaryVertice.GetY());
    allPV_z.push_back(primaryVertice.GetZ());

    allPV_mother_IP.push_back(motherParticle.GetDistanceFromVertex(primaryVertice));
    allPV_mother_IPchi2.push_back(motherParticle.GetDeviationFromVertex(primaryVertice));

    for (unsigned int j = 0; j < daughters.size(); ++j)
    {
      allPV_daughter_IP[j].push_back(daughters[j].GetDistanceFromVertex(primaryVertice));
      allPV_daughter_IPchi2[j].push_back(daughters[j].GetDeviationFromVertex(primaryVertice));
    }

    for (unsigned int j = 0; j < intermediates.size(); ++j)
    {
      allPV_intermediates_IP[j].push_back(intermediates[j].GetDistanceFromVertex(primaryVertice));
      allPV_intermediates_IPchi2[j].push_back(intermediates[j].GetDeviationFromVertex(primaryVertice));
    }
  }
}

void KFParticle_truthAndDetTools::clearVectors()
{
  for (int i = 0; i < m_num_tracks_nTuple; ++i)
  {
    // Truth vectors
    m_true_daughter_track_history_PDG_ID[i].clear();
    m_true_daughter_track_history_PDG_mass[i].clear();
    m_true_daughter_track_history_px[i].clear();
    m_true_daughter_track_history_py[i].clear();
    m_true_daughter_track_history_pz[i].clear();
    m_true_daughter_track_history_pE[i].clear();
    m_true_daughter_track_history_pT[i].clear();

    // Detector vectors
    residual_x[i].clear();
    residual_y[i].clear();
    residual_z[i].clear();
    detector_layer[i].clear();
    mvtx_staveID[i].clear();
    mvtx_chipID[i].clear();
    intt_ladderZID[i].clear();
    intt_ladderPhiID[i].clear();
    tpc_sectorID[i].clear();
    tpc_side[i].clear();

    // PV vectors
    allPV_daughter_IP[i].clear();
    allPV_daughter_IPchi2[i].clear();
  }

  allPV_x.clear();
  allPV_y.clear();
  allPV_z.clear();
  allPV_z.clear();

  allPV_mother_IP.clear();
  allPV_mother_IPchi2.clear();

  for (int i = 0; i < m_num_intermediate_states_nTuple; ++i)
  {
    allPV_intermediates_IP[i].clear();
    allPV_intermediates_IPchi2[i].clear();
  }
}
