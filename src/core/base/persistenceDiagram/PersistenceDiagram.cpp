#include <PersistenceDiagram.h>

ttk::CriticalType
  ttk::PersistenceDiagram::getNodeType(ttk::ftm::FTMTree_MT *tree,
                                       ttk::ftm::TreeType treeType,
                                       const SimplexId vertexId) const {
  using ttk::ftm::Node;
  using ttk::ftm::TreeType;
  const Node *node = tree->vertex2Node(vertexId);
  int upDegree{};
  int downDegree{};
  if(treeType == TreeType::Join or treeType == TreeType::Contour) {
    upDegree = node->getNumberOfUpSuperArcs();
    downDegree = node->getNumberOfDownSuperArcs();
  } else {
    upDegree = node->getNumberOfDownSuperArcs();
    downDegree = node->getNumberOfUpSuperArcs();
  }
  int degree = upDegree + downDegree;

  // saddle point
  if(degree > 1) {
    if(upDegree > 1)
      return CriticalType::Saddle2;
    else
      return CriticalType::Saddle1;
  }
  // local extremum
  else {
    if(upDegree)
      return CriticalType::Local_minimum;
    else
      return CriticalType::Local_maximum;
  }
}

void ttk::PersistenceDiagram::getValencesFromLink(
  const SimplexId vertexId,
  const std::vector<std::pair<polarity, polarity>> &vlp,
  DynamicTree &link,
  std::vector<polarity> &toPropageMin,
  std::vector<polarity> &toPropageMax,
  std::vector<std::vector<SimplexId>> &saddleCCMin,
  std::vector<std::vector<SimplexId>> &saddleCCMax) const {

  const auto nbCC = link.getNbCC();

  SimplexId downValence = 0, upValence = 0;
  saddleCCMin[vertexId].clear();
  saddleCCMax[vertexId].clear();

  if(nbCC > 2) {
    std::vector<size_t> CCIds;
    CCIds.reserve(nbCC);
    link.retrieveNbCC(CCIds);
    for(size_t i = 0; i < CCIds.size(); i++) {
      const SimplexId neighbor = CCIds[i];
      const polarity isUpper = vlp[neighbor].first;
      if(isUpper) {
        saddleCCMax[vertexId].emplace_back(neighbor);
        upValence++;
      } else {
        saddleCCMin[vertexId].emplace_back(neighbor);
        downValence++;
      }
    }

    if(downValence > 1) {
      toPropageMin[vertexId] = 255;
    } else {
      saddleCCMin[vertexId].clear();
      toPropageMin[vertexId] = 0;
    }
    if(upValence > 1) {
      toPropageMax[vertexId] = 255;
    } else {
      saddleCCMax[vertexId].clear();
      toPropageMax[vertexId] = 0;
    }

  } else { // not a saddle
    toPropageMax[vertexId] = 0;
    toPropageMin[vertexId] = 0;
  }
}

void ttk::PersistenceDiagram::buildVertexLinkByBoundary(
  const SimplexId vertexId, VLBoundaryType &vlbt) const {

  const auto bid = multiresTriangulation_.getVertexBoundaryIndex(vertexId);
  const auto nneigh = multiresTriangulation_.getVertexNeighborNumber(vertexId);
  vlbt[bid].reserve(nneigh);

  for(SimplexId i = 0; i < nneigh; i++) {
    SimplexId n0 = 0;
    multiresTriangulation_.getVertexNeighbor(vertexId, i, n0);
    for(SimplexId j = i + 1; j < nneigh; j++) {
      SimplexId n1 = 0;
      multiresTriangulation_.getVertexNeighbor(vertexId, j, n1);
      if(multiresTriangulation_.areVerticesNeighbors(n0, n1)) {
        vlbt[bid].emplace_back(i, j);
      }
    }
  }
}

void ttk::PersistenceDiagram::updateCriticalType(
  DynamicTree &link,
  std::vector<std::pair<polarity, polarity>> &vlp,
  std::vector<std::pair<SimplexId, SimplexId>> &vl) const {

  std::vector<SimplexId> monotony_changes_list{};

  for(size_t neighborId = 0; neighborId < vlp.size(); neighborId++) {
    const polarity isBroken = vlp[neighborId].second;
    if(isBroken != 0) {
      monotony_changes_list.emplace_back(neighborId);
    }
  }

  // loop on the link
  //   for each edge that shares n0
  //      if only one break and different polarity : remove
  //      else if only one break and same polarity : insert
  //      else : do nothing
  std::vector<std::vector<std::pair<SimplexId, SimplexId>>> edgesToInsertLater(
    vlp.size());
  std::vector<std::vector<std::pair<SimplexId, SimplexId>>> edgesToRemoveLater(
    vlp.size());

  for(size_t e = 0; e < vl.size(); e++) {
    const SimplexId n0 = vl[e].first;
    const SimplexId n1 = vl[e].second;
    const polarity isBroken0 = vlp[n0].second;
    const polarity isBroken1 = vlp[n1].second;

    if(isBroken0 != 0 and isBroken1 == 0) {
      if(vlp[n0].first != vlp[n1].first) {
        edgesToInsertLater[n0].emplace_back(n0, n1);
      } else {
        edgesToRemoveLater[n0].emplace_back(n0, n1);
      }
    } else if(isBroken0 == 0 and isBroken1 != 0) {
      if(vlp[n0].first != vlp[n1].first) {
        edgesToInsertLater[n1].emplace_back(n1, n0);
      } else {
        edgesToRemoveLater[n1].emplace_back(n1, n0);
      }
    }
  }

  // REMOVE EDGES:
  for(const auto brokenNode : monotony_changes_list) {
    vlp[brokenNode].first = ~vlp[brokenNode].first;
    vlp[brokenNode].second = 0;

    link.getNode(brokenNode)->evert();
    for(const auto &edge : edgesToRemoveLater[brokenNode]) {
      link.removeEdge(edge.first, edge.second);
    }
  }
  if(!edgesToRemoveLater.empty()) {
    // reconnect link
    for(const auto &edge : vl) {
      if(vlp[edge.first].first == vlp[edge.second].first) {
        link.insertEdge(edge.first, edge.second);
      }
    }
  }

  // INSERT EDGES
  for(const auto brokenNode : monotony_changes_list) {
    for(const auto &edge : edgesToInsertLater[brokenNode]) {
      link.insertEdge(edge.first, edge.second);
    }
  }
}

void ttk::PersistenceDiagram::getTripletsFromSaddles(
  const SimplexId vertexId,
  std::vector<triplet> &triplets,
  const std::vector<std::vector<SimplexId>> &vertexReps) const {

  const auto &reps = vertexReps[vertexId];
  const SimplexId m = reps[0];
#ifndef TTK_ENABLE_KAMIKAZE
  const auto &repsm = vertexReps[m];
  if(m == -1 || repsm.empty() || repsm[0] != m) {
    std::cout << "HERE PROBLEM" << std::endl;
  }
#endif // TTK_ENABLE_KAMIKAZE
  for(size_t i = 1; i < reps.size(); i++) {
    const SimplexId n = reps[i];
#ifndef TTK_ENABLE_KAMIKAZE
    const auto &repsn = vertexReps[n];
    if(n == -1 || repsn.empty() || repsn[0] != n) {
      std::cout << "HERE2 PROBLEM" << std::endl;
    }
#endif // TTK_ENABLE_KAMIKAZE
    triplets.emplace_back(vertexId, m, n);
  }
}

double ttk::PersistenceDiagram::predictNextIterationDuration(
  const double currItDuration, const size_t nCurrPairs) const {

  // number of vertices at current iteration
  const double nCurrVerts = multiresTriangulation_.getDecimatedVertexNumber();
  // prediction of duration at iteration n + 1 from iteration n
  // (linear regression, R^2 = 0.994)
  return -0.21 + 0.77 / (decimationLevel_ + 1) - 4.0 * nCurrPairs / nCurrVerts
         + currItDuration
             * (3.3 - 2.32 / nCurrPairs + 32.3 * nCurrPairs / nCurrVerts);
}

void ttk::PersistenceDiagram::stopComputationIf(const bool b) {
  if(b) {
    if(this->decimationLevel_ > this->stoppingDecimationLevel_) {
      std::cout << "Computation stopped at decimation level "
                << this->decimationLevel_ << std::endl;
    }
    this->stoppingDecimationLevel_ = this->decimationLevel_;
  }
}

void ttk::PersistenceDiagram::clearResumableState() {
  // force de-allocation
  vertexRepresentativesMin_ = {};
  vertexRepresentativesMax_ = {};
  vertexLinkPolarity_ = {};
  isNew_ = {};
  vertexLink_ = {};
  link_ = {};
  toProcess_ = {};
  toReprocess_ = {};
  saddleCCMin_ = {};
  saddleCCMax_ = {};
}
