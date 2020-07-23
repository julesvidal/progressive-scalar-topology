// #ifndef SCALARFIELDCRITICALPOINTS_INL
// #define SCALARFIELDCRITICALPOINTS_INL

#include <ScalarFieldCriticalPoints.h>

using std::cout;
using std::endl;
using std::make_pair;
using std::pair;
using std::string;
using std::vector;

template <class dataType>
ttk::ScalarFieldCriticalPoints<dataType>::ScalarFieldCriticalPoints() {

  dimension_ = 0;
  vertexNumber_ = 0;
  scalarValues_ = NULL;
  vertexLinkEdgeLists_ = NULL;
  criticalPoints_ = NULL;
  sosOffsets_ = NULL;
  triangulation_ = NULL;
  integrateTrajectories_ = false;

  forceNonManifoldCheck = false;

  //   threadNumber_ = 1;
}

template <class dataType>
ttk::ScalarFieldCriticalPoints<dataType>::~ScalarFieldCriticalPoints() {
}

template <class dataType>
int ttk::ScalarFieldCriticalPoints<dataType>::execute() {

  // check the consistency of the variables -- to adapt
#ifndef TTK_ENABLE_KAMIKAZE
  if((!dimension_) && ((!triangulation_) || (triangulation_->isEmpty())))
    return -1;
  if((!vertexNumber_) && ((!triangulation_) || (triangulation_->isEmpty())))
    return -2;
  if(!scalarValues_)
    return -3;
  if((!vertexLinkEdgeLists_)
     && ((!triangulation_) || (triangulation_->isEmpty())))
    return -4;
  if(!criticalPoints_)
    return -5;
#endif
  SimplexId decimatedVertexNumber = 0;
  if(triangulation_) {
    if(useMultiresTriangulation_) {
      cout << "Using the multires triangulation. Decimation : "
           << decimationLevel_ << endl;
      multiresTriangulation_.setTriangulation(triangulation_);
      multiresTriangulation_.setDebugLevel(debugLevel_);
      multiresTriangulation_.setDecimationLevel(decimationLevel_);
      decimatedVertexNumber = multiresTriangulation_.getDecimatedVertexNumber();
    }
    vertexNumber_ = triangulation_->getNumberOfVertices();
    dimension_ = triangulation_->getCellVertexNumber(0) - 1;
  }

  if(!sosOffsets_) {
    // let's use our own local copy
    sosOffsets_ = &localSosOffSets_;
  }
  if((SimplexId)sosOffsets_->size() != vertexNumber_) {
    Timer preProcess;
    sosOffsets_->resize(vertexNumber_);
    for(SimplexId i = 0; i < vertexNumber_; i++)
      (*sosOffsets_)[i] = i;

    {
      std::stringstream msg;
      msg << "[ScalarFieldCriticalPoints] Offset pre-processing done in "
          << preProcess.getElapsedTime() << " s. Go!" << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
  }

  Timer t;

  std::vector<char> vertexTypes(
    vertexNumber_, static_cast<char>(CriticalType::Regular));

  if(triangulation_) {
    if(useMultiresTriangulation_) {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
      for(SimplexId i = 0; i < (SimplexId)decimatedVertexNumber; i++) {
        // cout << "looking into " << i << "  (global : "
        //      << multiresTriangulation_.localToGlobalVertexId(i) << " over "
        //      << vertexNumber_ << endl;
        SimplexId globalId = multiresTriangulation_.localToGlobalVertexId(i);
        vertexTypes[globalId] = getCriticalType(globalId, triangulation_);
      }
    } else {
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
      for(SimplexId i = 0; i < (SimplexId)vertexNumber_; i++) {
        vertexTypes[i] = getCriticalType(i, triangulation_);
      }
    }

  } else if(vertexLinkEdgeLists_) {
    // legacy implementation
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(SimplexId i = 0; i < (SimplexId)vertexNumber_; i++) {

      vertexTypes[i] = getCriticalType(i, (*vertexLinkEdgeLists_)[i]);
    }
  }

  SimplexId minimumNumber = 0, maximumNumber = 0, saddleNumber = 0,
            oneSaddleNumber = 0, twoSaddleNumber = 0, monkeySaddleNumber = 0;

  // debug msg
  if(debugLevel_ >= Debug::infoMsg) {
    if(dimension_ == 3) {
      for(SimplexId i = 0; i < vertexNumber_; i++) {
        switch(vertexTypes[i]) {

          case static_cast<char>(CriticalType::Local_minimum):
            minimumNumber++;
            break;

          case static_cast<char>(CriticalType::Saddle1):
            oneSaddleNumber++;
            break;

          case static_cast<char>(CriticalType::Saddle2):
            twoSaddleNumber++;
            break;

          case static_cast<char>(CriticalType::Local_maximum):
            maximumNumber++;
            break;

          case static_cast<char>(CriticalType::Degenerate):
            monkeySaddleNumber++;
            break;
        }
      }
    } else if(dimension_ == 2) {
      for(SimplexId i = 0; i < vertexNumber_; i++) {
        switch(vertexTypes[i]) {

          case static_cast<char>(CriticalType::Local_minimum):
            minimumNumber++;
            break;

          case static_cast<char>(CriticalType::Saddle1):
            saddleNumber++;
            break;

          case static_cast<char>(CriticalType::Local_maximum):
            maximumNumber++;
            break;

          case static_cast<char>(CriticalType::Degenerate):
            monkeySaddleNumber++;
            break;
        }
      }
    }

    {
      std::stringstream msg;
      msg << "[ScalarFieldCriticalPoints] " << minimumNumber << " minima."
          << std::endl;
      if(dimension_ == 3) {
        msg << "[ScalarFieldCriticalPoints] " << oneSaddleNumber
            << " 1-saddle(s)." << std::endl;
        msg << "[ScalarFieldCriticalPoints] " << twoSaddleNumber
            << " 2-saddle(s)." << std::endl;
      }
      if(dimension_ == 2) {
        msg << "[ScalarFieldCriticalPoints] " << saddleNumber << " saddle(s)."
            << std::endl;
      }
      msg << "[ScalarFieldCriticalPoints] " << monkeySaddleNumber
          << " multi-saddle(s)." << std::endl;
      msg << "[ScalarFieldCriticalPoints] " << maximumNumber << " maxima."
          << std::endl;

      //       msg << "[ScalarFieldCriticalPoints] Euler characteristic
      // approximation:";
      //       if(monkeySaddleNumber){
      //         msg << " approximation";
      //       }
      //       msg << ": ";
      //       if(dimension_ == 3){
      //         msg
      //           << minimumNumber - oneSaddleNumber + twoSaddleNumber -
      // maximumNumber;
      //       }
      //       if(dimension_ == 2){
      //         msg << minimumNumber - saddleNumber + maximumNumber;
      //       }
      //       msg << std::endl;
      dMsg(std::cout, msg.str(), Debug::infoMsg);
    }
  }

  // prepare the output
  criticalPoints_->clear();
  criticalPoints_->reserve(vertexNumber_);
  for(SimplexId i = 0; i < vertexNumber_; i++) {
    if(vertexTypes[i] != static_cast<char>(CriticalType::Regular)) {
      criticalPoints_->emplace_back(i, vertexTypes[i]);
    }
  }

  {
    std::stringstream msg;
    msg << "[ScalarFieldCriticalPoints] Data-set (" << vertexNumber_
        << " vertices) processed in " << t.getElapsedTime() << " s. ("
        << threadNumber_ << " thread(s))." << std::endl;
    dMsg(std::cout, msg.str(), 2);
  }

  return 0;
}

template <class dataType>
int ttk::ScalarFieldCriticalPoints<dataType>::executeProgressive() {

  cout << "EXECUTING PROGRESSIVE" << endl;

#ifdef TTK_ENABLE_DYNAMIC_TREES
  cout << " Using Dynamic Trees" << endl;
#else
  cout << " Using Union Finds" << endl;
#endif

#if defined TTK_USE_PRECALCULATED_LINKS_FOR_UF \
  || defined TTK_ENABLE_DYNAMIC_TREES
  cout << " Using precalculated links" << endl;
#ifdef TTK_ENABLE_IMPLICIT_LINK_FOR_MULTIRES_CC
  cout << " with implicit link storage\n" << endl;
#else
  cout << " with explicit storage for each vertex\n" << endl;
#endif
#else
  cout << " No link storage !" << endl;
#endif

  int nb_of_reprocess = 0;
  int nb_of_process = 0;
  int nb_of_types_change = 0;
  int nb_of_reprocess_with_one_monotony_change = 0;
  int nb_of_types_change_with_one_monotony_change = 0;

// check the consistency of the variables -- to adapt
#ifndef TTK_ENABLE_KAMIKAZE
  if((!dimension_) && ((!triangulation_) || (triangulation_->isEmpty())))
    return -1;
  if((!vertexNumber_) && ((!triangulation_) || (triangulation_->isEmpty())))
    return -2;
  if(!scalarValues_)
    return -3;
  if((!vertexLinkEdgeLists_)
     && ((!triangulation_) || (triangulation_->isEmpty())))
    return -4;
  if(!criticalPoints_)
    return -5;
#endif
  SimplexId decimatedVertexNumber = 0;
  if(triangulation_) {
    cout << "Using the multires triangulation. Decimation : "
         << decimationLevel_ << endl;
    multiresTriangulation_.setTriangulation(triangulation_);
    multiresTriangulation_.setDecimationLevel(0);
  }

  if(!sosOffsets_) {
    // let's use our own local copy
    sosOffsets_ = &localSosOffSets_;
  }
  if((SimplexId)sosOffsets_->size() != vertexNumber_) {
    Timer preProcess;
    sosOffsets_->resize(vertexNumber_);
    for(SimplexId i = 0; i < vertexNumber_; i++)
      (*sosOffsets_)[i] = i;

    {
      std::stringstream msg;
      msg << "[ScalarFieldCriticalPoints] Offset pre-processing done in "
          << preProcess.getElapsedTime() << " s. Go!" << std::endl;
      dMsg(std::cout, msg.str(), timeMsg);
    }
  }

  std::vector<char> vertexTypes(
    vertexNumber_, static_cast<char>(CriticalType::Regular));

#ifdef TTK_ENABLE_IMPLICIT_LINK_FOR_MULTIRES_CC
  int dim = multiresTriangulation_.getDimensionality();
  if(dim == 3) {
    vertexLinkByBoundaryType_.resize(27);
  } else if(dim == 2) {
    vertexLinkByBoundaryType_.resize(9);
  } else {
    vertexLinkByBoundaryType_.resize(3);
  }
#else
  vertexLinkExplicit_.resize(vertexNumber_);
#endif
  Timer t;
  double allocation_time = t.getElapsedTime();
  vertexLinkPolarity_.resize(vertexNumber_);
  vertexLink_.resize(vertexNumber_);
  link_.resize(vertexNumber_);
  isNew_.resize(vertexNumber_, true);
  toReprocess_.resize(vertexNumber_, 0);
  toProcess_.resize(vertexNumber_, 0);

  const size_t maxNeigh = dim == 3 ? 14 : (dim == 2 ? 6 : 0);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i = 0; i < vertexNumber_; i++) {
    vertexLinkPolarity_[i].reserve(maxNeigh);
    link_[i].alloc(maxNeigh);
  }

  allocation_time = t.getElapsedTime() - allocation_time;
  cout << "\n Memory allocation/reservation in " << allocation_time << " s."
       << endl;

  startingDecimationLevel_ = decimationLevel_;

  if(triangulation_) {
    double time_link_computation = t.getElapsedTime();
#if defined TTK_ENABLE_DYNAMIC_TREES \
  || defined TTK_USE_PRECALCULATED_LINKS_FOR_UF
#ifdef TTK_ENABLE_IMPLICIT_LINK_FOR_MULTIRES_CC
    // Building link for each representative vertex:
    vector<SimplexId> boundaryRepresentatives;
    multiresTriangulation_.findBoundaryRepresentatives(boundaryRepresentatives);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(unsigned int i = 0; i < boundaryRepresentatives.size(); i++) {
      // cout<<" vertex "<<boundaryRepresentatives[i]<<" of boundary index
      // "<<i<<endl;
      if(boundaryRepresentatives[i] != -1) {
        buildVertexLinkByBoundary(boundaryRepresentatives[i]);
      }
    }
#endif
#endif
    double elapsed = t.getElapsedTime();
    if(debugLevel_ > 3)
      cout << "\n Link computation : " << elapsed - time_link_computation
           << "s. - total " << elapsed << " s.\n"
           << endl;

    multiresTriangulation_.setDecimationLevel(decimationLevel_);
    cout << "decimation level : " << decimationLevel_ << endl;
    decimatedVertexNumber = multiresTriangulation_.getDecimatedVertexNumber();
    vertexNumber_ = triangulation_->getNumberOfVertices();
    dimension_ = triangulation_->getCellVertexNumber(0) - 1;

    double time_initial_loop = t.getElapsedTime();
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(SimplexId i = 0; i < (SimplexId)decimatedVertexNumber; i++) {
      SimplexId globalId = multiresTriangulation_.localToGlobalVertexId(i);
      buildVertexLinkPolarity(globalId);
      toProcess_[globalId] = 255;
      // nb_of_process++;
      char new_type = -1;
#ifdef TTK_ENABLE_DYNAMIC_TREES
      new_type = processVertex(globalId, triangulation_);
#else
#ifdef TTK_ENABLE_IMPLICIT_LINK_FOR_MULTIRES_CC
      associateVertexLinkByBoundary(globalId);
#else
      buildVertexLink(globalId);
#endif
      new_type = processVertexWithUnionFinds(globalId);
#endif
      vertexTypes[globalId] = new_type;
      isNew_[globalId] = 0;
    } // end for openmp

    elapsed = t.getElapsedTime();
    if(debugLevel_ > 3)
      cout << "\n Initial loop : " << elapsed - time_initial_loop
           << "s. - total " << elapsed << " s." << endl;
    Timer while_loop_time;
    while(decimationLevel_ > stoppingDecimationLevel_) {
      decimationLevel_--;
      cout << "decimation level : " << decimationLevel_ << endl;
      multiresTriangulation_.setDecimationLevel(decimationLevel_);
      decimatedVertexNumber = multiresTriangulation_.getDecimatedVertexNumber();

      double time_elapsed = while_loop_time.getElapsedTime();
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
      for(SimplexId i = 0; i < (SimplexId)decimatedVertexNumber; i++) {
        SimplexId globalId = multiresTriangulation_.localToGlobalVertexId(i);
        if(isNew_[globalId] and decimationLevel_ > stoppingDecimationLevel_) {
          buildVertexLinkPolarity(globalId);
        } else if(!isNew_[globalId]) {
          getMonotonyChangeByOldPoint(globalId);
        }
      } // end for openmp
      //
      cout << "\n Time first loop - build polarity and get mono changes "
           << while_loop_time.getElapsedTime() - time_elapsed << endl;
      time_elapsed = while_loop_time.getElapsedTime();
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
      // second Loop on all new points: process needed points
      for(SimplexId i = 0; i < (SimplexId)decimatedVertexNumber; i++) {
        SimplexId globalId = multiresTriangulation_.localToGlobalVertexId(i);
        if(isNew_[globalId] and toProcess_[globalId] != 0) {
          // nb_of_process++;
          char new_type = -1;
#ifdef TTK_ENABLE_DYNAMIC_TREES
          new_type = processVertex(globalId, triangulation_);
#else
#ifdef TTK_ENABLE_IMPLICIT_LINK_FOR_MULTIRES_CC
          associateVertexLinkByBoundary(globalId);
#else
          buildVertexLink(globalId);
#endif
          new_type = processVertexWithUnionFinds(globalId);
#endif
          vertexTypes[globalId] = new_type;
        }
        isNew_[globalId] = false;
      } // end for openmp
      cout << "\n Time second loop - process new points"
           << while_loop_time.getElapsedTime() - time_elapsed << endl;
      time_elapsed = while_loop_time.getElapsedTime();

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
      // Loop on all former points that got impacted by a new point
      for(SimplexId i = 0; i < (SimplexId)decimatedVertexNumber; i++) {
        SimplexId globalId = multiresTriangulation_.localToGlobalVertexId(i);
        if(toReprocess_[globalId] != 0) {
          char new_type = -1;
          if(toProcess_[globalId] != 0) { // this one has been processed once
            // nb_of_reprocess++;
#ifdef TTK_ENABLE_DYNAMIC_TREES
            new_type = reProcess(globalId);
#else
            // updateLinkPolarity(globalId);
            new_type = processVertexWithUnionFinds(globalId);
#endif
          } else { // never been processed
            updateLinkPolarity(globalId);
#ifdef TTK_ENABLE_DYNAMIC_TREES
            new_type = processVertex(globalId, triangulation_);
#else
#ifdef TTK_ENABLE_IMPLICIT_LINK_FOR_MULTIRES_CC
            associateVertexLinkByBoundary(globalId);
#else
            buildVertexLink(globalId);
#endif
            new_type = processVertexWithUnionFinds(globalId);
#endif
            // nb_of_process++;
            toProcess_[globalId] = 255;
          }
          // char old_type = vertexTypes[globalId];
          vertexTypes[globalId] = new_type;
          toReprocess_[globalId] = 0;
        }
      } // end for openmp
      cout << "\n Time third loop - reprocess old points"
           << while_loop_time.getElapsedTime() - time_elapsed << endl;

    } // end while
    cout << "\n Time for while loop " << while_loop_time.getElapsedTime()
         << endl;
  } // end if triangulation_
  SimplexId minimumNumber = 0, maximumNumber = 0, saddleNumber = 0,
            oneSaddleNumber = 0, twoSaddleNumber = 0, monkeySaddleNumber = 0;

  // debug msg
  if(debugLevel_ >= Debug::infoMsg) {
    if(dimension_ == 3) {
      for(SimplexId i = 0; i < vertexNumber_; i++) {
        switch(vertexTypes[i]) {

          case static_cast<char>(CriticalType::Local_minimum):
            minimumNumber++;
            break;

          case static_cast<char>(CriticalType::Saddle1):
            oneSaddleNumber++;
            break;

          case static_cast<char>(CriticalType::Saddle2):
            twoSaddleNumber++;
            break;

          case static_cast<char>(CriticalType::Local_maximum):
            maximumNumber++;
            break;

          case static_cast<char>(CriticalType::Degenerate):
            monkeySaddleNumber++;
            break;
        }
      }
    } else if(dimension_ == 2) {
      for(SimplexId i = 0; i < vertexNumber_; i++) {
        switch(vertexTypes[i]) {

          case static_cast<char>(CriticalType::Local_minimum):
            minimumNumber++;
            break;

          case static_cast<char>(CriticalType::Saddle1):
            saddleNumber++;
            break;

          case static_cast<char>(CriticalType::Local_maximum):
            maximumNumber++;
            break;

          case static_cast<char>(CriticalType::Degenerate):
            monkeySaddleNumber++;
            break;
        }
      }
    }

    {
      std::stringstream msg;
      msg << "[ScalarFieldCriticalPoints] " << minimumNumber << " minima."
          << std::endl;
      if(dimension_ == 3) {
        msg << "[ScalarFieldCriticalPoints] " << oneSaddleNumber
            << " 1-saddle(s)." << std::endl;
        msg << "[ScalarFieldCriticalPoints] " << twoSaddleNumber
            << " 2-saddle(s)." << std::endl;
      }
      if(dimension_ == 2) {
        msg << "[ScalarFieldCriticalPoints] " << saddleNumber << " saddle(s)."
            << std::endl;
      }
      msg << "[ScalarFieldCriticalPoints] " << monkeySaddleNumber
          << " multi-saddle(s)." << std::endl;
      msg << "[ScalarFieldCriticalPoints] " << maximumNumber << " maxima."
          << std::endl;

      //       msg << "[ScalarFieldCriticalPoints] Euler characteristic
      // approximation:";
      //       if(monkeySaddleNumber){
      //         msg << " approximation";
      //       }
      //       msg << ": ";
      //       if(dimension_ == 3){
      //         msg
      //           << minimumNumber - oneSaddleNumber + twoSaddleNumber -
      // maximumNumber;
      //       }
      //       if(dimension_ == 2){
      //         msg << minimumNumber - saddleNumber + maximumNumber;
      //       }
      //       msg << std::endl;
      dMsg(std::cout, msg.str(), Debug::infoMsg);
    }
  }
  // cout<<"postproc of trajectories"<<endl;
  // postProcessTrajectories();
  // cout<<"done"<<endl;
  // prepare the output
  criticalPoints_->clear();
  criticalPoints_->reserve(vertexNumber_);
  criticalGenerationOutput_->clear();
  criticalGenerationOutput_->reserve(vertexNumber_);
  for(SimplexId i = 0; i < vertexNumber_; i++) {
    if(vertexTypes[i] != static_cast<char>(CriticalType::Regular)) {
      criticalPoints_->emplace_back(i, vertexTypes[i]);
      // criticalGenerationOutput_->emplace_back(criticalGeneration_[i]);
    }
  }

  {
    std::stringstream msg;
    msg << "[ScalarFieldCriticalPoints] Data-set (" << vertexNumber_
        << " vertices) processed in " << t.getElapsedTime() - allocation_time
        << " s. (" << threadNumber_ << " thread(s))." << std::endl;
    dMsg(std::cout, msg.str(), 2);
  }
  // std::ofstream file("/tmp/monotony_changes.txt");
  // for(int i=0; i<nbOfMonotonyChanges_.size(); i++){
  //   file<<nbOfMonotonyChanges_[i]<<endl;
  // }
  // file.close();
  cout << " NUMBER OF REPROCESSINGS " << nb_of_reprocess << " ( "
       << nb_of_reprocess_with_one_monotony_change << " with 1 edge change)"
       << endl;
  cout << " NUMBER OF PROCESSINGS " << nb_of_process << endl;
  cout << " NUMBER OF TYPECHANGES " << nb_of_types_change << " ( "
       << nb_of_types_change_with_one_monotony_change << " with 1 edge change)"
       << endl;
  cout << " time processing " << time_processing << endl;
  cout << " time reprocessing " << time_reprocessing << endl;
  cout << " time building link " << time_building_link << endl;
  cout << " time building link polarity " << time_building_link_polarity
       << endl;
  cout << " time get monotony changes " << time_get_monotony_changes << endl;
  multiresTriangulation_.getTimers();
  return 0;
}

template <class dataType>
std::pair<ttk::SimplexId, ttk::SimplexId>
  ttk::ScalarFieldCriticalPoints<dataType>::getNumberOfLowerUpperComponents(
    const SimplexId vertexId, Triangulation *triangulation) {

  SimplexId neighborNumber = triangulation->getVertexNeighborNumber(vertexId);
  std::vector<SimplexId> lowerNeighbors, upperNeighbors;

  for(SimplexId i = 0; i < neighborNumber; i++) {
    SimplexId neighborId = 0;
    triangulation->getVertexNeighbor(vertexId, i, neighborId);

    if(isSosLowerThan((*sosOffsets_)[neighborId], scalarValues_[neighborId],
                      (*sosOffsets_)[vertexId], scalarValues_[vertexId])) {

      lowerNeighbors.push_back(neighborId);
    }

    // upper link
    if(isSosHigherThan((*sosOffsets_)[neighborId], scalarValues_[neighborId],
                       (*sosOffsets_)[vertexId], scalarValues_[vertexId])) {

      upperNeighbors.push_back(neighborId);
    }
  }

  // shortcut, if min or max do not construct the complete star
  if(!forceNonManifoldCheck && lowerNeighbors.empty()) {
    // minimum
    return std::make_pair(0, 1);
  }

  if(!forceNonManifoldCheck && upperNeighbors.empty()) {
    // maximum
    return std::make_pair(1, 0);
  }

  // now do the actual work
  std::vector<UnionFind> lowerSeeds(lowerNeighbors.size());
  std::vector<UnionFind *> lowerList(lowerNeighbors.size());
  std::vector<UnionFind> upperSeeds(upperNeighbors.size());
  std::vector<UnionFind *> upperList(upperNeighbors.size());

  for(SimplexId i = 0; i < (SimplexId)lowerSeeds.size(); i++) {
    lowerList[i] = &(lowerSeeds[i]);
  }
  for(SimplexId i = 0; i < (SimplexId)upperSeeds.size(); i++) {
    upperList[i] = &(upperSeeds[i]);
  }

  SimplexId vertexStarSize = triangulation->getVertexStarNumber(vertexId);

  for(SimplexId i = 0; i < vertexStarSize; i++) {
    SimplexId cellId = 0;
    triangulation->getVertexStar(vertexId, i, cellId);

    SimplexId cellSize = triangulation->getCellVertexNumber(cellId);
    for(SimplexId j = 0; j < cellSize; j++) {
      SimplexId neighborId0 = -1;
      triangulation->getCellVertex(cellId, j, neighborId0);

      if(neighborId0 != vertexId) {
        // we are on the link

        bool lower0 = isSosLowerThan(
          (*sosOffsets_)[neighborId0], scalarValues_[neighborId0],
          (*sosOffsets_)[vertexId], scalarValues_[vertexId]);

        // connect it to everybody except himself and vertexId
        for(SimplexId k = j + 1; k < cellSize; k++) {

          SimplexId neighborId1 = -1;
          triangulation->getCellVertex(cellId, k, neighborId1);

          if((neighborId1 != neighborId0) && (neighborId1 != vertexId)) {

            bool lower1 = isSosLowerThan(
              (*sosOffsets_)[neighborId1], scalarValues_[neighborId1],
              (*sosOffsets_)[vertexId], scalarValues_[vertexId]);

            std::vector<SimplexId> *neighbors = &lowerNeighbors;
            std::vector<UnionFind *> *seeds = &lowerList;

            if(!lower0) {
              neighbors = &upperNeighbors;
              seeds = &upperList;
            }

            if(lower0 == lower1) {
              // connect their union-find sets!
              SimplexId lowerId0 = -1, lowerId1 = -1;
              for(SimplexId l = 0; l < (SimplexId)neighbors->size(); l++) {
                if((*neighbors)[l] == neighborId0) {
                  lowerId0 = l;
                }
                if((*neighbors)[l] == neighborId1) {
                  lowerId1 = l;
                }
              }
              if((lowerId0 != -1) && (lowerId1 != -1)) {
                (*seeds)[lowerId0]
                  = makeUnion((*seeds)[lowerId0], (*seeds)[lowerId1]);
                (*seeds)[lowerId1] = (*seeds)[lowerId0];
              }
            }
          }
        }
      }
    }
  }

  // let's remove duplicates now

  // update the UF if necessary
  for(SimplexId i = 0; i < (SimplexId)lowerList.size(); i++)
    lowerList[i] = lowerList[i]->find();
  for(SimplexId i = 0; i < (SimplexId)upperList.size(); i++)
    upperList[i] = upperList[i]->find();

  std::vector<UnionFind *>::iterator it;
  std::sort(lowerList.begin(), lowerList.end());
  it = unique(lowerList.begin(), lowerList.end());
  lowerList.resize(distance(lowerList.begin(), it));

  std::sort(upperList.begin(), upperList.end());
  it = unique(upperList.begin(), upperList.end());
  upperList.resize(distance(upperList.begin(), it));

  if(debugLevel_ >= Debug::advancedInfoMsg) {
    std::stringstream msg;
    msg << "[ScalarFieldCriticalPoints] Vertex #" << vertexId
        << ": lowerLink-#CC=" << lowerList.size()
        << " upperLink-#CC=" << upperList.size() << std::endl;

    dMsg(std::cout, msg.str(), Debug::advancedInfoMsg);
  }

  return std::make_pair(lowerList.size(), upperList.size());
}

template <class dataType>
char ttk::ScalarFieldCriticalPoints<dataType>::processVertexWithUnionFinds(
  const SimplexId &vertexId) {

  // Timer t;

  for(size_t neighborId = 0; neighborId < vertexLinkPolarity_[vertexId].size();
      neighborId++) {
    polarity isBroken = vertexLinkPolarity_[vertexId][neighborId].second;
    polarity oldPolarity = vertexLinkPolarity_[vertexId][neighborId].first;
    if(isBroken != 0) {
      vertexLinkPolarity_[vertexId][neighborId].first
        = (oldPolarity != 0) ? 0 : 255;
      vertexLinkPolarity_[vertexId][neighborId].second = 0;
    }
  }

  SimplexId downValence, upValence;
#ifdef TTK_USE_PRECALCULATED_LINKS_FOR_UF
  std::tie(downValence, upValence)
    = getNumberOfLowerUpperComponentsMultiresWithUnionFindsAndUsingLinks(
      vertexId);

#else

  std::tie(downValence, upValence)
    = getNumberOfLowerUpperComponentsMultiresWithUnionFinds(vertexId);

#endif

  char c = valencesToType(downValence, upValence);
  // time_processing += t.getElapsedTime();
  return c;
}

template <class dataType>
char ttk::ScalarFieldCriticalPoints<dataType>::processVertex(
  const SimplexId &vertexId, Triangulation * /*triangulation*/) {

  // Timer t;
  SimplexId downValence, upValence;

  std::tie(downValence, upValence) = getNumberOfLowerUpperComponentsMultires(
    vertexId, &multiresTriangulation_);

  char c = valencesToType(downValence, upValence);
  // time_processing += t.getElapsedTime();
  return c;
}

template <class dataType>
char ttk::ScalarFieldCriticalPoints<dataType>::getCriticalType(
  const SimplexId &vertexId, Triangulation *triangulation) {
  // cout << "vertex " << vertexId << endl;
  SimplexId downValence, upValence;
  if(useMultiresTriangulation_) {
    // cout<<"MULTIRES ALGO"<<endl;
    std::vector<int> gridDimension;
    triangulation_->getGridDimensions(gridDimension);
    // cout << "getting number of multi res" << endl;
    if(multiresTriangulation_.isInTriangulation(vertexId)) {
      // cout<<"getting number of components"<<endl;
      std::tie(downValence, upValence)
        = getNumberOfLowerUpperComponentsMultires(
          vertexId, &multiresTriangulation_);
    } else {
      downValence = 1;
      upValence = 1;
    }
  } else {
    std::tie(downValence, upValence)
      = getNumberOfLowerUpperComponents(vertexId, triangulation);
  }

  if(downValence == -1 && upValence == -1) {
    return -1;
  } else if(downValence == 0 && upValence == 1) {
    return static_cast<char>(CriticalType::Local_minimum);
  } else if(downValence == 1 && upValence == 0) {
    return static_cast<char>(CriticalType::Local_maximum);
  } else if(downValence == 1 && upValence == 1) {
    // regular point
    return static_cast<char>(CriticalType::Regular);
  } else {
    // saddles
    if(dimension_ == 2) {
      if((downValence == 2 && upValence == 1)
         || (downValence == 1 && upValence == 2)
         || (downValence == 2 && upValence == 2)) {
        // regular saddle
        return static_cast<char>(CriticalType::Saddle1);
      } else {
        // monkey saddle, saddle + extremum
        return static_cast<char>(CriticalType::Degenerate);
        // NOTE: you may have multi-saddles on the boundary in that
        // configuration
        // to make this computation 100% correct, one would need to
        // disambiguate boundary from interior vertices
      }
    } else if(dimension_ == 3) {
      if(downValence == 2 && upValence == 1) {
        return static_cast<char>(CriticalType::Saddle1);
      } else if(downValence == 1 && upValence == 2) {
        return static_cast<char>(CriticalType::Saddle2);
      } else {
        // monkey saddle, saddle + extremum
        return static_cast<char>(CriticalType::Degenerate);
        // NOTE: we may have a similar effect in 3D (TODO)
      }
    }
  }

  // -2: regular points
  return static_cast<char>(CriticalType::Regular);
}

template <class dataType>
char ttk::ScalarFieldCriticalPoints<dataType>::getCriticalType(
  const SimplexId &vertexId,
  const std::vector<std::pair<SimplexId, SimplexId>> &vertexLink) {

  std::map<SimplexId, SimplexId> global2LowerLink, global2UpperLink;
  std::map<SimplexId, SimplexId>::iterator neighborIt;

  SimplexId lowerCount = 0, upperCount = 0;

  for(SimplexId i = 0; i < (SimplexId)vertexLink.size(); i++) {

    SimplexId neighborId = vertexLink[i].first;

    // first vertex
    // lower link search
    if(isSosLowerThan((*sosOffsets_)[neighborId], scalarValues_[neighborId],
                      (*sosOffsets_)[vertexId], scalarValues_[vertexId])) {

      neighborIt = global2LowerLink.find(neighborId);
      if(neighborIt == global2LowerLink.end()) {
        // not in there, add it
        global2LowerLink[neighborId] = lowerCount;
        lowerCount++;
      }
    }

    // upper link
    if(isSosHigherThan((*sosOffsets_)[neighborId], scalarValues_[neighborId],
                       (*sosOffsets_)[vertexId], scalarValues_[vertexId])) {

      neighborIt = global2UpperLink.find(neighborId);
      if(neighborIt == global2UpperLink.end()) {
        // not in there, add it
        global2UpperLink[neighborId] = upperCount;
        upperCount++;
      }
    }

    // second vertex
    neighborId = vertexLink[i].second;

    // lower link search
    if(isSosLowerThan((*sosOffsets_)[neighborId], scalarValues_[neighborId],
                      (*sosOffsets_)[vertexId], scalarValues_[vertexId])) {

      neighborIt = global2LowerLink.find(neighborId);
      if(neighborIt == global2LowerLink.end()) {
        // not in there, add it
        global2LowerLink[neighborId] = lowerCount;
        lowerCount++;
      }
    }

    // upper link
    if(isSosHigherThan((*sosOffsets_)[neighborId], scalarValues_[neighborId],
                       (*sosOffsets_)[vertexId], scalarValues_[vertexId])) {

      neighborIt = global2UpperLink.find(neighborId);
      if(neighborIt == global2UpperLink.end()) {
        // not in there, add it
        global2UpperLink[neighborId] = upperCount;
        upperCount++;
      }
    }
  }

  if(debugLevel_ >= Debug::advancedInfoMsg) {
    std::stringstream msg;
    msg << "[ScalarFieldCriticalPoints] Vertex #" << vertexId << " lower link ("
        << lowerCount << " vertices)" << std::endl;

    msg << "[ScalarFieldCriticalPoints] Vertex #" << vertexId << " upper link ("
        << upperCount << " vertices)" << std::endl;

    dMsg(std::cout, msg.str(), Debug::advancedInfoMsg);
  }

  if(!lowerCount) {
    // minimum
    return static_cast<char>(CriticalType::Local_minimum);
  }
  if(!upperCount) {
    // maximum
    return static_cast<char>(CriticalType::Local_maximum);
  }

  // so far 40% of the computation, that's ok.

  // now enumerate the connected components of the lower and upper links
  // NOTE: a breadth first search might be faster than a UF
  // if so, one would need the one-skeleton data structure, not the edge list
  std::vector<UnionFind> lowerSeeds(lowerCount);
  std::vector<UnionFind> upperSeeds(upperCount);
  std::vector<UnionFind *> lowerList(lowerCount);
  std::vector<UnionFind *> upperList(upperCount);
  for(SimplexId i = 0; i < (SimplexId)lowerList.size(); i++)
    lowerList[i] = &(lowerSeeds[i]);
  for(SimplexId i = 0; i < (SimplexId)upperList.size(); i++)
    upperList[i] = &(upperSeeds[i]);

  for(SimplexId i = 0; i < (SimplexId)vertexLink.size(); i++) {

    SimplexId neighborId0 = vertexLink[i].first;
    SimplexId neighborId1 = vertexLink[i].second;

    // process the lower link
    if((isSosLowerThan((*sosOffsets_)[neighborId0], scalarValues_[neighborId0],
                       (*sosOffsets_)[vertexId], scalarValues_[vertexId]))
       && (isSosLowerThan((*sosOffsets_)[neighborId1],
                          scalarValues_[neighborId1], (*sosOffsets_)[vertexId],
                          scalarValues_[vertexId]))) {

      // both vertices are lower, let's add that edge and update the UF
      std::map<SimplexId, SimplexId>::iterator n0It
        = global2LowerLink.find(neighborId0);
      std::map<SimplexId, SimplexId>::iterator n1It
        = global2LowerLink.find(neighborId1);

      lowerList[n0It->second]
        = makeUnion(lowerList[n0It->second], lowerList[n1It->second]);
      lowerList[n1It->second] = lowerList[n0It->second];
    }

    // process the upper link
    if((isSosHigherThan((*sosOffsets_)[neighborId0], scalarValues_[neighborId0],
                        (*sosOffsets_)[vertexId], scalarValues_[vertexId]))
       && (isSosHigherThan((*sosOffsets_)[neighborId1],
                           scalarValues_[neighborId1], (*sosOffsets_)[vertexId],
                           scalarValues_[vertexId]))) {

      // both vertices are lower, let's add that edge and update the UF
      std::map<SimplexId, SimplexId>::iterator n0It
        = global2UpperLink.find(neighborId0);
      std::map<SimplexId, SimplexId>::iterator n1It
        = global2UpperLink.find(neighborId1);

      upperList[n0It->second]
        = makeUnion(upperList[n0It->second], upperList[n1It->second]);
      upperList[n1It->second] = upperList[n0It->second];
    }
  }

  // let's remove duplicates
  std::vector<UnionFind *>::iterator it;
  // update the UFs if necessary
  for(SimplexId i = 0; i < (SimplexId)lowerList.size(); i++)
    lowerList[i] = lowerList[i]->find();
  for(SimplexId i = 0; i < (SimplexId)upperList.size(); i++)
    upperList[i] = upperList[i]->find();

  sort(lowerList.begin(), lowerList.end());
  it = unique(lowerList.begin(), lowerList.end());
  lowerList.resize(distance(lowerList.begin(), it));

  sort(upperList.begin(), upperList.end());
  it = unique(upperList.begin(), upperList.end());
  upperList.resize(distance(upperList.begin(), it));

  if(debugLevel_ >= Debug::advancedInfoMsg) {
    std::stringstream msg;
    msg << "[ScalarFieldCriticalPoints] Vertex #" << vertexId
        << ": lowerLink-#CC=" << lowerList.size()
        << " upperLink-#CC=" << upperList.size() << std::endl;

    dMsg(std::cout, msg.str(), Debug::advancedInfoMsg);
  }

  if((lowerList.size() == 1) && (upperList.size() == 1))
    // regular point
    return static_cast<char>(CriticalType::Regular);
  else {
    // saddles
    if(dimension_ == 2) {
      if((lowerList.size() > 2) || (upperList.size() > 2)) {
        // monkey saddle
        return static_cast<char>(CriticalType::Degenerate);
      } else {
        // regular saddle
        return static_cast<char>(CriticalType::Saddle1);
        // NOTE: you may have multi-saddles on the boundary in that
        // configuration
        // to make this computation 100% correct, one would need to disambiguate
        // boundary from interior vertices
      }
    } else if(dimension_ == 3) {
      if((lowerList.size() == 2) && (upperList.size() == 1)) {
        return static_cast<char>(CriticalType::Saddle1);
      } else if((lowerList.size() == 1) && (upperList.size() == 2)) {
        return static_cast<char>(CriticalType::Saddle2);
      } else {
        // monkey saddle
        return static_cast<char>(CriticalType::Degenerate);
        // NOTE: we may have a similar effect in 3D (TODO)
      }
    }
  }

  // -2: regular points
  return static_cast<char>(CriticalType::Regular);
}

template <class dataType>
std::pair<ttk::SimplexId, ttk::SimplexId>
  ttk::ScalarFieldCriticalPoints<dataType>::
    getNumberOfLowerUpperComponentsMultiresWithUnionFindsAndUsingLinks(
      const SimplexId vertexId) {

  SimplexId neighborNumber = vertexLinkPolarity_[vertexId].size();
  // cout << neighborNumber << " neighbors" << endl;
  // multiresTriangulation->printInfos();

  std::vector<SimplexId> lowerNeighbors, upperNeighbors;
  // bool alreadyFilled = false;
  // if(heightOfNeighbors_[vertexId].size() > 0) {
  //   alreadyFilled = true;
  // }
  // cout<<"checking filled "<<alreadyFilled<<endl;
  // vector<bool> localNeighborsHeight;
  for(SimplexId i = 0; i < neighborNumber; i++) {
    if(vertexLinkPolarity_[vertexId][i].first == 0) {
      lowerNeighbors.push_back(i);
    }

    // upper link
    if(vertexLinkPolarity_[vertexId][i].first != 0) {
      upperNeighbors.push_back(i);
    }
  }

  // shortcut, if min or max do not construct the complete star
  if(!forceNonManifoldCheck && lowerNeighbors.empty()) {
    // minimum
    // cout<<"return simple min"<<endl;
    return std::make_pair(0, 1);
  }

  if(!forceNonManifoldCheck && upperNeighbors.empty()) {
    // maximum
    // cout<<"return simple max"<<endl;
    return std::make_pair(1, 0);
  }

  // now do the actual work
  std::vector<UnionFind> lowerSeeds(lowerNeighbors.size());
  std::vector<UnionFind *> lowerList(lowerNeighbors.size());
  std::vector<UnionFind> upperSeeds(upperNeighbors.size());
  std::vector<UnionFind *> upperList(upperNeighbors.size());

  for(SimplexId i = 0; i < (SimplexId)lowerSeeds.size(); i++) {
    lowerList[i] = &(lowerSeeds[i]);
  }
  for(SimplexId i = 0; i < (SimplexId)upperSeeds.size(); i++) {
    upperList[i] = &(upperSeeds[i]);
  }

  for(size_t edgeId = 0; edgeId < vertexLink_[vertexId]->size(); edgeId++) {
    SimplexId localId0 = vertexLink_[vertexId]->at(edgeId).first;
    SimplexId localId1 = vertexLink_[vertexId]->at(edgeId).second;

    polarity isUpper0 = vertexLinkPolarity_[vertexId][localId0].first;
    polarity isUpper1 = vertexLinkPolarity_[vertexId][localId1].first;

    std::vector<SimplexId> *neighbors = &lowerNeighbors;
    std::vector<UnionFind *> *seeds = &lowerList;

    if(isUpper0 != 0) {
      neighbors = &upperNeighbors;
      seeds = &upperList;
    }

    if(isUpper0 == isUpper1) {

      SimplexId lowerId0 = -1, lowerId1 = -1;
      for(SimplexId j = 0; j < (SimplexId)neighbors->size(); j++) {
        if((*neighbors)[j] == localId0) {
          lowerId0 = j;
        }
        if((*neighbors)[j] == localId1) {
          lowerId1 = j;
        }
      }
      if((lowerId0 != -1) && (lowerId1 != -1)) {
        (*seeds)[lowerId0] = makeUnion((*seeds)[lowerId0], (*seeds)[lowerId1]);
        (*seeds)[lowerId1] = (*seeds)[lowerId0];
      } else {
        cout << "SHOULDNT HAPPEN" << endl;
      }
    }
  }
  // for(SimplexId j = 0; j < neighborNumber; j++) {
  //   SimplexId neighborId0;
  //   multiresTriangulation_.getVertexNeighbor(vertexId, j, neighborId0);

  //   bool lower0
  //     = isSosLowerThan((*sosOffsets_)[neighborId0],
  //     scalarValues_[neighborId0],
  //                      (*sosOffsets_)[vertexId], scalarValues_[vertexId]);

  //   for(SimplexId k = j + 1; k < neighborNumber; k++) {
  //     SimplexId neighborId1;
  //     multiresTriangulation_.getVertexNeighbor(vertexId, k, neighborId1);
  //     // cout<<"  neighbor0 "<<j<<"  neighborId0 "<<neighborId0<<endl;
  //     // cout<<scalarValues_[vertexId]<<"  "<<scalarValues_[neighborId0]<<"
  //     lower ? "<<lower0<<endl;
  //     // cout << "  neighbor1 " << k << "  neighborId1 " << neighborId1 <<
  //     endl;

  //     bool lower1 = isSosLowerThan(
  //       (*sosOffsets_)[neighborId1], scalarValues_[neighborId1],
  //       (*sosOffsets_)[vertexId], scalarValues_[vertexId]);
  //   // cout<<scalarValues_[vertexId]<<"  "<<scalarValues_[neighborId1]<<"
  //   lower ? "<<lower1<<endl;

  //     if(lower0 == lower1) {
  //       bool areThoseNeighbors = multiresTriangulation_.areVerticesNeighbors(
  //         neighborId0, neighborId1);
  //       // cout<<"  are those neighbors ? "<<areThoseNeighbors<<endl;
  //       if(areThoseNeighbors) {
  //         SimplexId lowerId0 = -1, lowerId1 = -1;
  //         for(SimplexId l = 0; l < (SimplexId)neighbors->size(); l++) {
  //           if((*neighbors)[l] == neighborId0) {
  //             lowerId0 = l;
  //           }
  //           if((*neighbors)[l] == neighborId1) {
  //             lowerId1 = l;
  //           }
  //         }
  //         // connect their union-find sets!
  //         if((lowerId0 != -1) && (lowerId1 != -1)) {
  //           (*seeds)[lowerId0]
  //             = makeUnion((*seeds)[lowerId0], (*seeds)[lowerId1]);
  //           (*seeds)[lowerId1] = (*seeds)[lowerId0];
  //         } else {
  //           cout << "THIS SOULDNT HAPPEN" << endl;
  //         }
  //       }
  //     }
  //   }
  // }

  // let's remove duplicates now

  // update the UF if necessary
  for(SimplexId i = 0; i < (SimplexId)lowerList.size(); i++)
    lowerList[i] = lowerList[i]->find();
  for(SimplexId i = 0; i < (SimplexId)upperList.size(); i++)
    upperList[i] = upperList[i]->find();

  std::vector<UnionFind *>::iterator it;
  std::sort(lowerList.begin(), lowerList.end());
  it = unique(lowerList.begin(), lowerList.end());
  lowerList.resize(distance(lowerList.begin(), it));

  std::sort(upperList.begin(), upperList.end());
  it = unique(upperList.begin(), upperList.end());
  upperList.resize(distance(upperList.begin(), it));

  if(debugLevel_ >= Debug::advancedInfoMsg) {
    std::stringstream msg;
    msg << "[ScalarFieldCriticalPoints] Vertex #" << vertexId
        << ": lowerLink-#CC=" << lowerList.size()
        << " upperLink-#CC=" << upperList.size() << std::endl;

    dMsg(std::cout, msg.str(), Debug::advancedInfoMsg);
  }

  return std::make_pair(lowerList.size(), upperList.size());
}

template <class dataType>
std::pair<ttk::SimplexId, ttk::SimplexId>
  ttk::ScalarFieldCriticalPoints<dataType>::
    getNumberOfLowerUpperComponentsMultiresWithUnionFinds(
      const SimplexId vertexId) {

  SimplexId neighborNumber
    = multiresTriangulation_.getVertexNeighborNumber(vertexId);
  // cout << neighborNumber << " neighbors" << endl;
  // multiresTriangulation->printInfos();

  std::vector<SimplexId> lowerNeighbors, upperNeighbors;
  // bool alreadyFilled = false;
  // if(heightOfNeighbors_[vertexId].size() > 0) {
  //   alreadyFilled = true;
  // }
  // cout<<"checking filled "<<alreadyFilled<<endl;
  // vector<bool> localNeighborsHeight;
  for(SimplexId i = 0; i < neighborNumber; i++) {
    SimplexId neighborId = -1;
    multiresTriangulation_.getVertexNeighbor(vertexId, i, neighborId);

    if(isSosLowerThan((*sosOffsets_)[neighborId], scalarValues_[neighborId],
                      (*sosOffsets_)[vertexId], scalarValues_[vertexId])) {
      // localNeighborsHeight.push_back(false);
      lowerNeighbors.push_back(neighborId);
    }

    // upper link
    if(isSosHigherThan((*sosOffsets_)[neighborId], scalarValues_[neighborId],
                       (*sosOffsets_)[vertexId], scalarValues_[vertexId])) {
      // localNeighborsHeight.push_back(true);
      upperNeighbors.push_back(neighborId);
    }
  }

  // shortcut, if min or max do not construct the complete star
  if(!forceNonManifoldCheck && lowerNeighbors.empty()) {
    // minimum
    // cout<<"return simple min"<<endl;
    return std::make_pair(0, 1);
  }

  if(!forceNonManifoldCheck && upperNeighbors.empty()) {
    // maximum
    // cout<<"return simple max"<<endl;
    return std::make_pair(1, 0);
  }

  // if(alreadyFilled){
  //   // cout<<"entering 870"<<endl;
  //   // Check is there is a change
  //   bool hasChanged = checkChangeInNeighborsHeight(
  //     localNeighborsHeight, vertexId, neighborNumber);
  //   if(hasChanged == false) {
  //     return std::make_pair(-1,-1);
  //   }
  // }
  // else{
  //   fillInfoOnNeighborsHeight(localNeighborsHeight, vertexId,
  //   neighborNumber);
  // }

  // cout<<"entering 879"<<endl;

  // now do the actual work
  std::vector<UnionFind> lowerSeeds(lowerNeighbors.size());
  std::vector<UnionFind *> lowerList(lowerNeighbors.size());
  std::vector<UnionFind> upperSeeds(upperNeighbors.size());
  std::vector<UnionFind *> upperList(upperNeighbors.size());

  for(SimplexId i = 0; i < (SimplexId)lowerSeeds.size(); i++) {
    lowerList[i] = &(lowerSeeds[i]);
  }
  for(SimplexId i = 0; i < (SimplexId)upperSeeds.size(); i++) {
    upperList[i] = &(upperSeeds[i]);
  }

  for(SimplexId j = 0; j < neighborNumber; j++) {
    SimplexId neighborId0;
    multiresTriangulation_.getVertexNeighbor(vertexId, j, neighborId0);

    bool lower0
      = isSosLowerThan((*sosOffsets_)[neighborId0], scalarValues_[neighborId0],
                       (*sosOffsets_)[vertexId], scalarValues_[vertexId]);

    for(SimplexId k = j + 1; k < neighborNumber; k++) {
      SimplexId neighborId1;
      multiresTriangulation_.getVertexNeighbor(vertexId, k, neighborId1);
      // cout<<"  neighbor0 "<<j<<"  neighborId0 "<<neighborId0<<endl;
      // cout<<scalarValues_[vertexId]<<"  "<<scalarValues_[neighborId0]<<"
      // lower ? "<<lower0<<endl; cout << "  neighbor1 " << k << "  neighborId1
      // " << neighborId1 << endl;

      bool lower1 = isSosLowerThan(
        (*sosOffsets_)[neighborId1], scalarValues_[neighborId1],
        (*sosOffsets_)[vertexId], scalarValues_[vertexId]);
      // cout<<scalarValues_[vertexId]<<"  "<<scalarValues_[neighborId1]<<"
      // lower ? "<<lower1<<endl;

      std::vector<SimplexId> *neighbors = &lowerNeighbors;
      std::vector<UnionFind *> *seeds = &lowerList;

      if(!lower0) {
        neighbors = &upperNeighbors;
        seeds = &upperList;
      }
      if(lower0 == lower1) {
        bool areThoseNeighbors = multiresTriangulation_.areVerticesNeighbors(
          neighborId0, neighborId1);
        // cout<<"  are those neighbors ? "<<areThoseNeighbors<<endl;
        if(areThoseNeighbors) {
          SimplexId lowerId0 = -1, lowerId1 = -1;
          for(SimplexId l = 0; l < (SimplexId)neighbors->size(); l++) {
            if((*neighbors)[l] == neighborId0) {
              lowerId0 = l;
            }
            if((*neighbors)[l] == neighborId1) {
              lowerId1 = l;
            }
          }
          // connect their union-find sets!
          if((lowerId0 != -1) && (lowerId1 != -1)) {
            (*seeds)[lowerId0]
              = makeUnion((*seeds)[lowerId0], (*seeds)[lowerId1]);
            (*seeds)[lowerId1] = (*seeds)[lowerId0];
          } else {
            cout << "THIS SOULDNT HAPPEN" << endl;
          }
        }
      }
    }
  }

  // let's remove duplicates now

  // update the UF if necessary
  for(SimplexId i = 0; i < (SimplexId)lowerList.size(); i++)
    lowerList[i] = lowerList[i]->find();
  for(SimplexId i = 0; i < (SimplexId)upperList.size(); i++)
    upperList[i] = upperList[i]->find();

  std::vector<UnionFind *>::iterator it;
  std::sort(lowerList.begin(), lowerList.end());
  it = unique(lowerList.begin(), lowerList.end());
  lowerList.resize(distance(lowerList.begin(), it));

  std::sort(upperList.begin(), upperList.end());
  it = unique(upperList.begin(), upperList.end());
  upperList.resize(distance(upperList.begin(), it));

  if(debugLevel_ >= Debug::advancedInfoMsg) {
    std::stringstream msg;
    msg << "[ScalarFieldCriticalPoints] Vertex #" << vertexId
        << ": lowerLink-#CC=" << lowerList.size()
        << " upperLink-#CC=" << upperList.size() << std::endl;

    dMsg(std::cout, msg.str(), Debug::advancedInfoMsg);
  }

  return std::make_pair(lowerList.size(), upperList.size());
}

template <class dataType>
std::pair<ttk::SimplexId, ttk::SimplexId> ttk::ScalarFieldCriticalPoints<
  dataType>::getNumberOfLowerUpperComponentsMultires(const SimplexId vertexId,
                                                     const MultiresTriangulation
                                                       *multiresTriangulation) {

  // cout<<"processing vertex "<<vertexId<<endl;
  SimplexId neighborNumber
    = multiresTriangulation->getVertexNeighborNumber(vertexId);
  // if(debugLevel_>4)
  if(vertexLinkPolarity_[vertexId].empty()) {
    buildVertexLinkPolarity(vertexId);
  }
#ifdef TTK_ENABLE_IMPLICIT_LINK_FOR_MULTIRES_CC
  associateVertexLinkByBoundary(vertexId);
#else
  buildVertexLink(vertexId);
#endif

  link_[vertexId].alloc(neighborNumber);

  // cout<<" init "<<vertexId<<" Link : \n"<<link_[vertexId].print()<<endl;

  for(SimplexId edgeId = 0; edgeId < (SimplexId)vertexLink_[vertexId]->size();
      edgeId++) {
    SimplexId n0 = vertexLink_[vertexId]->at(edgeId).first;
    SimplexId n1 = vertexLink_[vertexId]->at(edgeId).second;

    polarity isUpper0 = vertexLinkPolarity_[vertexId][n0].first;
    polarity isUpper1 = vertexLinkPolarity_[vertexId][n1].first;
    if(debugLevel_ > 6) {
      cout << "looking at edge " << n0 << "(" << (int)isUpper0 << ")-" << n1
           << "(" << (int)isUpper1 << ") " << endl;
    }
    if(vertexLinkPolarity_[vertexId][n0].first
       == vertexLinkPolarity_[vertexId][n1].first) {
      // cout<<"  - inserting edge "<<n0<<" "<<n1<<endl;
      // the smallest id (n0) becomes the parent of n1
      link_[vertexId].insertEdge(n1, n0);
      // cout<<" current "<<vertexId<<" Link :
      // \n"<<link_[vertexId].print()<<endl;
    }
  }
  // cout<<"LINK \n"<<link_[vertexId].print()<<endl;

  SimplexId downValence = 0, upValence = 0;
  std::tie(downValence, upValence) = getValencesFromLink(vertexId);
  if(debugLevel_ >= Debug::advancedInfoMsg) {
    std::stringstream msg;
    msg << "[ScalarFieldCriticalPoints] Vertex #" << vertexId
        << ": lowerLink-#CC=" << downValence << " upperLink-#CC=" << upValence
        << std::endl;

    dMsg(std::cout, msg.str(), Debug::advancedInfoMsg);
  }

  return std::make_pair(downValence, upValence);
}

// template <class dataType>
// bool ttk::ScalarFieldCriticalPoints<dataType>::checkChangeInNeighborsHeight(
//   vector<bool> localNeighborsHeight,
//   SimplexId vertexId,
//   SimplexId neighborNumber) {
//   // cout<<"checking for changes"<<endl;
//   for(int i = 0; i < neighborNumber; i++) {
//     if(localNeighborsHeight[i] != heightOfNeighbors_[vertexId][i]) {
//       heightOfNeighbors_[vertexId][i] = localNeighborsHeight[i];
//       return true;
//     }
//   }
//   // cout<<"\tdone checking"<<endl;
//   return false;
// }

// template <class dataType>
// void ttk::ScalarFieldCriticalPoints<dataType>::fillInfoOnNeighborsHeight(
//   vector<bool> localNeighborsHeight,
//   SimplexId vertexId,
//   SimplexId neighborNumber) {
//   heightOfNeighbors_[vertexId].resize(neighborNumber);
//   for(int i = 0; i < neighborNumber; i++) {
//       heightOfNeighbors_[vertexId][i] = localNeighborsHeight[i];
//   }
// }
template <class dataType>
bool ttk::ScalarFieldCriticalPoints<dataType>::getMonotonyChangesSecondVersion(
  SimplexId vertexId) {
  Timer t;

  bool hasMonotonyChanged = false;
  vector<SimplexId> globalImpactedIndices(0);
  vector<int> neighborsThatAreNew(0);
  vector<int> neighborsThatAreOld(0);
  // multiresTriangulation_.printInfos();
  SimplexId neighborNumber
    = multiresTriangulation_.getVertexNeighborNumber(vertexId);

  for(SimplexId neighborId = 0; neighborId < neighborNumber; neighborId++) {
    SimplexId globalNeighborId = -1;

    multiresTriangulation_.getVertexNeighbor(
      vertexId, neighborId, globalNeighborId);

    globalImpactedIndices.push_back(globalNeighborId);
    if(isNew_[globalNeighborId]) {
      neighborsThatAreNew.push_back(neighborId);
    } else {
      neighborsThatAreOld.push_back(neighborId);
    }
  }

  // cout<<"loop on old points"<<endl;
  for(SimplexId i : neighborsThatAreOld) {
    SimplexId invertedNeighborId = -1;
    // vertexBackMap_[vertexId][i];
    // cout<<"getting inverted indices"<<endl;
    multiresTriangulation_.getInvertVertexNeighbor(
      globalImpactedIndices[i], vertexId, invertedNeighborId);
    // cout<<"done"<<endl;
    bool lower
      = isSosLowerThan((*sosOffsets_)[vertexId], scalarValues_[vertexId],
                       (*sosOffsets_)[globalImpactedIndices[i]],
                       scalarValues_[globalImpactedIndices[i]]);
    polarity isUpper = lower ? 0 : 255;
    polarity isUpperOld
      = vertexLinkPolarity_[globalImpactedIndices[i]][invertedNeighborId].first;
    if(isUpper != isUpperOld) { // change of monotony
      hasMonotonyChanged = true;
      // cout<<"  wanna reprocess "<<globalImpactedIndices[i]<<" ?
      // "<<(int)(toReprocess_[globalImpactedIndices[i]])<<endl;
      toReprocess_[globalImpactedIndices[i]] = 255;
      // cout<<"  wanna reprocess "<<globalImpactedIndices[i]<<" ?
      // "<<(int)(toReprocess_[globalImpactedIndices[i]])<<endl;
      vertexLinkPolarity_[globalImpactedIndices[i]][invertedNeighborId].second
        = 255;
    }
  }
  // cout<<"loop on old points done"<<endl;
  if(hasMonotonyChanged) { // neighbors that are new must now be processed
    toProcess_[vertexId] = 255;
    // cout<<"loop on new points"<<endl;
    for(SimplexId i : neighborsThatAreNew) {
      toProcess_[globalImpactedIndices[i]] = 255;
    }
  }
  if(neighborsThatAreOld.size() > 2) {
    toProcess_[vertexId] = 255;
  }

  // time_get_monotony_changes+= t.getElapsedTime();
  return hasMonotonyChanged;
}

template <class dataType>
bool ttk::ScalarFieldCriticalPoints<dataType>::getMonotonyChanges(
  SimplexId vertexId, vector<char> & /*vertexTypes*/) {

  if(debugLevel_ > 9) {
    cout << " get monotony change" << endl;
  }
  bool hasMonotonyChanged = false;
  SimplexId v0[3];
  SimplexId v1[3];
  // multiresTriangulation_.printInfos();
  if(debugLevel_ > 9) {
    cout << "getting impacted vertices" << endl;
  }
  multiresTriangulation_.getImpactedVertices(vertexId, v0, v1);
  if(debugLevel_ > 9) {
    cout << "getting impacted vertices done" << endl;
  }

  SimplexId localNeighborId0 = v0[0];
  SimplexId localNeighborId1 = v1[0];
  SimplexId neighborId0 = v0[1];
  SimplexId neighborId1 = v1[1];
  SimplexId invertedLocalNeighborId0 = v0[2];
  SimplexId invertedLocalNeighborId1 = v1[2];
  if(debugLevel_ > 10) {
    cout << "Filled values : " << localNeighborId0 << " " << neighborId0 << " "
         << invertedLocalNeighborId0 << " " << localNeighborId1 << " "
         << neighborId1 << " " << invertedLocalNeighborId1 << " " << endl;
  }
  if(localNeighborId0 != -1 and localNeighborId1 != -1) {
    // check if the monotony of the edge has changed
    // cout<<"checking for a change "<<endl;

    // cout<<invertedLocalNeighborId0<<"
    // "<<vertexLinkPolarity_[neighborId0].size()<<"
    // "<<invertedLocalNeighborId1<<"
    // "<<vertexLinkPolarity_[neighborId1].size()<<endl;
    polarity isUpperOld0
      = vertexLinkPolarity_[neighborId0][invertedLocalNeighborId0].first;
    polarity isUpperOld1
      = vertexLinkPolarity_[neighborId1][invertedLocalNeighborId1].first;
    // cout<<"no issues on inverted neighbors"<<endl;

    bool lower0
      = isSosLowerThan((*sosOffsets_)[vertexId], scalarValues_[vertexId],
                       (*sosOffsets_)[neighborId0], scalarValues_[neighborId0]);
    polarity isUpper0 = lower0 ? 0 : 255;
    bool lower1
      = isSosLowerThan((*sosOffsets_)[vertexId], scalarValues_[vertexId],
                       (*sosOffsets_)[neighborId1], scalarValues_[neighborId1]);
    polarity isUpper1 = lower1 ? 0 : 255;
    // cout<<"   isupperOld0 "<<(int)isUpperOld0<<"   isupperOld1
    // "<<(int)isUpperOld1<<"   isupperd0 "<<(int)isUpper0<<"   isupper1
    // "<<(int)isUpper1<<endl;
    if(isUpperOld0 != isUpper0) { // the monotony has changed for v[0]
      hasMonotonyChanged = true;
      // if(enableTracking_ and vertexTypes[neighborId0]
      // !=static_cast<char>(CriticalType::Regular)){
      //   suspectsForTracking[neighborId0].push_back(vertexId);
      // }
      // cout<<"  -- "<<vertexId<<" breaks "<<neighborId0<<" on
      // "<<invertedLocalNeighborId0<<endl;
      toReprocess_[neighborId0] = 255;
      vertexLinkPolarity_[neighborId0][invertedLocalNeighborId0].second
        = 255; // BREAK of edge
    }

    if(isUpperOld1 != isUpper1) { // the monotony has changed for v[1]
      hasMonotonyChanged = true;
      // cout<<"  -- "<<vertexId<<" breaks "<<neighborId1<<" on
      // "<<invertedLocalNeighborId1<<endl;
      toReprocess_[neighborId1] = 255;
      vertexLinkPolarity_[neighborId1][invertedLocalNeighborId1].second
        = 255; // BREAK of edge
    }
  } else {
    cout
      << "[ScalarFieldCC, function getMonotonyChanges] How did you get there ? "
      << neighborId0 << " " << neighborId1 << endl;
    return 0;
  }
  if(debugLevel_ > 9) {
    cout << " all done got monotony change" << endl;
  }
  return hasMonotonyChanged;
}

// bool ScalarFieldCriticalPoints::getLinksBorder(SimplexId vertexId){
//   SimplexId neighborNumber =
//   multiresTriangulation_.getVertexNeighborNumber(vertexId); for(int
//   neighborId=0; neighborId<neighborNumber; neighborId++  ){

//   }
//   vector<SimplexId> commonNeighbors;
// }

// void constructNeighborWhiteList(SimplexId vertexId, MultiresTriangulation*
// multiresTriangulation){
//   SimplexId neighborNumber
//     = multiresTriangulation->getVertexNeighborNumber(vertexId);
//   for(int n0 = 0; n0 < neighborNumber; n0++) {
//     for(int n1 = 0; n1 < neighborNumber; n1++) {

//     }
//   }
template <class dataType>
string
  ttk::ScalarFieldCriticalPoints<dataType>::criticalTypeToString(char type) {
  switch(type) {
    case static_cast<char>(CriticalType::Local_minimum):
      return "minimum";
      break;

    case static_cast<char>(CriticalType::Saddle1):
      return "saddle1";
      break;

    case static_cast<char>(CriticalType::Saddle2):
      return "saddle2";
      break;

    case static_cast<char>(CriticalType::Local_maximum):
      return "maximum";
      break;

    case static_cast<char>(CriticalType::Degenerate):
      return "degenerate";
      break;
    case static_cast<char>(CriticalType::Regular):
      return "regular";
      break;
    default:
      return "";
      break;
  }
}

// template <class dataType>
// int ttk::ScalarFieldCriticalPoints<dataType>::checkSingleEdgeChange(SimplexId
// vertexId){
//   int nb_of_changes;
//   int i0;
//   for(int i = 0; i < vertexLinkPolarity_[vertexId].size(); i++){
//     i0 = i;
//     char changed = vertexLinkPolarity_[vertexId][i].second;
//     if(changed==255){
//       nb_of_changes++;
//       vertexLinkPolarity_[vertexId][i].second = 0;
//     }
//   }
//   if(nb_of_changes == 1)
//     return i0;
//   else
//     return -1;;

// }

// template <class dataType>
// int ttk::ScalarFieldCriticalPoints<dataType>::isEdgeOnLinkBoundary(
//   int edgeId, SimplexId globalId) {
//   for(int i = 0; i < vertexLink_[globalId].size(); i++) {
//     std::pair<char,char> linkEdge = vertexLink_[globalId][i];
//     if(linkEdge.first == edgeId){
//       if(vertexLinkPolarity_[globalId][linkEdge.first].first
//          != vertexLinkPolarity_[globalId][linkEdge.second].first) // border
//          of link
//         return true;
//     }
//   }
//   return false;
// }
//
template <class dataType>
char ttk::ScalarFieldCriticalPoints<dataType>::valencesToType(
  SimplexId downValence, SimplexId upValence) {

  if(downValence == -1 && upValence == -1) {
    return -1;
  } else if(downValence == 0 && upValence == 1) {
    return static_cast<char>(CriticalType::Local_minimum);
  } else if(downValence == 1 && upValence == 0) {
    return static_cast<char>(CriticalType::Local_maximum);
  } else if(downValence == 1 && upValence == 1) {
    // regular point
    return static_cast<char>(CriticalType::Regular);
  } else {
    // saddles
    if(dimension_ == 2) {
      if((downValence == 2 && upValence == 1)
         || (downValence == 1 && upValence == 2)
         || (downValence == 2 && upValence == 2)) {
        // regular saddle
        return static_cast<char>(CriticalType::Saddle1);
      } else {
        // monkey saddle, saddle + extremum
        return static_cast<char>(CriticalType::Degenerate);
        // NOTE: you may have multi-saddles on the boundary in that
        // configuration
        // to make this computation 100% correct, one would need to
        // disambiguate boundary from interior vertices
      }
    } else if(dimension_ == 3) {
      if(downValence == 2 && upValence == 1) {
        return static_cast<char>(CriticalType::Saddle1);
      } else if(downValence == 1 && upValence == 2) {
        return static_cast<char>(CriticalType::Saddle2);
      } else {
        // monkey saddle, saddle + extremum
        return static_cast<char>(CriticalType::Degenerate);
        // NOTE: we may have a similar effect in 3D (TODO)
      }
    }
  }

  // -2: regular points
  return static_cast<char>(CriticalType::Regular);
}

template <class dataType>
void ttk::ScalarFieldCriticalPoints<dataType>::buildVertexLinkPolarity(
  SimplexId vertexId) {
  // Timer t;
  SimplexId neighborNumber
    = multiresTriangulation_.getVertexNeighborNumber(vertexId);
  vertexLinkPolarity_[vertexId].resize(neighborNumber);

  for(SimplexId i = 0; i < neighborNumber; i++) {
    SimplexId neighborId0 = 0;
    multiresTriangulation_.getVertexNeighbor(vertexId, i, neighborId0);

    bool lower0
      = isSosLowerThan((*sosOffsets_)[neighborId0], scalarValues_[neighborId0],
                       (*sosOffsets_)[vertexId], scalarValues_[vertexId]);
    polarity isUpper0 = lower0 ? 0 : 255;
    vertexLinkPolarity_[vertexId][i] = std::make_pair(isUpper0, 0);
  }
  // time_building_link_polarity += t.getElapsedTime();
}

template <class dataType>
void ttk::ScalarFieldCriticalPoints<dataType>::buildVertexLink(
  SimplexId vertexId) {
  // Timer t;
  SimplexId neighborNumber
    = multiresTriangulation_.getVertexNeighborNumber(vertexId);
  vertexLinkExplicit_[vertexId].resize(0);

  for(SimplexId i = 0; i < neighborNumber; i++) {
    SimplexId neighborId0 = 0;
    multiresTriangulation_.getVertexNeighbor(vertexId, i, neighborId0);
    for(SimplexId j = i + 1; j < neighborNumber; j++) {
      SimplexId neighborId1 = 0;
      multiresTriangulation_.getVertexNeighbor(vertexId, j, neighborId1);
      if(multiresTriangulation_.areVerticesNeighbors(
           neighborId0, neighborId1)) {
        if(i < j)
          vertexLinkExplicit_[vertexId].push_back({i, j});
        else
          vertexLinkExplicit_[vertexId].push_back({j, i});
      }
    }
  }
  // time_building_link += t.getElapsedTime();
  vertexLink_[vertexId] = &(vertexLinkExplicit_[vertexId]);
}

template <class dataType>
void ttk::ScalarFieldCriticalPoints<dataType>::buildVertexLinkByBoundary(
  SimplexId vertexId) {
  // Timer t;
  int boundaryIndex = multiresTriangulation_.getVertexBoundaryIndex(vertexId);
  // cout << "constructing link for boundary index " << boundaryIndex
  //      << " (called by vertex " << vertexId << ")" << endl;
  SimplexId neighborNumber
    = multiresTriangulation_.getVertexNeighborNumber(vertexId);

  for(SimplexId i = 0; i < neighborNumber; i++) {
    SimplexId neighborId0 = 0;
    multiresTriangulation_.getVertexNeighbor(vertexId, i, neighborId0);
    for(SimplexId j = i + 1; j < neighborNumber; j++) {
      SimplexId neighborId1 = 0;
      multiresTriangulation_.getVertexNeighbor(vertexId, j, neighborId1);
      if(multiresTriangulation_.areVerticesNeighbors(
           neighborId0, neighborId1)) {
        if(i < j)
          vertexLinkByBoundaryType_[boundaryIndex].push_back({i, j});
        else
          vertexLinkByBoundaryType_[boundaryIndex].push_back({j, i});
      }
    }
  }
  // vertexLink_[vertexId] = &(vertexLinkByBoundaryType_[boundaryIndex]);
  // time_building_link += t.getElapsedTime();
}

template <class dataType>
void ttk::ScalarFieldCriticalPoints<dataType>::associateVertexLinkByBoundary(
  SimplexId vertexId) {
  // Timer t;
  int boundaryIndex = multiresTriangulation_.getVertexBoundaryIndex(vertexId);
  vertexLink_[vertexId] = &(vertexLinkByBoundaryType_[boundaryIndex]);
  // time_building_link += t.getElapsedTime();
}

template <class dataType>
std::pair<ttk::SimplexId, ttk::SimplexId>
  ttk::ScalarFieldCriticalPoints<dataType>::getValencesFromLink(
    SimplexId vertexId) {
  std::vector<size_t> CCIds;
  link_[vertexId].retrieveNbCC(CCIds);
  if(debugLevel_ > 4) {
    cout << "getting valences for " << CCIds.size() << "composants" << endl;
  }
  SimplexId downValence = 0, upValence = 0;

  for(size_t i = 0; i < CCIds.size(); i++) {
    SimplexId neighbor = CCIds[i];
    polarity isUpper = vertexLinkPolarity_[vertexId][neighbor].first;
    if(debugLevel_ > 4) {
      cout << "neighbor " << neighbor << " is a root. is it upper? "
           << (int)isUpper << endl;
    }
    if(isUpper != 0) {
      upValence++;
    } else {
      downValence++;
    }
  }
  // cout << "valences " << downValence << " " << upValence << endl;
  return std::make_pair(downValence, upValence);
}

template <class dataType>
char ttk::ScalarFieldCriticalPoints<dataType>::reProcess(SimplexId vertexId) {
  // loop on the link polarity

  // for each broken edge
  //  - modify the polarity
  // Timer t;
  if(debugLevel_ > 4) {
    cout << "reprocessing vertex " << vertexId << endl;
  }
  vector<SimplexId> monotony_changes_list(0);
  // cout<<" old "<<vertexId<<" Link : \n"<<link_[vertexId].print()<<endl;
  for(size_t neighborId = 0; neighborId < vertexLinkPolarity_[vertexId].size();
      neighborId++) {
    polarity isBroken = vertexLinkPolarity_[vertexId][neighborId].second;
    // polarity oldPolarity = vertexLinkPolarity_[vertexId][neighborId].first;
    if(isBroken != 0) {
      // cout<<"  - vertex "<<vertexId<<" is broken on Lneighbor
      // "<<neighborId<<endl;
      monotony_changes_list.push_back(neighborId);
      // link_[vertexId].getNode(neighborId)->evert();
    }
  }

  // loop on the link
  //   for each edge that shares n0
  //      if only one break and different polarity : remove
  //      else if only one break and same polarity : insert
  //      else : do nothing
  vector<vector<pair<SimplexId, SimplexId>>> edgesToInsertLater(
    vertexLinkPolarity_[vertexId].size());
  vector<vector<pair<SimplexId, SimplexId>>> edgesToRemoveLater(
    vertexLinkPolarity_[vertexId].size());

  for(size_t e = 0; e < vertexLink_[vertexId]->size(); e++) {
    SimplexId n0 = vertexLink_[vertexId]->at(e).first;
    SimplexId n1 = vertexLink_[vertexId]->at(e).second;
    polarity isBroken0 = vertexLinkPolarity_[vertexId][n0].second;
    polarity isBroken1 = vertexLinkPolarity_[vertexId][n1].second;

    if(isBroken0 != 0 and isBroken1 == 0) {
      if(vertexLinkPolarity_[vertexId][n0].first
         != vertexLinkPolarity_[vertexId][n1].first) {
        edgesToInsertLater[n0].push_back({n0, n1});
      } else {
        edgesToRemoveLater[n0].push_back({n0, n1});
      }
    }

    if(isBroken0 == 0 and isBroken1 != 0) {
      if(vertexLinkPolarity_[vertexId][n0].first
         != vertexLinkPolarity_[vertexId][n1].first) {
        edgesToInsertLater[n1].push_back({n1, n0});
      } else {
        edgesToRemoveLater[n1].push_back({n1, n0});
      }
    }
  }

  // cout<<" current "<<vertexId<<" Link : \n"<<link_[vertexId].print()<<endl;

  // REMOVE EDGES:
  for(SimplexId brokenNode : monotony_changes_list) {

    vertexLinkPolarity_[vertexId][brokenNode].first
      = ~vertexLinkPolarity_[vertexId][brokenNode].first;
    vertexLinkPolarity_[vertexId][brokenNode].second = 0;

    link_[vertexId].getNode(brokenNode)->evert();
    for(pair<SimplexId, SimplexId> edge : edgesToRemoveLater[brokenNode]) {
      link_[vertexId].removeEdge(edge.first, edge.second);
    }
  }
  if(edgesToRemoveLater.size() > 0)
    reConnectLink(vertexId);

  // INSERT EDGES
  for(SimplexId brokenNode : monotony_changes_list) {
    for(pair<SimplexId, SimplexId> edge : edgesToInsertLater[brokenNode]) {
      // cout<<"  - inserting edge "<<edge.first<<" "<<edge.second<<endl;
      link_[vertexId].insertEdge(edge.first, edge.second);
      // cout<<" current "<<vertexId<<" Link :
      // \n"<<link_[vertexId].print()<<endl;
    }
  }
  // int ret = -1;
  // if(vertexLinkPolarity_[vertexId][n0].second == 255){ //n0 broke
  //   SimplexId oldCC = link_[vertexId].getCCFromNode(n1);
  //   ret = link_[vertexId].removeEdge(n0, n1);
  //   cout<<"ret = "<<ret<<endl;
  //   if(ret == 2){ //the non-broken node was removed
  //     findReplacingEdge(vertexId, n1, n0);
  //   }
  // }else{ //n1 broke
  //   SimplexId oldCC = link_[vertexId].getCCFromNode(n0);
  //  ret = link_[vertexId].removeEdge(n1, n0);
  //   cout<<"ret = "<<ret<<endl;
  //   if(ret == 2){ //the non-broken node was removed
  //     findReplacingEdge(vertexId, n0, n1);
  //   }
  // }
  // }
  // }
  // }
  // for(int i = 0; i < (int)edgesToInsertLater.size(); i++) {
  // SimplexId n0 = edgesToInsertLater[i].first;
  // SimplexId n1 = edgesToInsertLater[i].second;
  // link_[vertexId].insertEdge(n1, n0, 1);
  // }
  // for(int neighborId = 0; neighborId<vertexLinkPolarity_[vertexId].size();
  // neighborId ++){ vertexLinkPolarity_[vertexId][neighborId].second = 0;
  // }

  SimplexId downValence = 0, upValence = 0;
  std::tie(downValence, upValence) = getValencesFromLink(vertexId);
  if(debugLevel_ > 4) {
    // cout<<" done"<<endl;
  }

  // cout<<" new "<<vertexId<<" Link : \n"<<link_[vertexId].print()<<endl;
  // time_reprocessing += t.getElapsedTime();
  char c = valencesToType(downValence, upValence);
  return c;
}

template <class dataType>
void ttk::ScalarFieldCriticalPoints<dataType>::findReplacingEdge(
  SimplexId vertexId, SimplexId removedNode, SimplexId oldCC) {
  for(size_t iEdge = 0; iEdge < vertexLink_[vertexId]->size(); iEdge++) {
    SimplexId n0 = vertexLink_[vertexId]->at(iEdge).first;
    SimplexId n1 = vertexLink_[vertexId]->at(iEdge).second;
    polarity neighborPolarity0 = vertexLinkPolarity_[vertexId][n0].first;
    polarity removedNodePolarity
      = vertexLinkPolarity_[vertexId][removedNode].first;
    polarity neighborPolarity1 = vertexLinkPolarity_[vertexId][n1].first;

    // cout << " searching replacement edge for " << removedNode << " : " << n0
    // << "-" << n1 << endl; cout << "polarity "<<(int)neighborPolarity0<<" -
    // "<<(int)neighborPolarity1<<"  ,removed node
    // "<<(int)removedNodePolarity<<endl; cout << "roots
    // "<<link_[vertexId].getCCFromNode(n0) <<" -
    // "<<link_[vertexId].getCCFromNode(n1) <<endl;
    if(link_[vertexId].getCCFromNode(n0) == static_cast<size_t>(removedNode)
       and neighborPolarity1 == removedNodePolarity) {
      // cout << " edge linked to removed node and good polarity " <<
      // (int)neighborPolarity1 << endl; cout << "composants " <<
      // link_[vertexId].getCCFromNode(n0) << " and " <<
      // link_[vertexId].getCCFromNode(n1) << ". " << oldCC << " wanted" <<
      // endl;
      if(static_cast<size_t>(oldCC) == link_[vertexId].getCCFromNode(n1)) {
        // cout << " replacement : inserting " << removedNode << "-" << n1 <<
        // endl;
        link_[vertexId].insertEdge(n0, n1);
        return;
      }
    } else if(link_[vertexId].getCCFromNode(n1)
                == static_cast<size_t>(removedNode)
              and neighborPolarity0 == removedNodePolarity) {
      // cout << " edge linked to removed node and good polarity " <<
      // (int)neighborPolarity0 << endl; cout << "composants " <<
      // link_[vertexId].getCCFromNode(n0) << " and " <<
      // link_[vertexId].getCCFromNode(n1) << ". " << oldCC << " wanted" <<
      // endl;
      if(static_cast<size_t>(oldCC) == link_[vertexId].getCCFromNode(n0)) {
        // cout << " replacement : inserting " << n0 << "-" << removedNode <<
        // endl;
        link_[vertexId].insertEdge(n1, n0);
        return;
      }
    }
  }
}

template <class dataType>
void ttk::ScalarFieldCriticalPoints<dataType>::reConnectLink(
  SimplexId vertexId) {
  for(pair<SimplexId, SimplexId> edge : *(vertexLink_[vertexId])) {
    if(vertexLinkPolarity_[vertexId][edge.first].first
       == vertexLinkPolarity_[vertexId][edge.second].first) {
      link_[vertexId].insertEdge(edge.first, edge.second);
    }
  }
}

// template <class dataType>
// void ttk::ScalarFieldCriticalPoints<dataType>::trackFormerPosition(SimplexId
// vertexId, char old_type, vector<char>& vertexTypes){
//   vector<SimplexId> extendedNeighbors =
//   multiresTriangulation_.getExtendedStar(vertexId); SimplexId
//   newCriticalPoint = -1; int count = 0; for(SimplexId neighbor :
//   extendedNeighbors){
//     cout<<" vertexId "<<vertexId<<" of type "<<(int)old_type<<" ,  neighbor
//     "<<neighbor<<" of type "<<(int)(vertexTypes[neighbor])<<endl;
//     if(vertexTypes[neighbor] == old_type){
//       count++;
//       newCriticalPoint = neighbor;
//     }
//   }
//    if(count < 1) {
//     cout<< "critical point of type  "<<old_type<<" has vanished"<<endl;
//   } else{
//     if(count>1){
//       cout << "Not able to determine the displacement for cc of type = "
//            << old_type << " (count " << count
//            << ") \n Several candidates, took the last one" << endl;
//     }
//     cout << "got the trajectory ! " << vertexId << " -> " << newCriticalPoint
//          << endl;
//     tempTrajectories_[decimationLevel_ - stoppingDecimationLevel_ ]
//       .push_back(
//         {decimationLevel_, vertexId, newCriticalPoint, old_type, vertexId});
//     criticalGeneration_[newCriticalPoint] = criticalGeneration_[vertexId];
//     // cout<<"pushed back tuple"<<endl;
//   }
//   criticalGeneration_[vertexId] = -1;
// }

template <class dataType>
dataType ttk::ScalarFieldCriticalPoints<dataType>::integrateToMaximum(
  SimplexId vertexId,
  vector<char> &vertexTypes,
  vector<SimplexId> *integralPath) {
  if(vertexTypes[vertexId]
     == static_cast<char>(ttk::CriticalType::Local_maximum)) { // found a max
    // tempTrajectories_[decimationLevel_ - stoppingDecimationLevel_].push_back(
    //   {decimationLevel_, integralPath,
    //    static_cast<char>(ttk::CriticalType::Local_maximum),
    //    integralPath[0]});
    integralPath->push_back(vertexId);
    return scalarValues_[vertexId];
  }
  if(vertexTypes[vertexId] == static_cast<char>(ttk::CriticalType::Saddle1)
     or vertexTypes[vertexId]
          == static_cast<char>(ttk::CriticalType::Saddle2)) {
    if(debugLevel_ > 4) {
      cout
        << " WATCH OUT, on our way to a max ,we integrated to a saddle after "
        << integralPath->size() << " STEPS" << endl;
    }
    integralPath->push_back(vertexId);
    // cout<<" sized switched to "<<integralPath->size()<<endl;
    return integrateToMaximumFromSaddle(vertexId, vertexTypes, integralPath);

  } else {
    // find highest neighbor
    // using the link
    SimplexId neighborNumber
      = multiresTriangulation_.getVertexNeighborNumber(vertexId);
    SimplexId neighborId0 = -1;
    multiresTriangulation_.getVertexNeighbor(vertexId, 0, neighborId0);
    SimplexId maxNeighbor = neighborId0;
    for(SimplexId i = 1; i < neighborNumber; i++) {
      SimplexId neighborId = -1;
      multiresTriangulation_.getVertexNeighbor(vertexId, i, neighborId);
      // dataType value = scalarValues_[neighborId];
      if(isSosHigherThan((*sosOffsets_)[neighborId], scalarValues_[neighborId],
                         (*sosOffsets_)[vertexId], scalarValues_[vertexId])) {
        maxNeighbor = neighborId;
      }
    }
    integralPath->push_back(vertexId);
    return integrateToMaximum(maxNeighbor, vertexTypes, integralPath);
  }
}

template <class dataType>
dataType ttk::ScalarFieldCriticalPoints<dataType>::integrateToMaximumFromSaddle(
  SimplexId vertexId,
  vector<char> &vertexTypes,
  vector<SimplexId> *integralPath) {
  std::vector<size_t> CCIds;
  link_[vertexId].retrieveNbCC(CCIds);
  dataType maxValue = std::numeric_limits<dataType>::min();
  vector<SimplexId> maxPath;
  for(size_t i = 0; i < CCIds.size(); i++) {
    SimplexId neighbor = CCIds[i];
    polarity isUpper = vertexLinkPolarity_[vertexId][neighbor].first;
    if(isUpper != 0) {
      SimplexId neighborId;
      multiresTriangulation_.getVertexNeighbor(vertexId, neighbor, neighborId);
      // cout<<"saddle "<<vertexId<<" - launching integration from
      // "<<neighborId<<endl;
      vector<SimplexId> path(0);
      cout << " appel from saddle max " << ++COUNT << " size "
           << integralPath->size() << endl;
      // if(COUNT>0) return std::numeric_limits<double>::min();
      dataType value = integrateToMaximum(neighborId, vertexTypes, &path);
      if(value > maxValue) {
        maxValue = value;
        maxPath.resize(path.size());
        maxPath.swap(path);
      }
    }
  }
  for(SimplexId id : maxPath) {
    integralPath->push_back(id);
  }
  return maxValue;
}

template <class dataType>
dataType ttk::ScalarFieldCriticalPoints<dataType>::integrateToMinimum(
  SimplexId vertexId,
  vector<char> &vertexTypes,
  vector<SimplexId> *integralPath) {

  if(vertexTypes[vertexId]
     == static_cast<char>(ttk::CriticalType::Local_minimum)) { // found a max
    // tempTrajectories_[decimationLevel_ - stoppingDecimationLevel_].push_back(
    //   {decimationLevel_, integralPath,
    //    static_cast<char>(ttk::CriticalType::Local_maximum),
    //    integralPath[0]});
    integralPath->push_back(vertexId);
    return scalarValues_[vertexId];
  }
  if(vertexTypes[vertexId] == static_cast<char>(ttk::CriticalType::Saddle1)
     or vertexTypes[vertexId]
          == static_cast<char>(ttk::CriticalType::Saddle2)) {
    if(debugLevel_ > 4) {
      cout
        << " WATCH OUT, on our way to a min ,we integrated to a saddle after "
        << integralPath->size() << " STEPS" << endl;
    }
    integralPath->push_back(vertexId);
    //
    return integrateToMinimumFromSaddle(vertexId, vertexTypes, integralPath);

  } else {
    // find highest neighbor
    // using the link
    SimplexId neighborNumber
      = multiresTriangulation_.getVertexNeighborNumber(vertexId);
    SimplexId neighborId0 = -1;
    multiresTriangulation_.getVertexNeighbor(vertexId, 0, neighborId0);
    SimplexId minNeighbor = neighborId0;
    for(SimplexId i = 1; i < neighborNumber; i++) {
      SimplexId neighborId = -1;
      multiresTriangulation_.getVertexNeighbor(vertexId, i, neighborId);
      // dataType value = scalarValues_[neighborId];
      // if(value < minValue) {
      if(isSosLowerThan((*sosOffsets_)[neighborId], scalarValues_[neighborId],
                        (*sosOffsets_)[vertexId], scalarValues_[vertexId])) {
        minNeighbor = neighborId;
      }
    }
    integralPath->push_back(vertexId);
    return integrateToMinimum(minNeighbor, vertexTypes, integralPath);
  }
}
template <class dataType>
dataType ttk::ScalarFieldCriticalPoints<dataType>::integrateToMinimumFromSaddle(
  SimplexId vertexId,
  vector<char> &vertexTypes,
  vector<SimplexId> *integralPath) {

  std::vector<size_t> CCIds;
  link_[vertexId].retrieveNbCC(CCIds);
  dataType minValue = std::numeric_limits<dataType>::max();
  vector<SimplexId> minPath;
  for(size_t i = 0; i < CCIds.size(); i++) {
    SimplexId neighbor = CCIds[i];
    polarity isUpper = vertexLinkPolarity_[vertexId][neighbor].first;
    if(isUpper == 0) {
      SimplexId neighborId;
      multiresTriangulation_.getVertexNeighbor(vertexId, neighbor, neighborId);
      // cout<<"saddle "<<vertexId<<" - launching integration from
      // "<<neighborId<<endl;
      vector<SimplexId> path(0);
      cout << " appel from saddle min " << ++COUNT << " size "
           << integralPath->size() << endl;
      // if(COUNT>0) return std::numeric_limits<double>::max();
      dataType value = integrateToMinimum(neighborId, vertexTypes, &path);
      if(value < minValue) {
        minValue = value;
        minPath.resize(path.size());
        minPath.swap(path);
      }
    }
  }
  for(SimplexId id : minPath) {
    integralPath->push_back(id);
  }
  return minValue;
}
// template <class dataType>
// dataType
// ttk::ScalarFieldCriticalPoints<dataType>::integrateToMinimumIterative(
//  SimplexId vertexId,
//  vector<char> &vertexTypes,
//  vector<SimplexId> *integralPath) {
//  cout<<"integrating from vertex "<<vertexId<< " step
//  "<<integralPath->size()<<endl;

//  SimplexId currentId = vertexId;
//  integralPath->push_back(currentId);

//  while(vertexTypes[currentId] !=
//  static_cast<char>(ttk::CriticalType::Local_minimum)){ // found a max
//    //find lowest edge
//    if(vertexTypes[currentId] == static_cast<char>(ttk::CriticalType::Saddle1)
//    or vertexTypes[currentId] ==
//    static_cast<char>(ttk::CriticalType::Saddle2)) {
//      integrateToMinimumFromSaddleIterative(currentId, vertexTypes,
//      integralPath);
//    }
//    SimplexId neighborNumber =
//    multiresTriangulation_.getVertexNeighborNumber(currentId); SimplexId
//    neighborId0 = -1; multiresTriangulation_.getVertexNeighbor(currentId, 0,
//    neighborId0); dataType minValue = scalarValues_[neighborId0]; SimplexId
//    minNeighborId = 0; SimplexId minNeighbor = neighborId0;

//    for(unsigned int i = 1; i<neighborNumber; i++){
//      // cout<<" request "<<i<<" outta "<<neighborNumber<<endl;
//      SimplexId neighborId = -1;
//      multiresTriangulation_.getVertexNeighbor(currentId, i, neighborId);
//      // cout<<"fetching value"<<endl;
//      dataType value = scalarValues_[neighborId];
//      // cout<<"fetched value"<<endl;
//      if(value < minValue) {
//        minNeighborId = i;
//        minValue = value;
//        minNeighbor = neighborId;
//      }
//    }
//    currentId = minNeighbor;
//    integralPath->push_back(currentId);
//  }
//  return scalarValues_[currentId];
//  if(vertexTypes[vertexId] == static_cast<char>(ttk::CriticalType::Saddle1) or
//  vertexTypes[vertexId] == static_cast<char>(ttk::CriticalType::Saddle2)){
//    if(debugLevel_>4){
//    cout<<" WATCH OUT, on our way to a min ,we integrated to a saddle after
//    "<<integralPath->size()<<" STEPS"<<endl;
//    }
//    integralPath->push_back(vertexId);
//    cout<<" sized switched to "<<integralPath->size()<<endl;
//    //
//    return integrateToMinimumFromSaddle(vertexId, vertexTypes, integralPath);

//  } else {
//    //find highest neighbor
//    // using the link
//    cout<<"pushing back"<<endl;
//    integralPath->push_back(vertexId);
//    cout<<"pushed back"<<endl;
//    return integrateToMinimum(minNeighbor, vertexTypes, integralPath);
//  }
//}

// template <class dataType>
// void ttk::ScalarFieldCriticalPoints<dataType>::findCriticalSon(SimplexId
// vertexId, vector<char>& vertexTypes){
//   char old_type = vertexTypes[vertexId];
//   vector<SimplexId> extendedNeighbors =
//   multiresTriangulation_.getExtendedStar(vertexId); SimplexId
//   newCriticalPoint = -1; int count = 0; for(SimplexId neighbor :
//   extendedNeighbors){
//     cout<<"neighbor "<<neighbor<<endl;
//     if(vertexTypes[neighbor] == old_type){
//       count++;
//       newCriticalPoint = neighbor;
//     }
//   }
//   if(count>1){
//     cout << "Not able to determine the son for cc of type = "
//       << old_type << " (count " << count << ")"<< endl;
//   } else if(count < 1) {
//     cout<< "critical point of type  "<<old_type<<" has no son !"<<endl;
//   } else{
//     cout << "got the son ! " << vertexId << " -> " << newCriticalPoint
//          << endl;
//     criticalGeneration_[newCriticalPoint] = criticalGeneration_[vertexId] +
//     1; criticalGeneration_[vertexId] = 0;
//   }
// }

// template <class dataType>
// void ttk::ScalarFieldCriticalPoints<dataType>::postProcessTrajectories(){
//   int start = startingDecimationLevel_ - stoppingDecimationLevel_ - 1;
//   int end  = 0;
//   for(int d = start; d > end; d--) {
//     for(int i = 0; i < tempTrajectories_[d].size(); i++) {
//       trajTuple ti = tempTrajectories_[d][i];
//       for(int j = 0; j < tempTrajectories_[d - 1].size(); j++) {
//         trajTuple tj = tempTrajectories_[d - 1][j];
//         cout<<" size "<<get<1>(tj).size()<<endl;
//         if(get<1>(ti).back() == get<1>(tj).at(0)) {
//           get<3>(tj) = get<3>(ti);
//           tempTrajectories_[d - 1][j] = tj;
//         }
//       }
//     }
//   }
//   for(int d = start; d >= end; d--) {
//     for(int i = 0; i < tempTrajectories_[d].size(); i++) {
//       trajectories_->push_back(tempTrajectories_[d][i]);
//     }
//   }
// }
template <class dataType>
void ttk::ScalarFieldCriticalPoints<dataType>::printTrajectories() {
  cout << "Printing trajectories : " << endl;
  for(size_t i = 0; i < trajectoriesMax_->size(); i++) {
    cout << " - traj " << i << " : ";
    for(size_t j = 0; j < trajectoriesMax_->at(i).size(); j++)
      cout << " - " << trajectoriesMax_->at(i)[j];
  }
  cout << "\n" << endl;
}

template <class dataType>
void ttk::ScalarFieldCriticalPoints<dataType>::updateLinkPolarity(
  SimplexId vertexId) {
  for(unsigned int i = 0; i < vertexLinkPolarity_[vertexId].size(); i++) {
    SimplexId neighborId = -1;
    multiresTriangulation_.getVertexNeighbor(vertexId, i, neighborId);
    bool lower
      = isSosLowerThan((*sosOffsets_)[neighborId], scalarValues_[neighborId],
                       (*sosOffsets_)[vertexId], scalarValues_[vertexId]);
    polarity isUpper = lower ? 0 : 255;
    vertexLinkPolarity_[vertexId][i] = make_pair(isUpper, 0);
  }
}

template <class dataType>
void ttk::ScalarFieldCriticalPoints<dataType>::printPointInfos(
  SimplexId vertexId) {
  if(vertexId < vertexNumber_) {
    SimplexId neighborNumber
      = multiresTriangulation_.getVertexNeighborNumber(vertexId);
    cout << " === INFOS === vertex " << vertexId << endl;
    SimplexId p[3];
    multiresTriangulation_.vertexToPosition(vertexId, p);
    cout << " == position " << p[0] << " " << p[1] << " " << p[2] << endl;
    cout << " == scalar value " << scalarValues_[vertexId] << endl;
    cout << " == " << neighborNumber << " neighbors" << endl;
    for(int i = 0; i < neighborNumber; i++) {
      SimplexId neighborId = -1;
      multiresTriangulation_.getVertexNeighbor(vertexId, i, neighborId);
      cout << "   " << i << " - " << neighborId << " of value "
           << scalarValues_[neighborId] << endl;
    }
    // cout << " == link infos :" << endl;
    // for(int i = 0; i < vertexLink_[vertexId]->size(); i++) {
    //   SimplexId n0;
    //   SimplexId n1;
    //   multiresTriangulation_.getVertexNeighbor(
    //     vertexId, vertexLink_[vertexId]->at(i).first, n0);
    //   multiresTriangulation_.getVertexNeighbor(
    //     vertexId, vertexLink_[vertexId]->at(i).second, n1);
    //   cout << "   edge " << n0 << " - " << n1 << endl;
    // }
  }
}

template <class dataType>
void ttk::ScalarFieldCriticalPoints<dataType>::buildVertexBackMap(
  const SimplexId vertexId) {
  SimplexId neighborNumber
    = multiresTriangulation_.getVertexNeighborNumber(vertexId);
  vertexBackMap_[vertexId].resize(neighborNumber);
  for(SimplexId i = 0; i < neighborNumber; i++) {
    SimplexId globalNeighborId = -1;
    SimplexId inverseLocalNeighbor = -1;
    multiresTriangulation_.getVertexNeighbor(vertexId, i, globalNeighborId);
    multiresTriangulation_.getInvertVertexNeighbor(
      globalNeighborId, vertexId, inverseLocalNeighbor);
    vertexBackMap_[vertexId][i] = inverseLocalNeighbor;
  }
}

template <class dataType>
bool ttk::ScalarFieldCriticalPoints<dataType>::getMonotonyChangeByOldPoint(
  SimplexId vertexId) {
  bool hasMonotonyChanged = false;
  SimplexId neighborNumber
    = multiresTriangulation_.getVertexNeighborNumber(vertexId);
  for(SimplexId i = 0; i < neighborNumber; i++) {
    SimplexId neighborId = -1;
    multiresTriangulation_.getVertexNeighbor(vertexId, i, neighborId);
    // check for monotony changes

    bool lower
      = isSosLowerThan((*sosOffsets_)[neighborId], scalarValues_[neighborId],
                       (*sosOffsets_)[vertexId], scalarValues_[vertexId]);
    polarity isUpper = lower ? 0 : 255;
    polarity isUpperOld = vertexLinkPolarity_[vertexId][i].first;

    if(isUpper != isUpperOld) { // change of monotony
      hasMonotonyChanged = true;

      toReprocess_[vertexId] = 255;
      toProcess_[neighborId] = 255;
      SimplexId neighborNumberNew
        = multiresTriangulation_.getVertexNeighborNumber(neighborId);
      for(SimplexId j = 0; j < neighborNumberNew; j++) {
        SimplexId neighborIdNew = -1;
        multiresTriangulation_.getVertexNeighbor(neighborId, j, neighborIdNew);
        if(isNew_[neighborIdNew])
          toProcess_[neighborIdNew] = 255;
      }

      vertexLinkPolarity_[vertexId][i].second = 255;
    }
  }
  return hasMonotonyChanged;
}

template class ttk::ScalarFieldCriticalPoints<double>;
template class ttk::ScalarFieldCriticalPoints<float>;
template class ttk::ScalarFieldCriticalPoints<int>;
template class ttk::ScalarFieldCriticalPoints<char>;
template class ttk::ScalarFieldCriticalPoints<unsigned int>;
template class ttk::ScalarFieldCriticalPoints<unsigned char>;
template class ttk::ScalarFieldCriticalPoints<signed char>;
template class ttk::ScalarFieldCriticalPoints<short>;
template class ttk::ScalarFieldCriticalPoints<unsigned short>;
template class ttk::ScalarFieldCriticalPoints<long>;
template class ttk::ScalarFieldCriticalPoints<unsigned long>;
template class ttk::ScalarFieldCriticalPoints<long long>;
template class ttk::ScalarFieldCriticalPoints<unsigned long long>;
// #endif /* end of include guard: SCALARFIELDCRITICALPOINTS_INL */
