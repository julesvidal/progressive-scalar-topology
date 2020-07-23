/// \ingroup base
/// \class ttk::ScalarFieldCriticalPoints
///
/// \brief TTK processing package for the computation of critical points in PL
/// scalar fields defined on PL manifolds.
///
/// This class computes the list of critical points of the input scalar field
/// and classify them according to their type.
///
/// \param dataType Data type of the input scalar field (char, float,
/// etc.).
///
///
/// \sa ttkScalarFieldCriticalPoints.cpp %for a usage example.

#ifndef _SCALARFIELDCRITICALPOINTS_H
#define _SCALARFIELDCRITICALPOINTS_H

#include <map>

// base code includes
#include <DynamicTree.h>
#include <MultiresTriangulation.h>
#include <Triangulation.h>
#include <UnionFind.h>
#include <Wrapper.h>

namespace ttk {

  using polarity = unsigned char;
  using trajTuple = std::tuple<int, std::vector<SimplexId>, char, SimplexId>;

  template <class dataType>
  class ScalarFieldCriticalPoints : public Debug {

  public:
    ScalarFieldCriticalPoints();

    ~ScalarFieldCriticalPoints();

    /// Execute the package.
    /// \param argment Dummy integer argument.
    /// \return Returns 0 upon success, negative values otherwise.
    int execute();
    int executeProgressive();

    std::pair<SimplexId, SimplexId>
      getNumberOfLowerUpperComponents(const SimplexId vertexId,
                                      Triangulation *triangulation);

    std::pair<SimplexId, SimplexId>
      getNumberOfLowerUpperComponentsMultiresWithUnionFinds(
        const SimplexId vertexId);
    std::pair<SimplexId, SimplexId>
      getNumberOfLowerUpperComponentsMultiresWithUnionFindsAndUsingLinks(
        const SimplexId vertexId);

    std::pair<SimplexId, SimplexId> getNumberOfLowerUpperComponentsMultires(
      const SimplexId vertexId,
      const MultiresTriangulation *multiresTriangulation);

    bool getMonotonyChanges(SimplexId vertexId_, std::vector<char> &);
    bool checkChangeInNeighborsHeight(std::vector<bool>, SimplexId, SimplexId);
    void fillInfoOnNeighborsHeight(std::vector<bool>, SimplexId, SimplexId);

    char getCriticalType(const SimplexId &vertexId) {

      return getCriticalType(vertexId, triangulation_);
    }

    char getCriticalType(const SimplexId &vertexId,
                         Triangulation *triangulation);

    char processVertex(const SimplexId &vertexId, Triangulation *triangulation);

    char getCriticalType(
      const SimplexId &vertexId,
      const std::vector<std::pair<SimplexId, SimplexId>> &vertexLinkEdgeList);

    static bool isSosHigherThan(const SimplexId &offset0,
                                const dataType &value0,
                                const SimplexId &offset1,
                                const dataType &value1) {

      return ((value0 > value1) || ((value0 == value1) && (offset0 > offset1)));
    }

    static bool isSosLowerThan(const SimplexId &offset0,
                               const dataType &value0,
                               const SimplexId &offset1,
                               const dataType &value1) {

      return ((value0 < value1) || ((value0 == value1) && (offset0 < offset1)));
    }

    int setDomainDimension(const int &dimension) {

      dimension_ = dimension;

      return 0;
    }

    int setOutput(std::vector<std::pair<SimplexId, char>> *criticalPoints) {

      criticalPoints_ = criticalPoints;

      return 0;
    }

    int setOutputTrajectoriesMax(
      std::vector<std::vector<SimplexId>> *trajectories,
      std::vector<int> *trajectoriesDecimation) {
      trajectoriesMax_ = trajectories;
      trajectoriesMaxDecimation_ = trajectoriesDecimation;
      return 0;
    }

    int setOutputTrajectoriesMin(
      std::vector<std::vector<SimplexId>> *trajectories,
      std::vector<int> *trajectoriesDecimation) {
      trajectoriesMin_ = trajectories;
      trajectoriesMinDecimation_ = trajectoriesDecimation;
      return 0;
    }

    int
      setOutputCriticalGeneration(std::vector<int> *criticalGenerationOutput) {
      criticalGenerationOutput_ = criticalGenerationOutput;
      return 0;
    }

    int setupTriangulation(Triangulation *triangulation) {

      triangulation_ = triangulation;

      // pre-condition functions
      if(triangulation_) {
        triangulation_->preprocessVertexNeighbors();
        triangulation_->preprocessVertexStars();
      }

      return 0;
    }

    int setScalarValues(const void *data) {

      scalarValues_ = (const dataType *)data;

      return 0;
    }

    int setSosOffsets(std::vector<SimplexId> *offsets) {

      sosOffsets_ = offsets;

      return 0;
    }

    int setVertexLinkEdgeLists(
      const std::vector<std::vector<std::pair<SimplexId, SimplexId>>>
        *edgeList) {

      vertexLinkEdgeLists_ = edgeList;

      return 0;
    }

    /// Set the number of vertices in the scalar field.
    /// \param vertexNumber Number of vertices in the data-set.
    /// \return Returns 0 upon success, negative values otherwise.
    int setVertexNumber(const SimplexId &vertexNumber) {
      vertexNumber_ = vertexNumber;
      return 0;
    }

    void setNonManifold(const bool b) {
      forceNonManifoldCheck = b;
    }

    void setUseMultiresTriangulation(const bool b) {
      useMultiresTriangulation_ = b;
    }

    void setIntegrateTrajectories(const bool b) {
      integrateTrajectories_ = b;
    }

    // void setUseProgressive(const bool b){
    //   useProgressive_ = b;
    // }

    void setStoppingDecimationLevel(const int stoppingDecimationLevel) {
      stoppingDecimationLevel_ = stoppingDecimationLevel;
    }

    void setDecimationLevel(const int decimationLevel) {
      decimationLevel_ = decimationLevel;
    }

    std::string criticalTypeToString(char);

    MultiresTriangulation *getMultiresTriangulation() {
      return &multiresTriangulation_;
    }

    std::vector<int> *getProcessingsPerVertex() {
      return &processingsPerVertex_;
    }

    // int checkSingleEdgeChange(SimplexId vertexId);

    // int ScalarFieldCriticalPoints<dataType>::isEdgeOnLinkBoundary(int,
    // SimplexId);
    char valencesToType(SimplexId downValence, SimplexId upValence);
    void buildVertexLink(SimplexId vertexId);
    void buildVertexLinkByBoundary(SimplexId vertexId);
    void associateVertexLinkByBoundary(SimplexId vertexId);
    void buildVertexLinkPolarity(SimplexId vertexId);
    void findReplacingEdge(SimplexId vertexId, SimplexId n, SimplexId cc);

    std::pair<ttk::SimplexId, ttk::SimplexId>
      getValencesFromLink(SimplexId vertexId);

    char reProcess(SimplexId vertexId);
    void reConnectLink(SimplexId vertexId);
    // void trackFormerPosition(SimplexId vertexId,
    //                          char old_type,
    //                          vector<char> &vertexTypes);
    void postProcessTrajectories();

    dataType integrateToMaximum(SimplexId vertexId,
                                std::vector<char> &vertexTypes,
                                std::vector<SimplexId> *integralPath);
    dataType integrateToMaximumFromSaddle(SimplexId vertexId,
                                          std::vector<char> &vertexTypes,
                                          std::vector<SimplexId> *integralPath);
    dataType integrateToMinimum(SimplexId vertexId,
                                std::vector<char> &vertexTypes,
                                std::vector<SimplexId> *integralPath);
    dataType integrateToMinimumFromSaddle(SimplexId vertexId,
                                          std::vector<char> &vertexTypes,
                                          std::vector<SimplexId> *integralPath);
    void printTrajectories();
    void printPointInfos(SimplexId);
    bool getMonotonyChangesSecondVersion(SimplexId vertexId);
    bool getMonotonyChangeByOldPoint(SimplexId vertexId);
    void updateLinkPolarity(SimplexId);
    void buildVertexBackMap(const SimplexId vertexId);
    char processVertexWithUnionFinds(const SimplexId &vertexId);

  protected:
    bool useMultiresTriangulation_;
    bool enableTracking_;
    int COUNT;
    std::vector<trajTuple> *trajectories_;
    std::vector<std::vector<trajTuple>> tempTrajectories_;
    std::vector<std::vector<SimplexId>> *trajectoriesMin_;
    std::vector<int> *trajectoriesMinDecimation_;
    std::vector<std::vector<SimplexId>> *trajectoriesMax_;
    std::vector<int> *trajectoriesMaxDecimation_;
    std::vector<int> criticalGeneration_;
    std::vector<int> *criticalGenerationOutput_;
    std::vector<SimplexId> criticalFather_;
    std::vector<std::vector<SimplexId>> vertexBackMap_;
    std::vector<bool> verticesToExamine_;
    std::vector<bool> isNew_;
    std::vector<unsigned char> toReprocess_;
    std::vector<unsigned char> toProcess_;
    std::vector<unsigned char> toCheck_;
    std::vector<int> nbOfMonotonyChanges_;
    std::vector<int> processingsPerVertex_;
    std::vector<DynamicTree> link_;

    MultiresTriangulation multiresTriangulation_;
    std::vector<std::vector<std::pair<SimplexId, SimplexId>> *> vertexLink_;
    std::vector<std::vector<std::pair<SimplexId, SimplexId>>>
      vertexLinkByBoundaryType_;
    std::vector<std::vector<std::pair<SimplexId, SimplexId>>>
      vertexLinkExplicit_;

    // nb_points * nb_neighbors * < previously upper link ? , currently upper
    // link >
    std::vector<std::vector<std::pair<polarity, polarity>>> vertexLinkPolarity_;

    int decimationLevel_;
    int stoppingDecimationLevel_;
    int startingDecimationLevel_;
    bool integrateTrajectories_;

    int dimension_;
    SimplexId vertexNumber_;
    const dataType *scalarValues_;
    const std::vector<std::vector<std::pair<SimplexId, SimplexId>>>
      *vertexLinkEdgeLists_;
    std::vector<std::pair<SimplexId, char>> *criticalPoints_;
    std::vector<SimplexId> *sosOffsets_;
    std::vector<SimplexId> localSosOffSets_;
    Triangulation *triangulation_;

    bool forceNonManifoldCheck;
    double time_processing;
    double time_reprocessing;
    double time_building_link;
    double time_building_link_polarity;
    double time_get_monotony_changes;
  };
} // namespace ttk

// if the package is not a template, comment the following line
// #include <ScalarFieldCriticalPoints.inl>

#endif // SCALARFIELDCRITICALPOINTS_H
