#include <ttkScalarFieldCriticalPoints.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkScalarFieldCriticalPoints)

  ttkScalarFieldCriticalPoints::ttkScalarFieldCriticalPoints() {

  // init
  ForceInputOffsetScalarField = false;
  VertexBoundary = true;
  VertexIds = true;
  VertexScalars = true;

  ScalarFieldId = 0;
  OffsetFieldId = -1;
  OffsetField = ttk::OffsetScalarFieldName;

  UseMultiresTriangulation = false;
  UseProgressive = false;
  IntegrateTrajectories = false;
  DecimationLevel = 0;
  StoppingDecimationLevel = 0;
  UseAllCores = true;
  SetNumberOfOutputPorts(2);
}

ttkScalarFieldCriticalPoints::~ttkScalarFieldCriticalPoints() {
}

template <typename VTK_TT>
int ttkScalarFieldCriticalPoints::dispatch(Triangulation *triangulation,
                                           void *scalarValues,
                                           const SimplexId vertexNumber) {
  ScalarFieldCriticalPoints<VTK_TT> criticalPoints;
  criticalPoints.setupTriangulation(triangulation);
  int domainDimension = triangulation->getCellVertexNumber(0) - 1;

  criticalPoints.setWrapper(this);
  criticalPoints.setDomainDimension(domainDimension);
  // set up input
  criticalPoints.setUseMultiresTriangulation(UseMultiresTriangulation);
  criticalPoints.setDecimationLevel(DecimationLevel);
  criticalPoints.setStoppingDecimationLevel(StoppingDecimationLevel);
  // 1 -- vertex values
  criticalPoints.setScalarValues(scalarValues);
  criticalPoints.setVertexNumber(vertexNumber);

  // 2 -- set offsets (here, let the baseCode class fill it for us)
  criticalPoints.setSosOffsets(&sosOffsets_);

  // 3 -- set the connectivity
  criticalPoints.setupTriangulation(triangulation);

  // set up output
  criticalPoints.setOutput(&criticalPoints_);
  if(UseProgressive) {
    // criticalPoints.setUseProgressive(UseProgressive);
    criticalPoints.setIntegrateTrajectories(IntegrateTrajectories);
    criticalPoints.setOutputTrajectoriesMax(&trajectoriesMax_, &trajectoriesMaxDecimation_);
    criticalPoints.setOutputTrajectoriesMin(&trajectoriesMin_, &trajectoriesMinDecimation_);
    criticalPoints.setOutputCriticalGeneration(&criticalGenerationOutput_);
    criticalPoints.executeProgressive();

    
    // vtkSmartPointer<vtkImageData> imData = vtkSmartPointer<vtkImageData>::New();
    // std::vector<int> dimensions(3);
    // dimensions
    //   = criticalPoints.getMultiresTriangulation()->getGridDimensions();
    // imData->SetDimensions(dimensions[0], dimensions[1], dimensions[2]);
    // cout<<dimensions[0]<<" "<<dimensions[1]<<" "<<dimensions[2]<<endl;
    // vtkSmartPointer<vtkIntArray> data_array = vtkSmartPointer<vtkIntArray>::New();
    // int vertexNumbers = criticalPoints.getMultiresTriangulation()->getVertexNumber();
    // data_array->SetNumberOfComponents(1);
    // data_array->SetNumberOfTuples(vertexNumbers);
    // for(int v = 0; v < vertexNumbers; v++) {
    //   data_array->SetTuple1(v, criticalPoints.getProcessingsPerVertex()->at(v));
    // }
    // imData->GetPointData()->SetScalars(data_array);

    // std::vector<float> spacing = criticalPoints.getMultiresTriangulation()
    //                                ->getTriangulation()
    //                                ->getImplicitSpacing();
    // imData->SetSpacing(spacing[0],spacing[1],spacing[2]);

    // std::vector<float> origin = criticalPoints.getMultiresTriangulation()
    //                                ->getTriangulation()
    //                                ->getOrigin();
    // imData->SetOrigin(origin[0],origin[1],origin[2]);

    // vtkSmartPointer<vtkXMLImageDataWriter> writer
    //   = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    // writer->SetFileName("./numberOfReprocess.vti");
    // writer->SetInputData(imData);
    // writer->Write();

  } else {
    criticalPoints.execute();
  }

  return 0;
}

int ttkScalarFieldCriticalPoints::doIt(vector<vtkDataSet *> &inputs,
                                       vector<vtkDataSet *> &outputs) {
  Memory m;
  Timer t;

#ifndef TTK_ENABLE_KAMIKAZE
  if(!inputs.size()) {
    cerr
      << "[ttkScalarFieldCriticalPoints] Error: not enough input information."
      << endl;
    return -1;
  }
#endif

  vtkDataSet *input = inputs[0];
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(outputs[0]);
  vtkUnstructuredGrid *output_trajectories = vtkUnstructuredGrid::SafeDownCast(outputs[1]);

#ifndef TTK_ENABLE_KAMIKAZE
  if(!input) {
    cerr << "[ttkScalarFieldCriticalPoints] Error: input pointer is NULL."
         << endl;
    return -1;
  }

  if(!input->GetNumberOfPoints()) {
    cerr << "[ttkScalarFieldCriticalPoints] Error: input has no point." << endl;
    return -1;
  }
#endif

  Triangulation *triangulation = ttkTriangulation::getTriangulation(input);
#ifndef TTK_ENABLE_KAMIKAZE
  if(!triangulation) {
    cerr << "[ttkScalarFieldCriticalPoints] Error: input triangulation is NULL."
         << endl;
    return -1;
  }
#endif

  if(VertexBoundary)
    triangulation->preprocessBoundaryVertices();

  // in the following, the target scalar field of the input is replaced in the
  // variable 'output' with the result of the computation.
  // if your wrapper produces an output of the same type of the input, you
  // should proceed in the same way.
  vtkDataArray *inputScalarField = NULL;
  vtkDataArray *offsetField = NULL;

  if(ScalarField.length()) {
    inputScalarField = input->GetPointData()->GetArray(ScalarField.data());
  } else {
    inputScalarField = input->GetPointData()->GetArray(ScalarFieldId);
  }

  if(!inputScalarField)
    return -1;

  {
    stringstream msg;
    msg << "[ttkScalarFieldCriticalPoints] Starting computation on field `"
        << inputScalarField->GetName() << "'..." << endl;
    dMsg(cout, msg.str(), infoMsg);
  }

  if(OffsetFieldId != -1) {
    offsetField = input->GetPointData()->GetArray(OffsetFieldId);
    if(offsetField) {
      ForceInputOffsetScalarField = true;
      OffsetField = offsetField->GetName();
    }
  }

  if(ForceInputOffsetScalarField) {
    if(OffsetField.length()) {

      offsetField = input->GetPointData()->GetArray(OffsetField.data());
      // not good... in the future, we want to use the pointer itself...
      sosOffsets_.resize(offsetField->GetNumberOfTuples());
      for(SimplexId i = 0; i < offsetField->GetNumberOfTuples(); i++) {
        SimplexId offset = 0;
        offset = offsetField->GetTuple1(i);
        sosOffsets_[i] = offset;
      }
    }
  } else if(input->GetPointData()->GetArray(ttk::OffsetScalarFieldName)) {
    offsetField = input->GetPointData()->GetArray(OffsetScalarFieldName);

    // not good... in the future, we want to use the pointer itself...
    sosOffsets_.resize(offsetField->GetNumberOfTuples());
    for(SimplexId i = 0; i < offsetField->GetNumberOfTuples(); i++) {
      SimplexId offset = 0;
      offset = offsetField->GetTuple1(i);
      sosOffsets_[i] = offset;
    }
  }

  switch(inputScalarField->GetDataType()) {
    vtkTemplateMacro(dispatch<VTK_TT>(triangulation,
                                      inputScalarField->GetVoidPointer(0),
                                      input->GetNumberOfPoints()));
  }
  // dispatch<double>(triangulation, inputScalarField->GetVoidPointer(0), input->GetNumberOfPoints());

  // allocate the output
  vtkSmartPointer<vtkCharArray> vertexTypes
    = vtkSmartPointer<vtkCharArray>::New();

  vertexTypes->SetNumberOfComponents(1);
  vertexTypes->SetNumberOfTuples(criticalPoints_.size());
  vertexTypes->SetName("CriticalType");

  vtkSmartPointer<vtkPoints> pointSet = vtkSmartPointer<vtkPoints>::New();
  pointSet->SetNumberOfPoints(criticalPoints_.size());
  double p[3];
  for(SimplexId i = 0; i < (SimplexId)criticalPoints_.size(); i++) {
    input->GetPoint(criticalPoints_[i].first, p);
    pointSet->SetPoint(i, p);
    vertexTypes->SetTuple1(i, (float)criticalPoints_[i].second);
  }
  output->SetPoints(pointSet);
  output->GetPointData()->AddArray(vertexTypes);

  if(UseProgressive){
    vtkSmartPointer<vtkIntArray> criticalGenerationArray = vtkSmartPointer<vtkIntArray>::New();
    criticalGenerationArray->SetNumberOfTuples(criticalPoints_.size());
    criticalGenerationArray->SetNumberOfComponents(1);
    criticalGenerationArray->SetName("CriticalGeneration");

    for(SimplexId i = 0; i < (SimplexId)criticalPoints_.size(); i++) {
      criticalGenerationArray->SetTuple1(i, criticalGenerationOutput_[i]);
    }
    output->GetPointData()->AddArray(criticalGenerationArray);
  }

  if(VertexBoundary) {
    vtkSmartPointer<vtkCharArray> vertexBoundary
      = vtkSmartPointer<vtkCharArray>::New();
    vertexBoundary->SetNumberOfComponents(1);
    vertexBoundary->SetNumberOfTuples(criticalPoints_.size());
    vertexBoundary->SetName("IsOnBoundary");

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
    for(SimplexId i = 0; i < (SimplexId)criticalPoints_.size(); i++) {
      vertexBoundary->SetTuple1(
        i, (char)triangulation->isVertexOnBoundary(criticalPoints_[i].first));
    }

    output->GetPointData()->AddArray(vertexBoundary);
  } else {
    output->GetPointData()->RemoveArray("IsOnBoundary");
  }

  if(VertexIds) {
    vtkSmartPointer<ttkSimplexIdTypeArray> vertexIds
      = vtkSmartPointer<ttkSimplexIdTypeArray>::New();
    vertexIds->SetNumberOfComponents(1);
    vertexIds->SetNumberOfTuples(criticalPoints_.size());
    vertexIds->SetName(ttk::VertexScalarFieldName);

    for(SimplexId i = 0; i < (SimplexId)criticalPoints_.size(); i++) {
      vertexIds->SetTuple1(i, criticalPoints_[i].first);
    }

    output->GetPointData()->AddArray(vertexIds);
  } else {
    output->GetPointData()->RemoveArray(ttk::VertexScalarFieldName);
  }

  if(VertexScalars) {
    for(SimplexId i = 0; i < input->GetPointData()->GetNumberOfArrays(); i++) {

      vtkDataArray *scalarField = input->GetPointData()->GetArray(i);
      vtkSmartPointer<vtkDataArray> scalarArray;

      auto copyToScalarArray = [&]() {
        scalarArray->SetNumberOfComponents(
          scalarField->GetNumberOfComponents());
        scalarArray->SetNumberOfTuples(criticalPoints_.size());
        scalarArray->SetName(scalarField->GetName());
        std::vector<double> value(scalarField->GetNumberOfComponents());
        for(SimplexId j = 0; j < (SimplexId)criticalPoints_.size(); j++) {
          scalarField->GetTuple(criticalPoints_[j].first, value.data());
          scalarArray->SetTuple(j, value.data());
        }
        output->GetPointData()->AddArray(scalarArray);
      };

      switch(scalarField->GetDataType()) {
        case VTK_CHAR:
          scalarArray = vtkSmartPointer<vtkCharArray>::New();
          copyToScalarArray();
          break;
        case VTK_DOUBLE:
          scalarArray = vtkSmartPointer<vtkDoubleArray>::New();
          copyToScalarArray();
          break;
        case VTK_FLOAT:
          scalarArray = vtkSmartPointer<vtkFloatArray>::New();
          copyToScalarArray();
          break;
        case VTK_INT:
          scalarArray = vtkSmartPointer<vtkIntArray>::New();
          copyToScalarArray();
          break;
        case VTK_ID_TYPE:
          scalarArray = vtkSmartPointer<vtkIdTypeArray>::New();
          copyToScalarArray();
          break;
        default: {
          stringstream msg;
          msg << "[ttkScalarFieldCriticalPoints] Scalar attachment: "
              << "unsupported data type :(" << endl;
          dMsg(cerr, msg.str(), detailedInfoMsg);
        } break;
      }
    }
  } else {
    for(SimplexId i = 0; i < input->GetPointData()->GetNumberOfArrays(); i++) {
      output->GetPointData()->RemoveArray(
        input->GetPointData()->GetArray(i)->GetName());
    }
  }

  //fill trajectories output
  // cout<<"filling trajectories output"<<endl;
  vtkSmartPointer<vtkPoints> trajectoriesPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkIntArray> lineDecimation = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkCharArray> lineType = vtkSmartPointer<vtkCharArray>::New();
  vtkSmartPointer<vtkIntArray> linePointDecimation = vtkSmartPointer<vtkIntArray>::New();

  lineDecimation->SetNumberOfComponents(1);
  // lineDecimation->SetNumberOfTuples(trajectories_.size());
  lineDecimation->SetName("DecimationLevel");

  lineType->SetNumberOfComponents(1);
  // lineType->SetNumberOfTuples(trajectories_.size());
  lineType->SetName("CriticalType");

  linePointDecimation->SetNumberOfComponents(1);
  // linePointDecimation->SetNumberOfTuples(2 * trajectories_.size());
  linePointDecimation->SetName("decimationLevel");
  output_trajectories->Reset();
  output_trajectories->Allocate();
  // list<unsigned long int> pointList;
  int count = 0;
  for(size_t j = 0; j<trajectoriesMax_.size(); j++){
    for(unsigned int k = 0; k<trajectoriesMax_[j].size()-1; k++){ // fillin max trajectories
      double p1[3];
      double p2[3];
      input->GetPoint(trajectoriesMax_[j][k], p1);
      input->GetPoint(trajectoriesMax_[j][k + 1], p2);
      int id1 = trajectoriesPoints->InsertNextPoint(p1);
      int id2 = trajectoriesPoints->InsertNextPoint(p2);
      vtkIdType lineIds[2];
      lineIds[0] = id1;
      lineIds[1] = id2;
      output_trajectories->InsertNextCell(VTK_LINE, 2, lineIds);
      lineDecimation->InsertNextTuple1(trajectoriesMaxDecimation_[j]);
      lineType->InsertNextTuple1(static_cast<char>(ttk::CriticalType::Local_maximum));
      count++;
    }
  }
  for(size_t j = 0; j<trajectoriesMin_.size(); j++){
    for(unsigned int k = 0; k < trajectoriesMin_[j].size() - 1;
        k++) { // fillin min trajectories
      double p1[3];
      double p2[3];
      input->GetPoint(trajectoriesMin_[j][k], p1);
      input->GetPoint(trajectoriesMin_[j][k + 1], p2);
      int id1 = trajectoriesPoints->InsertNextPoint(p1);
      int id2 = trajectoriesPoints->InsertNextPoint(p2);
      vtkIdType lineIds[2];
      lineIds[0] = id1;
      lineIds[1] = id2;
      output_trajectories->InsertNextCell(VTK_LINE, 2, lineIds);
      lineDecimation->InsertNextTuple1(trajectoriesMinDecimation_[j]);
      lineType->InsertNextTuple1(static_cast<char>(ttk::CriticalType::Local_minimum));
      count++;
    }
  }
  output_trajectories->SetPoints(trajectoriesPoints);
  // output_trajectories->GetPointData()->AddArray(linePointDecimation);
  output_trajectories->GetCellData()->AddArray(lineType);
  output_trajectories->GetCellData()->AddArray(lineDecimation);
  // cout<<"done"<<endl;
  // pointList.sort();
  // pointList.unique();
  // for(int i_point=0; i_point<pointList.size(); i_point++){
  // }

  {
    stringstream msg;
    msg << "[ttkScalarFieldCriticalPoints] Memory usage: "
        << m.getElapsedUsage() << " MB." << endl;
    dMsg(cout, msg.str(), 2);
  }

  return 0;
}
