/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkProfileFilter.cxx,v $
=========================================================================*/
#include "vtkProfileFilter.h"
#include "AstroVizHelpersLib/AstroVizHelpers.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkCellData.h"
#include "vtkTable.h"
#include "vtkSortDataArray.h"
#include "vtkMath.h"
#include "vtkInformationDataObjectKey.h"
#include "vtkPointSet.h" 
#include "vtkMultiProcessController.h"
#include "vtkSmartPointer.h"
#include "vtkPointData.h"
#include "vtkLine.h"
#include "vtkPlane.h"
#include "vtkTimerLog.h"
#include <cmath>
#include <algorithm>
using vtkstd::string;

vtkCxxRevisionMacro(vtkProfileFilter, "$Revision: 1.72 $");
vtkStandardNewMacro(vtkProfileFilter);
vtkCxxSetObjectMacro(vtkProfileFilter, Controller, vtkMultiProcessController);
//----------------------------------------------------------------------------
vtkProfileFilter::vtkProfileFilter()
{
  this->SetNumberOfInputPorts(2);
  this->SetInputArrayToProcess(
    0,
    0,
    0,
    vtkDataObject::FIELD_ASSOCIATION_POINTS_THEN_CELLS,
    vtkDataSetAttributes::SCALARS);
  // Defaults for quantities which will be computed based on user's
  // later input
  this->MaxR=1.0;
  this->Delta=1; 
  this->BinNumber=30;
  this->Controller = NULL;
  this->SetController(vtkMultiProcessController::GetGlobalController());
}

//----------------------------------------------------------------------------
vtkProfileFilter::~vtkProfileFilter()
{
   this->SetController(0);
}

//----------------------------------------------------------------------------
void vtkProfileFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  os << indent << "bin number: " << this->BinNumber << "\n";
}

//----------------------------------------------------------------------------
void vtkProfileFilter::SetSourceConnection(vtkAlgorithmOutput* algOutput)
{
  this->SetInputConnection(1, algOutput);
}

//----------------------------------------------------------------------------
int vtkProfileFilter::FillInputPortInformation (int port, 
  vtkInformation *info)
{
  this->Superclass::FillInputPortInformation(port, info);
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  return 1;
}

//----------------------------------------------------------------------------
int vtkProfileFilter::RequestData(vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  vtkSmartPointer<vtkTimerLog> timer = vtkSmartPointer<vtkTimerLog>::New();
  timer->StartTimer();

  // Get the input with which we want to work
  vtkPointSet* input = vtkPointSet::GetData(inputVector[0]);
  // Will set the center based upon the selection in the GUI
  vtkDataSet* centerInfo = vtkDataSet::GetData(inputVector[1]);
  // Get name of data array containing mass
  vtkDataArray* massArray = this->GetInputArrayToProcess(0, inputVector);
  if (!massArray)
    {
    vtkErrorMacro("Failed to locate mass array");
    return 0;
    }
  vtkTable* const output = vtkTable::GetData(outputVector,0);
  output->Initialize();

  this->ProfileQuantities.clear();
  this->AdditionalProfileQuantities.clear();

  //
  // Initialize all arrays/profiles locally. 
  // Much faster than broadcasting.
  //

  //
  // create default profile arrays
  //
  this->ProfileQuantities.push_back(
    ProfileElement("bin radius", 1,   
    TOTAL));

  this->ProfileQuantities.push_back(
    ProfileElement("bin radius min", 1,   
    TOTAL));

  this->ProfileQuantities.push_back(
    ProfileElement("number in bin", 1,   
    TOTAL));

  this->ProfileQuantities.push_back(
    ProfileElement("number in bin", 1,   
    CUMULATIVE));

  //
  // Create 3 profile objects for each point array
  //
  vtkPointData *pd = input->GetPointData();
  for (int i=0; i<pd->GetNumberOfArrays(); ++i)
  {
    vtkSmartPointer<vtkDataArray> nextArray = pd->GetArray(i);
    int numComponents = nextArray->GetNumberOfComponents();
    const char *baseName = nextArray->GetName();
    // TOTAL, AVERAGE, CUMULATIVE
    this->ProfileQuantities.push_back(
      ProfileElement(baseName, numComponents,   
      TOTAL));
    this->ProfileQuantities.push_back(
      ProfileElement(baseName, numComponents,   
      AVERAGE));
    this->ProfileQuantities.push_back(
      ProfileElement(baseName, numComponents,   
      CUMULATIVE));
  }

  // Create Additional profiles
  // Note : make sure we reserve space before getting profile refs
  // because vector::push_back() can reallocate space and 
  // invalidate earlier ones when size changes
  this->AdditionalProfileQuantities.reserve(15);

  // DoubleFunction0
  this->AdditionalProfileQuantities.push_back(
    ProfileElement("angular momentum",3,   
    &ComputeAngularMomentum, 
    AVERAGE));
  this->AdditionalProfileQuantities.push_back(
    ProfileElement("radial velocity",3,
    &ComputeRadialVelocity, 
    AVERAGE));
  this->AdditionalProfileQuantities.push_back(
    ProfileElement("tangential velocity",3,
    &ComputeTangentialVelocity, 
    AVERAGE));

  // DoubleFunction2
  this->AdditionalProfileQuantities.push_back(
    ProfileElement("velocity squared",1,
    &ComputeVelocitySquared,
    AVERAGE));

  // DoubleFunction1
  this->AdditionalProfileQuantities.push_back(
    ProfileElement("radial velocity squared",1,
    &ComputeRadialVelocitySquared,
    AVERAGE));
  this->AdditionalProfileQuantities.push_back(
    ProfileElement("tangential velocity squared",1,
    &ComputeTangentialVelocitySquared,
    AVERAGE));

  //
  // These profiles refer to other profiles, so we use 
  // references to the ones we've already created for fast access.
  // We need to search the profile arrays because
  // the order is not defined. Fetch the profile now
  // and we can directly use the arrays it holds later.
  //
  ProfileElement::PIterator mass = std::find_if(
    this->ProfileQuantities.begin(), this->ProfileQuantities.end(),
    ProfileElement::FindProfile(massArray->GetName(),CUMULATIVE)
  );
  ProfileElement::PIterator binradius = std::find_if(
    this->ProfileQuantities.begin(), this->ProfileQuantities.end(),
    ProfileElement::FindProfile("bin radius",TOTAL)
  );
  ProfileElement::PIterator velocity = std::find_if(
    this->ProfileQuantities.begin(), this->ProfileQuantities.end(),
    ProfileElement::FindProfile("velocity",AVERAGE)
  );
  ProfileElement::PIterator velocity2 = std::find_if(
    this->AdditionalProfileQuantities.begin(), this->AdditionalProfileQuantities.end(),
    ProfileElement::FindProfile("velocity squared",AVERAGE)
  );
  ProfileElement::PIterator tangential = std::find_if(
    this->AdditionalProfileQuantities.begin(), this->AdditionalProfileQuantities.end(),
    ProfileElement::FindProfile("tangential velocity",AVERAGE)
  );
  ProfileElement::PIterator tangential2 = std::find_if(
    this->AdditionalProfileQuantities.begin(), this->AdditionalProfileQuantities.end(),
    ProfileElement::FindProfile("tangential velocity squared",AVERAGE)
  );
  ProfileElement::PIterator radial = std::find_if(
    this->AdditionalProfileQuantities.begin(), this->AdditionalProfileQuantities.end(),
    ProfileElement::FindProfile("radial velocity",AVERAGE)
  );
  ProfileElement::PIterator radial2 = std::find_if(
    this->AdditionalProfileQuantities.begin(), this->AdditionalProfileQuantities.end(),
    ProfileElement::FindProfile("radial velocity squared",AVERAGE)
  );

  // PostProcessFunc
  this->AdditionalProfileQuantities.push_back(
    ProfileElement("velocity dispersion",3,
    &ComputeVelocityDispersion,
    velocity2, velocity));
  this->AdditionalProfileQuantities.push_back(
    ProfileElement("tangential velocity dispersion",3,
    &ComputeVelocityDispersion,
    tangential2, tangential));
  this->AdditionalProfileQuantities.push_back(
    ProfileElement("radial velocity dispersion",3,
    &ComputeVelocityDispersion,
    radial2, radial));
  this->AdditionalProfileQuantities.push_back(
    ProfileElement("circular velocity",1,
    &ComputeCircularVelocity,
    mass, binradius));
  this->AdditionalProfileQuantities.push_back(
    ProfileElement("density",1,
    &ComputeDensity,
    mass, binradius));

  //
  // Allocate vtk Arrays and add to output vtkTable
  //
  for (int i=0; i<this->ProfileQuantities.size(); i++) {
    vtkFloatArray *farray = this->ProfileQuantities[i].AllocateBinArray(this->BinNumber);
    output->AddColumn(farray);
  }

  for (int i=0; i<this->AdditionalProfileQuantities.size(); ++i) {
    vtkFloatArray *farray = this->AdditionalProfileQuantities[i].AllocateBinArray(this->BinNumber);
    output->AddColumn(farray);
  }

  // setting the bin radii in the output
  ProfileElement &binradii = this->ProfileQuantities[BinRadius];
  ProfileElement &binradiimin = this->ProfileQuantities[BinRadiusMin];
  this->SetBoundsAndBinExtents(input,centerInfo); 
  for(int binNum = 0; binNum < this->BinNumber; ++binNum)
  {
    double radiusmax = (binNum+1)*this->BinSpacing;
    double radiusmin =  binNum*this->BinSpacing;
    this->UpdateBin(binNum, SET, binradii,    radiusmax);
    this->UpdateBin(binNum, SET, binradiimin, radiusmin);     
  }  

  //
  // Compute statistics for all data on this process
  //
  this->UpdateStatistics(input,output);

  // runs in parallel, syncing class member data, if necessary, if not
  // functions in serial
  if (RunInParallel(this->Controller))
    {
    int procId=this->Controller->GetLocalProcessId();
    int numProc=this->Controller->GetNumberOfProcesses();
    if (procId==0)
      {
      // Receive computations from each process and merge the tables into
      // the output of process 0
      for (int proc = 1; proc < numProc; ++proc)
        {
        vtkSmartPointer<vtkTable> recLocalTable = vtkSmartPointer<vtkTable>::New();
        recLocalTable->Initialize();
        this->Controller->Receive(recLocalTable,proc,DATA_TABLE);
        this->MergeTables(input,output,recLocalTable);
        }
      // Perform final computations
      // Updating averages and doing relevant postprocessing
      this->BinAveragesAndPostprocessing(input,output);
      }
    else
      {
      // send result to root
      this->Controller->Send(output,0,DATA_TABLE);
      // final merged result is on process 0
      // so all other processes should output empty vtkTable data 
      output->Initialize();
      }
    }  
  else
    {
    // Updating averages and doing relevant postprocessing
    this->BinAveragesAndPostprocessing(input,output);
    }

  //
  // we don't need to hold onto all this data, so release it.
  //
  this->ProfileQuantities.clear();
  this->AdditionalProfileQuantities.clear();

  timer->StopTimer();
  if (this->Controller->GetLocalProcessId()==0) { 
    std::cout << "ProfileFilter : " << timer->GetElapsedTime() << " seconds" << std::endl;
  }
  return 1;
}

//----------------------------------------------------------------------------
void vtkProfileFilter::MergeTables(vtkPointSet* input,
  vtkTable* originalTable, vtkTable* tableToMerge)
{
  assert(originalTable->GetNumberOfRows()==tableToMerge->GetNumberOfRows());
  //
  // The table we has received should have identical columns to our local table
  // so we can loop over all our (local) profiles and simply sum data from
  // the remote table using a simple lookup

  // We skip Bin Radius and Bin Radius min because they do not need to be merged.

  // Loop over arrays first, then elements so that we can cache the data pointers
  // and do a whole array at once.
  for (int i=NumberInBinTotal; i<this->ProfileQuantities.size(); i++) {
    ProfileElement &profile = this->ProfileQuantities[i];
    //
    // Get the data pointers for local and (remote) data to be merged
    //
    vtkFloatArray *updateArray = vtkFloatArray::SafeDownCast(
      tableToMerge->GetColumnByName(profile.GetColumnName()));

    if (!updateArray || updateArray->GetNumberOfTuples()!=this->BinNumber) {
      vtkErrorMacro(<<"Fatal : The columns being merged are not correct");
      return;
    }
    //
    // We've now got the two array pointers. Loop aver bins doing the merge
    //
    for (int binNum = 0; binNum<this->BinNumber; ++binNum) {
      float *updateData = updateArray->GetPointer(binNum*profile.NumberComponents);
      this->UpdateBin<float>(binNum, ADD, profile, updateData);
    }
  }

  //
  // Cut and paste from above using the additional profiles array (from zero)
  //
  for (int i=0; i<this->AdditionalProfileQuantities.size(); ++i) {
    ProfileElement &profile = this->AdditionalProfileQuantities[i];
    //
    // Get the data pointers for local and (remote) data to be merged
    //
    vtkFloatArray *updateArray = vtkFloatArray::SafeDownCast(
      tableToMerge->GetColumnByName(profile.GetColumnName())
    );
    if (!updateArray || updateArray->GetNumberOfTuples()!=this->BinNumber) {
      vtkErrorMacro(<<"Fatal : The columns being merged are not correct");
      return;
    }
    //
    // We've now got the two array pointers. Loop aver bins doing the merge
    //
    for (int binNum = 0; binNum<this->BinNumber; ++binNum) {
      float *updateData = updateArray->GetPointer(binNum*profile.NumberComponents);
      this->UpdateBin<float>(binNum, ADD, profile, updateData);
    }
  }
}

//----------------------------------------------------------------------------
void vtkProfileFilter::SetBoundsAndBinExtents(vtkPointSet* input, 
  vtkDataSet* source)
{
  if(RunInParallel(this->Controller))
    {
    int procId=this->Controller->GetLocalProcessId();
    int numProc=this->Controller->GetNumberOfProcesses();
    if(procId==0)
      {
      source->GetCenter(this->Center);
      // Syncronizing the centers
      this->Controller->Broadcast(this->Center,3,0);      
      }
    else
      {
      // Syncronizing the centers
      this->Controller->Broadcast(this->Center,3,0);
      }
    // calculating the max R
    this->MaxR=ComputeMaxRadiusInParallel(this->Controller,
      input,this->Center);
    }
  else
    {
    // we aren't using MPI or have only one process
    source->GetCenter(this->Center);
    //calculating the the max R
    this->MaxR=ComputeMaxR(input,this->Center);      
    }
  // this->MaxR, this->BinNumber are already set/synced
  // whether we are in parallel or serial, and each process can perform
  // computation on its own
  this->BinSpacing=this->CalculateBinSpacing(this->MaxR,this->BinNumber);
}

//----------------------------------------------------------------------------
int vtkProfileFilter::GetBinNumber(double x[])
{
  double distanceToCenter;
  if(this->ProfileAxis[0]==0 && this->ProfileAxis[1]==0 && this->ProfileAxis[2]==0) {// spherical profile
    distanceToCenter = sqrt(vtkMath::Distance2BetweenPoints(this->Center,x));
  }  
  else {
    distanceToCenter = vtkLine::DistanceToLine(x, this->Center, this->ProfileAxis);
  }
  
  if(this->ProfileHeight!=0) {
    if(this->ProfileAxis[0]==0 && this->ProfileAxis[1]==0 && this->ProfileAxis[2]==0){ //spherical profile
      if(distanceToCenter > this->ProfileHeight){
        return -1;
      }
    }
    else {
      double projectedPoint[3];
      vtkPlane::ProjectPoint(x, this->Center, this->ProfileAxis, projectedPoint); // TODO: this ProfileAxis is not the correct normal vector 
      if( sqrt(vtkMath::Distance2BetweenPoints(x,projectedPoint)) > this->ProfileHeight) {
        return -1;
      }
    }
    
  }  
    
  return floor(distanceToCenter/this->BinSpacing);
}

//----------------------------------------------------------------------------
double vtkProfileFilter::CalculateBinSpacing(double maxR,int binNumber)
{
  return maxR/binNumber;
}

//----------------------------------------------------------------------------
void vtkProfileFilter::UpdateStatistics(vtkPointSet* input, vtkTable* output)
{
  vtkPoints      *points = input->GetPoints();
  vtkIdType            N = points->GetNumberOfPoints();
  vtkDataArray *velocity = input->GetPointData()->GetArray("velocity");
  for(vtkIdType nextPointId = 0; nextPointId < N; ++nextPointId)
    {
      this->UpdateBinStatistics(input,points, velocity, nextPointId, output);
    }
}

//----------------------------------------------------------------------------
void vtkProfileFilter::UpdateBinStatistics(vtkPointSet* input, 
  vtkPoints *points, vtkDataArray *velocity,
  vtkIdType pointGlobalId, vtkTable* output)
{
  double r[3], v[3], x[3];
	points->GetPoint(pointGlobalId,x);
  //
  int binNum=this->GetBinNumber(x);
  if (binNum < 0){
    // This indicates the point is not to be included.
    return;
  }
  // As we bin by radius always need
  vtkMath::Subtract(x,this->Center, r);
  // Many of the quantities explicitely require the velocity
  velocity->GetTuple(pointGlobalId, v);

  //
  // Update bin counts
  //
  double value = 1.0;
  ProfileElement &numBinTotal = this->ProfileQuantities[NumberInBinTotal];
  ProfileElement &numBinCumul = this->ProfileQuantities[NumberInBinCumulative];
  this->UpdateBin(binNum, ADD, numBinTotal, value);  
  this->UpdateCumulativeBins(binNum,ADD, numBinCumul, value);  

  //
  // Updating quanties for the input data arrays
  //
  double value3[3];
  vtkPointData *pd = input->GetPointData();
  for(int i=0; i<pd->GetNumberOfArrays(); ++i)
  {
    // *3 because TOTAL, AVERAGE, CUMULATIVE
    ProfileElement &profileT = this->ProfileQuantities[PointDataArrays+i*3+0];
    ProfileElement &profileA = this->ProfileQuantities[PointDataArrays+i*3+1];
    ProfileElement &profileC = this->ProfileQuantities[PointDataArrays+i*3+2];

    vtkDataArray *nextArray = pd->GetArray(i);
    nextArray->GetTuple(pointGlobalId,value3);
    this->UpdateBin(binNum,ADD,profileT,value3);  
    this->UpdateBin(binNum,ADD,profileA,value3);  
    this->UpdateCumulativeBins(binNum,ADD,profileC,value3);
  }

  //
  // Updating additional profiles
  //
  double newValues[3];
  for(int i=0; i<this->AdditionalProfileQuantities.size(); ++i)
  {
    ProfileElement &profile=this->AdditionalProfileQuantities[i];
    switch (profile.FuntionType) {
      case FUNC0_TYPE:
        profile.func0(v,r,newValues);
        break;
      case FUNC1_TYPE:
        newValues[0] = profile.func1(v,r);
        break;
      case FUNC2_TYPE:
        newValues[0] = profile.func2(v);
        break;
    }
    if (profile.FuntionType!=POSTPROCESS_TYPE) {
      this->UpdateBin(binNum,ADD,profile,newValues);
    }
  }
}

//----------------------------------------------------------------------------
void vtkProfileFilter::BinAveragesAndPostprocessing(
  vtkPointSet* input,vtkTable* output)
{
  vtkPointData *pd = input->GetPointData();
  ProfileElement &numBinTotal = this->ProfileQuantities[NumberInBinTotal];
  //
  for(int binNum = 0; binNum < this->BinNumber; ++binNum) 
  { 
    double binTotal = numBinTotal.DataPointer[binNum];
    // For each input array, update its average column by getting
    // the total from the total column then dividing by 
    // the number in the bin. Only do this if the number in the bin
    // is greater than zero
    if (binTotal>0)
    {
      vtkSmartPointer<vtkDataArray> nextArray;
      double factor = 1.0/binTotal;
      for(int i=0; i<pd->GetNumberOfArrays(); ++i)
      {
        // *3 because TOTAL, AVERAGE, CUMULATIVE
        // ProfileElement &profileT = this->ProfileQuantities[PointDataArrays+i*3+0];
        ProfileElement &profileA = this->ProfileQuantities[PointDataArrays+i*3+1];
        // ProfileElement &profileC = this->ProfileQuantities[PointDataArrays+i*3+2];
        //
        this->UpdateBin(binNum, MULTIPLY, profileA, factor);
      }

      // For each additional quantity we also divide by N in the average
      // if requested
      for(int i=0; i<this->AdditionalProfileQuantities.size(); ++i)
      {
        ProfileElement &profile = this->AdditionalProfileQuantities[i];
        if (profile.ProfileColumnType==AVERAGE)
        {
          this->UpdateBin(binNum, MULTIPLY, profile, factor);
        }
      }
    }

    // Finally post processing those items which are marked as such, don't
    // do this within the average loop as these quantities are allowed
    // to depend on averages themselves

    for (int i=0; i<this->AdditionalProfileQuantities.size(); ++i)
    {
      ProfileElement &profile = this->AdditionalProfileQuantities[i];
      if (profile.FuntionType==POSTPROCESS_TYPE)
      {
        ProfileElement::PIterator profile1 = profile.ArgOne;
        ProfileElement::PIterator profile2 = profile.ArgTwo;
        double arg1[3], arg2[3], newval[3];
        (*profile.ArgOne).Data->GetTuple(binNum, arg1);
        (*profile.ArgTwo).Data->GetTuple(binNum, arg2);
        profile.funcP(arg1,arg2, newval);     
        this->UpdateBin(binNum, SET, profile, newval);
      }
    }
  }
}

//----------------------------------------------------------------------------
std::string vtkProfileFilter::GetColumnNameSlow(const char *baseName, 
  ColumnType columnType)
{
  switch(columnType)
    {
    case AVERAGE:
      return std::string(baseName)+"_average";
    case TOTAL:
      return std::string(baseName)+"_total";
    case CUMULATIVE:
      return std::string(baseName)+"_cumulative";
    default:
      vtkDebugMacro("columnType not found for function GetColumnName, returning error string");
      return "error";
    }
}

//----------------------------------------------------------------------------
void vtkProfileFilter::MergeBins(int binNum, BinUpdateType updateType,
   const char *baseName, ColumnType columnType, vtkTable* originalTable,
  vtkTable* tableToMerge)
{
/*
  vtkVariant originalData = this->GetData(binNum,baseName,columnType,
    originalTable);
  vtkVariant mergeData = this->GetData(binNum,baseName,columnType,
    tableToMerge);
  if(originalData.IsArray() && mergeData.IsArray())
    {
    this->UpdateArrayBin(binNum,updateType,baseName,columnType,
      mergeData.ToArray(),originalData.ToArray(),originalTable);
    }
  else
    {
    assert(originalData.IsNumeric() && mergeData.IsNumeric());
    this->UpdateDoubleBin(binNum,updateType,baseName,columnType,
      mergeData.ToDouble(),originalData.ToDouble(),originalTable);
    }
*/
}

//----------------------------------------------------------------------------
template<typename T>
void vtkProfileFilter::UpdateBin(int binNum, BinUpdateType updateType,
   ProfileElement &profile, T *updateData)
{
  float *tableData = &profile.DataPointer[binNum*profile.NumberComponents];
  UpdateBinN<T>(updateType, profile, tableData, updateData);
}
//----------------------------------------------------------------------------
void vtkProfileFilter::UpdateBin(int binNum, BinUpdateType updateType,
   ProfileElement &profile, double updateData)
{
  float *tableData = &profile.DataPointer[binNum*profile.NumberComponents];
  UpdateBin1<double>(updateType, profile, tableData, updateData);
}
//----------------------------------------------------------------------------
template<typename T>
void vtkProfileFilter::UpdateBinN(BinUpdateType updateType,
   ProfileElement &profile, float *tableData, T *newData)
{
  for (int i=0; i<profile.NumberComponents; i++) {
    switch (updateType) {
      case ADD:
        tableData[i] += newData[i];
        break;
      case MULTIPLY:
        tableData[i] *= newData[i];
        break;
      case SET:
        tableData[i] = newData[i];
        break;
    }
  }
}

//----------------------------------------------------------------------------
template<typename T>
void vtkProfileFilter::UpdateBin1(BinUpdateType updateType,
   ProfileElement &profile, float *tableData, T newData)
{
  for (int i=0; i<profile.NumberComponents; i++) {
    switch (updateType) {
      case ADD:
        tableData[i] += newData;
        break;
      case MULTIPLY:
        tableData[i] *= newData;
        break;
      case SET:
        tableData[i] = newData;
        break;
    }
  }
}

//----------------------------------------------------------------------------
void vtkProfileFilter::UpdateCumulativeBins(int binNum, BinUpdateType     
  updateType, ProfileElement &profile, double* dataToAdd)
{
  for(int bin = binNum; bin < this->BinNumber; ++bin)
    {
    this->UpdateBin(bin,updateType,profile,dataToAdd);
    }
}

//----------------------------------------------------------------------------
void vtkProfileFilter::UpdateCumulativeBins(int binNum, BinUpdateType     
  updateType, ProfileElement &profile, double dataToAdd)
{
  for(int bin = binNum; bin < this->BinNumber; ++bin)
    {
    this->UpdateBin(bin,updateType,profile,dataToAdd);
    }
}

//----------------------------------------------------------------------------
vtkProfileFilter::ProfileElement::ProfileElement(const char *baseName, 
  int numberComponents, DoubleFunction0 funcPtr,
  ColumnType columnType)
{
  this->BaseName          = baseName;
  this->NumberComponents  = numberComponents;
  this->ProfileColumnType = columnType;
  this->func0             = funcPtr;
  this->FuntionType       = FUNC0_TYPE;  
  this->CreateColumnNames();;
}

//----------------------------------------------------------------------------
vtkProfileFilter::ProfileElement::ProfileElement(const char *baseName, 
  int numberComponents, DoubleFunction1 funcPtr,
  ColumnType columnType)
{
  this->BaseName          = baseName;
  this->NumberComponents  = numberComponents;
  this->ProfileColumnType = columnType;
  this->func1             = funcPtr;
  this->FuntionType       = FUNC1_TYPE;  
  this->CreateColumnNames();;
}

//----------------------------------------------------------------------------
vtkProfileFilter::ProfileElement::ProfileElement(const char *baseName, 
  int numberComponents, DoubleFunction2 funcPtr,
  ColumnType columnType)
{
  this->BaseName          = baseName;
  this->NumberComponents  = numberComponents;
  this->ProfileColumnType = columnType;
  this->func2             = funcPtr;
  this->FuntionType       = FUNC2_TYPE;  
  this->CreateColumnNames();;
}

//----------------------------------------------------------------------------
vtkProfileFilter::ProfileElement::ProfileElement(const char *baseName, 
  int numberComponents, ColumnType columnType)
{
  this->BaseName          = baseName;
  this->NumberComponents  = numberComponents;
  this->ProfileColumnType = columnType;
  this->FuntionType       = NULL_FUNC;
  this->CreateColumnNames();;
}

//----------------------------------------------------------------------------
vtkProfileFilter::ProfileElement::ProfileElement(const char *baseName, 
	int numberComponents, DoubleFunction0 funcPtr,
			PIterator argOne, PIterator argTwo)
{
  this->BaseName          = baseName;
  this->NumberComponents  = numberComponents;
  this->ProfileColumnType = TOTAL;
  this->funcP             = funcPtr;
  this->FuntionType       = POSTPROCESS_TYPE;
  this->ArgOne            = argOne;
  this->ArgTwo            = argTwo;
  this->CreateColumnNames();
}

//----------------------------------------------------------------------------
vtkProfileFilter::ProfileElement::~ProfileElement()
{
  
}
//----------------------------------------------------------------------------
void vtkProfileFilter::ProfileElement::CreateColumnNames()
{
  switch(this->ProfileColumnType)
    {
    case AVERAGE:
      this->DerivedName = this->BaseName + "_average";
      break;
    case TOTAL:
      this->DerivedName = this->BaseName + "_total";
      break;
    case CUMULATIVE:
      this->DerivedName = this->BaseName + "_cumulative";
      break;
    default:
      vtkGenericWarningMacro("columnType not found for function CreateColumnNames, returning error string");
      this->DerivedName = "error";
    }
}
//----------------------------------------------------------------------------
const char *vtkProfileFilter::ProfileElement::GetColumnName() const
{
  return this->DerivedName.c_str();
}
//----------------------------------------------------------------------------
vtkFloatArray *vtkProfileFilter::ProfileElement::AllocateBinArray(vtkIdType numtuples)
{
	this->Data = vtkSmartPointer<vtkFloatArray>::New();
	this->Data->SetNumberOfComponents(this->NumberComponents);
	this->Data->SetNumberOfTuples(numtuples);
	this->Data->SetName(this->GetColumnName());
	this->DataPointer = this->Data->GetPointer(0);
	//initializes everything to zero
	for (vtkIdType i=0; i<this->NumberComponents*numtuples; ++i)
		{
		this->DataPointer[i] = 0.0;
		}
  return this->Data;
}
