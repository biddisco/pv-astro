/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkVirialRadiusFilter.cxx,v $
=========================================================================*/
#include "vtkVirialRadiusFilter.h"
#include "AstroVizHelpersLib/AstroVizHelpers.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkCellData.h"
#include "vtkSortDataArray.h"
#include "vtkDataSet.h"
#include "vtkMath.h"
#include "vtkInformationDataObjectKey.h"
#include "vtkPointSet.h" 
#include "vtkPointLocator.h"
#include "vtkMultiProcessController.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPolyData.h"
#include <cmath>
using vtkstd::string;

vtkCxxRevisionMacro(vtkVirialRadiusFilter, "$Revision: 1.72 $");
vtkStandardNewMacro(vtkVirialRadiusFilter);
vtkCxxSetObjectMacro(vtkVirialRadiusFilter,Controller,
	vtkMultiProcessController);
	
//----------------------------------------------------------------------------
vtkVirialRadiusFilter::vtkVirialRadiusFilter()
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
	this->Delta=0.0;
	this->Controller = NULL;
	this->SetController(vtkMultiProcessController::GetGlobalController());
}

//----------------------------------------------------------------------------
vtkVirialRadiusFilter::~vtkVirialRadiusFilter()
{
	this->SetController(NULL);
}

//----------------------------------------------------------------------------
void vtkVirialRadiusFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  os << indent << "overdensity: " << this->Delta << "\n"
		<< "softening :" << this->Softening << "\n";
}

//----------------------------------------------------------------------------
void vtkVirialRadiusFilter::SetSourceConnection(vtkAlgorithmOutput* algOutput)
{
  this->SetInputConnection(1, algOutput);
}

//----------------------------------------------------------------------------
int vtkVirialRadiusFilter::FillInputPortInformation (int port, 
	vtkInformation *info)
{
  this->Superclass::FillInputPortInformation(port, info);
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  return 1;
}

//----------------------------------------------------------------------------
int vtkVirialRadiusFilter::FillOutputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  // now add our info
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
  return 1;
}

//----------------------------------------------------------------------------
int vtkVirialRadiusFilter::RequestData(vtkInformation *request,
																	vtkInformationVector **inputVector,
																	vtkInformationVector *outputVector)
{
	// Now we can get the input with which we want to work
 	vtkPointSet* input = vtkPointSet::GetData(inputVector[0]);
	// Setting the center based upon the selection in the GUI
	vtkDataSet* pointInfo = vtkDataSet::GetData(inputVector[1]);
	// Output
	vtkUnstructuredGrid* output = vtkUnstructuredGrid::GetData(outputVector);
 	output->ShallowCopy(input);
	// Get name of data array containing mass
	vtkDataArray* massArray = this->GetInputArrayToProcess(0, inputVector);
  if (!massArray)
    {
    vtkErrorMacro("Failed to locate mass array");
    return 0;
    }
	this->CalculateAndSetBounds(output,pointInfo);
	
	// Building the point locator and the struct to use as an 
	// input to the rootfinder.
	// 1. Building the point locator
	vtkPointLocator* locator = vtkPointLocator::New();
		locator->SetDataSet(output);
		locator->BuildLocator();
	// Will communicate with other processes if necessary
	VirialRadiusInfo virialRadiusInfo = \
	 	ComputeVirialRadius(this->GetController(),
		locator,massArray->GetName(),this->Softening,
		this->Delta,this->MaxR,this->Center);	
	// note that if there was an error finding the virialRadius the 
	// radius returned is < 0
	if(virialRadiusInfo.virialRadius>0)
		{
		// This is causing problems, segfault
		//setting the dataSet to this newInput
		vtkPolyData* newDataSet = \
			GetDatasetWithinVirialRadius(virialRadiusInfo);	
		// resetting output, then copying into it the new data set
		output->Initialize();
		output->DeepCopy(newDataSet);
		output->SetCells(VTK_VERTEX,newDataSet->GetVerts());
		newDataSet->Delete();
		}
	else	
		{
		vtkErrorMacro("Unable to find virial radius: considering changing your delta or selecting a different point around which to search. For now simply copying input");
		}
	return 1;	
}

//----------------------------------------------------------------------------
void vtkVirialRadiusFilter::CalculateAndSetBounds(vtkPointSet* input, 
	vtkDataSet* source)
{
	if(RunInParallel(this->GetController()))
		{
		int procId=this->GetController()->GetLocalProcessId();
		int numProc=this->GetController()->GetNumberOfProcesses();
		if(procId==0)
			{
      source->GetCenter(this->Center);
			// Syncronizing the centers
			this->GetController()->Broadcast(this->Center,3,0);			
			}
		else
			{
			// Syncronizing the centers
			this->GetController()->Broadcast(this->Center,3,0);
			}
		// calculating the max R
		this->MaxR=ComputeMaxRadiusInParallel(this->GetController(),
			input,this->Center);
		}
	else
		{
		// we aren't using MPI or have only one process
    source->GetCenter(this->Center);
		//calculating the the max R
		this->MaxR=ComputeMaxR(input,this->Center);			
		}
}

