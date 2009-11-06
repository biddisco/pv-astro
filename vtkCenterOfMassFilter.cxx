/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkCenterOfMassFilter.cxx,v $
=========================================================================*/
#include "vtkCenterOfMassFilter.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStringArray.h"
#include "vtkSphereSource.h"
#include "vtkMultiProcessController.h"
#include "astrovizhelpers/DataSetHelpers.h"
#include "astrovizhelpers/ProfileHelpers.h"


vtkCxxRevisionMacro(vtkCenterOfMassFilter, "$Revision: 1.72 $");
vtkStandardNewMacro(vtkCenterOfMassFilter);

//----------------------------------------------------------------------------
vtkCenterOfMassFilter::vtkCenterOfMassFilter()
{
	this->Overdensity = 0; 
	this->Softening=1e-6f;
}

//----------------------------------------------------------------------------
vtkCenterOfMassFilter::~vtkCenterOfMassFilter()
{
}

//----------------------------------------------------------------------------
void vtkCenterOfMassFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
int vtkCenterOfMassFilter::FillInputPortInformation(int, vtkInformation* info)
{
  // This filter uses the vtkDataSet cell traversal methods so it
  // suppors any data set type as input.
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  return 1;
}

//----------------------------------------------------------------------------
int vtkCenterOfMassFilter::RequestData(vtkInformation*,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector)
{
  // Get input and output data.
  vtkPointSet* input = vtkPointSet::GetData(inputVector[0]);
  vtkPolyData* output = vtkPolyData::GetData(outputVector);
	// TODO: playing around with Parallel PV, checking if serial or parallel
	if (vtkMultiProcessController::GetGlobalController() == NULL)
		{
		cout << "\nSERIAL\n";
		}
	else
		{
		cout << "\nPARALLEL\n";
		}
	// we will create one point in the output: the center of mass point
	output->SetPoints(vtkSmartPointer<vtkPoints>::New());
	output->SetVerts(vtkSmartPointer<vtkCellArray>::New()); 
	double* dbCenterOfMass = ComputeCOM(input);
	float* centerOfMass = new float[3];
	for(int i = 0; i < 3; ++i)
		{
		centerOfMass[i]=static_cast<float>(dbCenterOfMass[i]);
		}
	// if the Overdensity is non zero and we are able to find a
	// virial radius then we set the output to the sphere
	// around the COM at the virial radius.
	if(this->Overdensity>0)
		{
			double maxR=ComputeMaxR(input,dbCenterOfMass);
			VirialRadiusInfo virialRadiusInfo=\
			ComputeVirialRadius(input,this->Softening,
				this->Overdensity,maxR,dbCenterOfMass);
			if(virialRadiusInfo.virialRadius>0)
				{
				//Here is where we create the sphere around the COM to display
				vtkWarningMacro("the virial radius is " 
												<< virialRadiusInfo.virialRadius);
				// Creating the sphere
				CreateSphere(output,\
										virialRadiusInfo.virialRadius,dbCenterOfMass);
				}
			else
				{
				vtkWarningMacro("unable to find the virial radius from over density you specified. Perhaps it is too high. For now displaying only the center of mass");
				// Placing the point's data in the output
				SetPointValue(output,centerOfMass);					
				}
		}
	else
		{
			// Placing the point's data in the output
			SetPointValue(output,centerOfMass);
		}
	// finally, some memory management
	delete [] dbCenterOfMass;
	delete [] centerOfMass;
  return 1;
}
