/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkFriendsOfFriendsHaloFinder.cxx,v $
=========================================================================*/
#include "vtkFriendsOfFriendsHaloFinder.h"
#include "vtkIdList.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"
#include "vtkPKdTree.h"
#include "vtkPointSet.h"
#include "vtkPoints.h"
#include "vtkGenericPointIterator.h"
#include "vtkDataArray.h"
#include "vtkMath.h"
#include "vtkMultiProcessController.h"
#include "vtkPKdTree.h"
#include "vtkDistributedDataFilter.h"
#include "vtkCallbackCommand.h"
#include <vtkstd/vector>
#include <vtkstd/map>
#include "astrovizhelpers/DataSetHelpers.h"


vtkCxxRevisionMacro(vtkFriendsOfFriendsHaloFinder, "$Revision: 1.72 $");
vtkStandardNewMacro(vtkFriendsOfFriendsHaloFinder);
vtkCxxSetObjectMacro(vtkFriendsOfFriendsHaloFinder,Controller,
	vtkMultiProcessController);

//----------------------------------------------------------------------------
vtkFriendsOfFriendsHaloFinder::vtkFriendsOfFriendsHaloFinder()
{
  this->LinkingLength = 1e-6; //default
	this->MinimumNumberOfParticles = 50; // default
	this->Controller = NULL;
  this->SetController(vtkMultiProcessController::GetGlobalController());
}

//----------------------------------------------------------------------------
vtkFriendsOfFriendsHaloFinder::~vtkFriendsOfFriendsHaloFinder()
{
  this->SetController(NULL);
}

//----------------------------------------------------------------------------
void vtkFriendsOfFriendsHaloFinder::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "Linking Length: " << this->LinkingLength 
		<<	indent << "Minimum Number Of Particles: " 
		<<  this->MinimumNumberOfParticles << "\n";
}

//----------------------------------------------------------------------------
int vtkFriendsOfFriendsHaloFinder::FillInputPortInformation(int, vtkInformation* info)
{
  // This filter uses the vtkDataSet cell traversal methods so it
  // suppors any data set type as input.
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  return 1;
}

//----------------------------------------------------------------------------
int vtkFriendsOfFriendsHaloFinder::FindHaloes(vtkPointSet* input,
	vtkPointSet* output)
{
	if(this->MinimumNumberOfParticles < 2)
		{
		vtkWarningMacro("setting minimum number of particles to 2, a minimum number of particles below this makes no sense.");
		this->MinimumNumberOfParticles=2;
		}
	// Building the Kd tree
	vtkSmartPointer<vtkPKdTree> pointTree; 
	if(this->Controller!=NULL)
		{
		// Run D3
		
		// get PKdTree from D3's output
			
		}
	else
		{
		pointTree	= vtkSmartPointer<vtkPKdTree>::New();
				pointTree->BuildLocatorFromPoints(input);
		}
	// calculating the initial haloes- yes it really is done in just this 
	// one line.
	vtkSmartPointer<vtkIdTypeArray> haloIdArray = \
		pointTree->BuildMapForDuplicatePoints(this->LinkingLength);
	haloIdArray->SetNumberOfComponents(1);
	haloIdArray->SetNumberOfTuples(input->GetPoints()->GetNumberOfPoints());
	haloIdArray->SetName("halo ID");
	// Now assign halos, if this point has at least one other pair,
	// it is a halo, if not it is not (set to 0)
	// first building map of id to count of that id, O(N)
	vtkstd::map<vtkIdType,int> haloCount;
	int uniqueId=1;
	// use negatives to figure differentiate between unique id assignment 
	// (positive) and count (negative) in the same map
	for(int nextHaloId = 0;
		nextHaloId < haloIdArray->GetNumberOfTuples();
	 	++nextHaloId)
		{
		vtkIdType haloId = haloIdArray->GetValue(nextHaloId);
		if(haloCount[haloId]==-1*this->MinimumNumberOfParticles)
			{
			// we have seen the id minimum number of particles times,
			// so we assign it a unique halo id, considering it a halo
			haloCount[haloId]=uniqueId;
			uniqueId+=1;
			}
		else if(haloCount[haloId]<1)
			{
			// this counts each time we see the id, until we reach the
			// minimum number of particles to count as a halo
			haloCount[haloId]-=1;
			}
		}
	// finally setting to zero points which have 
	// count < this->MinimumNumberOfParticles, O(N), and
	// assigning those we have seen more the requisite number
	// of times to their unique id
	for(int nextHaloId = 0;
		nextHaloId < haloIdArray->GetNumberOfTuples();
	 	++nextHaloId)
		{
		vtkIdType haloId = haloIdArray->GetValue(nextHaloId);
		if(haloCount[haloId]<1)
			{
			// we only saw it less than requisite number of times
			haloIdArray->SetValue(nextHaloId,0);
			}
		else
			{
			// we saw it more than once, assign it to its unique id
			haloIdArray->SetValue(nextHaloId,haloCount[haloId]);
			}
		}
	output->GetPointData()->AddArray(haloIdArray);
}

//----------------------------------------------------------------------------
int vtkFriendsOfFriendsHaloFinder::RequestData(vtkInformation*,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector)
{
  // Get input and output data.
  vtkPointSet* input = vtkPointSet::GetData(inputVector[0]);
  vtkPointSet* output = vtkPointSet::GetData(outputVector);
  output->ShallowCopy(input);
	this->FindHaloes(input,output);
	// Finally, some memory management
  output->Squeeze();
  return 1;
}
