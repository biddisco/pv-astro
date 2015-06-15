/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkStrangeAttractors.cxx,v $
=========================================================================*/
#include "vtkStrangeAttractors.h"
#ifdef HAS_GMP
#include "gmp.h" // Requires http://gmplib.org/
#endif
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h" 
#include "vtkIntArray.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkDistributedDataFilter.h"
#include "vtkMultiProcessController.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkSmartPointer.h"
#include "vtkDataArraySelection.h"
#include "vtkCharArray.h"
#include <cmath>
#include <assert.h>
#include <sys/stat.h>
#include <iostream>
#include <deque>
vtkStandardNewMacro(vtkStrangeAttractors);


//----------------------------------------------------------------------------
vtkStrangeAttractors::vtkStrangeAttractors()
{
  this->FileName          = 0;
  this->SetNumberOfInputPorts(0); 
  this->Positions     = NULL;
  this->Vertices      = NULL;
  this->ParticleIndex = 0;
  this->PrimesOnly  = 0;
  this->Primes   = NULL;
  this->UpdatePiece       = 0;
  this->UpdateNumPieces   = 0;
}

//----------------------------------------------------------------------------
vtkStrangeAttractors::~vtkStrangeAttractors()
{
  this->SetFileName(0);
}

//----------------------------------------------------------------------------
void vtkStrangeAttractors::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "FileName: "
     << (this->FileName ? this->FileName : "(none)") << "\n";
}

//----------------------------------------------------------------------------
bool vtkStrangeAttractors::RunInParallel(vtkMultiProcessController* controller)
{
  return (controller != NULL && controller->GetNumberOfProcesses() > 1);
}


//----------------------------------------------------------------------------
void vtkStrangeAttractors::AllocateAllVariableArrays(vtkIdType numBodies,
	vtkPolyData* output)
{
  // Allocate objects to hold points and vertex cells. 
  this->Positions = vtkSmartPointer<vtkPoints>::New();
  this->Positions->SetDataTypeToFloat();
  this->Positions->SetNumberOfPoints(numBodies);
  //
  this->Vertices  = vtkSmartPointer<vtkCellArray>::New();
  vtkIdType *cells = this->Vertices->WritePointer(numBodies, numBodies*2);
  for(vtkIdType i=0; i<numBodies; ++i) {
    cells[i*2]   = 1;
    cells[i*2+1] = i;
  }

  // Storing the points and cells in the output data object.
  output->SetPoints(this->Positions);
  output->SetVerts(this->Vertices); 
#ifdef HAS_GMP
  this->Primes = vtkSmartPointer<vtkIntArray>::New();
  this->Primes->SetName("primes");
  this->Primes->SetNumberOfComponents(1);
  this->Primes->SetNumberOfTuples(numBodies);
  this->Primes->FillComponent(0, 0.0);
  output->GetPointData()->AddArray(this->Primes);
#endif		
  this->Velocity=vtkSmartPointer<vtkFloatArray>::New();
  this->Velocity->SetName("velocity");
  this->Velocity->SetNumberOfComponents(3);
  this->Velocity->SetNumberOfTuples(numBodies);
  for(int i=0; i<3; i++) {
    this->Velocity->FillComponent(i,0.0);
  }
  output->GetPointData()->AddArray(this->Velocity);
}

//----------------------------------------------------------------------------
int vtkStrangeAttractors::RequestInformation(
	vtkInformation* vtkNotUsed(request),
	vtkInformationVector** vtkNotUsed(inputVector),
	vtkInformationVector* outputVector)
{
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  // means that the data set can be divided into an arbitrary number of pieces
//  outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(),-1);
  return 1;
}

//----------------------------------------------------------------------------
int vtkStrangeAttractors::RequestData(vtkInformation*,
	vtkInformationVector**,vtkInformationVector* outputVector)
{
	// Make sure we have a file to read.
  if(!this->FileName) {
    vtkErrorMacro("A FileName must be specified.");
    return 0;
  }
  // Get output information
  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  // get the output polydata
  vtkPolyData *output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkSmartPointer<vtkPolyData> tipsyReadInitialOutput = vtkSmartPointer<vtkPolyData>::New();
    
  // get this->UpdatePiece information
  this->UpdatePiece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  this->UpdateNumPieces =outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
  
  // reset counter before reading
  this->ParticleIndex = 0;
  // Finding size of file
  FILE *filep = fopen(this->FileName,"r");
  fseek(filep, 0, SEEK_END);
  unsigned long long fileSize = ftell(filep);
  fseek(filep, 0, SEEK_SET);
	
  unsigned long pieceSize = floor(fileSize*1./this->UpdateNumPieces);
  unsigned long beginIndex = this->UpdatePiece*this->UpdateNumPieces;
  unsigned long endIndex = (this->UpdatePiece == this->UpdateNumPieces - 1) ? fileSize : (this->UpdatePiece+1)*pieceSize;

  // TODO: prime reading not supported in parallel, this will require a sync;
    
  vtkErrorMacro("Reading a file of size " << fileSize);
							
  // Counting the number of primes if we are going to be
  // only displaying primes. Useful for allocating arrays
  // TODO: this crashes for big files (bigger than int)
  // think I have to limit reading in to one portion of the file at a time
  std::vector<char> buf;
  buf.assign(fileSize,0);	
  fread(&buf[0],1,fileSize,filep); 
  // TODO: add back in
  // Removing for debugging /*
#ifdef HAS_GMP
  unsigned long long numLocalPrimes=0;
  unsigned long long numPrimes=0;
  mpz_t p;
  mpz_init(p);
  if(this->PrimesOnly) {
    // TODO: here is where I parallelize this
    for(unsigned long long i=beginIndex; i < endIndex; i++) {
      mpz_set_ui(p,buf[i]);
      int is_prime=mpz_probab_prime_p(p,5);			
      if(is_prime) {
	numLocalPrimes++;
      }
    }
    // An MPI sum of this
    if(RunInParaller(this->GetController()) {
      this->Controller->AllReduce(numLocalPrimes, numPrimes, 1, vtkCommunicator::SUM_OP); 	
    }
    else {
      numPrimes=numLocalPrimes;
    }
  }
  // Allocating arrays and preparing for read
  this->PrimesOnly ? this->AllocateAllVariableArrays(numPrimes+1,output) : \
    this->AllocateAllVariableArrays(fileSize,output);
  if(RunInParallel(this->GetController()) {
      vtkErrorMacro("warning-parallel is not supported in prime detection mmode at this point");
      // enforcing this
      beginIndex=0;
      endIndex=fileSize;
    }
#endif
#ifndef HAS_GMP
   this->AllocateAllVariableArrays(fileSize,output);
   vtkErrorMacro("Warning, prime detection is disabled.");
#endif	
    std::deque<char> delayedCoords(7,0);
    float pos[3];
    float vel[3];
    unsigned long long particleIdx=0;
    for(unsigned long long i=beginIndex; i<endIndex; i++) {
      unsigned long long is_prime=0;
#ifdef HAS_GMP
      // determining if buf[i] is prime
      mpz_set_ui(p,buf[i]);
      is_prime=mpz_probab_prime_p(p,5);
#endif
      if(!this->PrimesOnly || is_prime){
	if(!this->PrimesOnly || is_prime){

	  particleIdx++;
	  delayedCoords.push_front(buf[i]);
	  delayedCoords.pop_back();
	  // we only start computing coords once we have enough in our buffer
	  pos[0]=delayedCoords[3]-delayedCoords[2];
	  pos[1]=delayedCoords[2]-delayedCoords[1];
	  pos[2]=delayedCoords[1]-delayedCoords[0];			
	  
	  vel[0]=delayedCoords[6]-delayedCoords[5];
	  vel[1]=delayedCoords[5]-delayedCoords[4];
	  vel[2]=delayedCoords[4]-delayedCoords[3];
			
	  this->Primes->SetTuple1(particleIdx,is_prime);
	  this->Positions->SetPoint(particleIdx,pos);
	  this->Velocity->SetTuple(particleIdx,vel);
	}
      }
#ifdef HAS_GMP
      mpz_clear(p);
#endif

  vtkDebugMacro("Read " << output->GetPoints()->GetNumberOfPoints() \
		<< " points.");
  // release memory smartpointers - just to play safe.
  this->Vertices    = NULL;
  this->Primes   = NULL;
  this->Positions = NULL;
  this->Velocity   = NULL;
// */
 	return 1;
}
}

