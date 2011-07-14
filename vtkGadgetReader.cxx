/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkGadgetReader.cxx,v $
  Author: Christine Corbett Moran
=========================================================================*/
#include "vtkGadgetReader.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkDoubleArray.h" 
#include "vtkIntArray.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkDistributedDataFilter.h"
#include "vtkMultiProcessController.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkSmartPointer.h"
#include "vtkDataArraySelection.h"
#include "vtkMultiProcessController.h"
#include <cmath>
#include <assert.h>
#include <string>
#include <vector>
#include <cstdio>
#include <cfloat>
#include <iostream>
#include <functional>
#include <algorithm>
#include "assert.h"
#include "gadget/utility.h"
vtkCxxRevisionMacro(vtkGadgetReader, "$Revision: 1.0 $");
vtkStandardNewMacro(vtkGadgetReader);





/*=============================================================================
 *                              COMMON DEFINES
 *=============================================================================*/

double       GADGET_LUNIT;
double       GADGET_MUNIT;
int          LGADGET;
int          DGADGET;
unsigned int blklen;

long         IDmin =  4294967296;
long         IDmax = -4294967296; //simply for debugging purposes

#define MAXSTRING          2048 
#define GADGET_SKIP        ReadUInt(icfile,&blklen,this->Swap);
#define SIZEOFGADGETHEADER 256
#define MZERO             (1e-10)
#define X                  0
#define Y                  1
#define Z                  2


/*=============================================================================
 *                                STRUCTURES
 *=============================================================================*/
struct info_gadget
{
  int      no_gadget_files;
  int      i_gadget_file;
  long   *(np[6]);
  long     nall;
  
  struct io_gadget_header
  {
    int      np[6];
    double   massarr[6];
    double   expansion;
    double   redshift;
    int      flagsfr;
    int      flagfeedback;
    int      nall[6];
    int      flagcooling;
    int      NumFiles;
    double   BoxSize;
    double   Omega0;
    double   OmegaLambda;
    double   HubbleParam;
    char     unused[SIZEOFGADGETHEADER- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];
  } header; 
  
} gadget;

struct particle_data 
{
  double     Pos[3];       /* particle position   */  
  double     Vel[3];       /* particle velocity   */  
  double     Mass;         /* particle mass       */
  double     u;            /* gas internal energy */
  long       ID;           /* particle IDs        */
} *Part;


//----------------------------------------------------------------------------
vtkSmartPointer<vtkDoubleArray> AllocateGadgetDataArray(
  vtkDataSet *output, const char* arrayName, int numComponents, unsigned long numTuples)
{
  vtkSmartPointer<vtkDoubleArray> dataArray=vtkSmartPointer<vtkDoubleArray>::New();
	dataArray->SetNumberOfComponents(numComponents);
	dataArray->SetNumberOfTuples(numTuples);
	dataArray->SetName(arrayName);
	// initializes everything to zero
	for (int i=0; i < numComponents; ++i) {
		dataArray->FillComponent(i, 0.0);
  }
  output->GetPointData()->AddArray(dataArray);
  return dataArray;
}

//----------------------------------------------------------------------------
vtkGadgetReader::vtkGadgetReader()
{
	srand((unsigned)time(0));
  this->FileName          = 0;
  this->UpdatePiece       = 0;
  this->UpdateNumPieces   = 0;
  this->SetNumberOfInputPorts(0); 
  //
  this->PointDataArraySelection  = vtkDataArraySelection::New();
  //
  this->Positions     = NULL;
  this->Vertices      = NULL;
  this->GlobalIds     = NULL;
  this->ParticleIndex = 0;
  this->Potential   = NULL;
  this->Mass        = NULL;
  this->EPS         = NULL;
  this->RHO         = NULL;
  this->Hsmooth     = NULL;
  this->Temperature = NULL;
  this->Metals      = NULL;
  this->Tform       = NULL;
  this->Velocity    = NULL;
  this->Controller = NULL;
  this->Controller=vtkMultiProcessController::GetGlobalController();
  
}

//----------------------------------------------------------------------------
vtkGadgetReader::~vtkGadgetReader()
{
  this->SetFileName(0);
  this->PointDataArraySelection->Delete();
}

//----------------------------------------------------------------------------
void vtkGadgetReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "FileName: "
     << (this->FileName ? this->FileName : "(none)") << "\n";
}

		
//----------------------------------------------------------------------------
void vtkGadgetReader::AllocateAllGadgetVariableArrays(vtkIdType numBodies,
	vtkPolyData* output)
{
  // Allocate objects to hold points and vertex cells. 
  this->Positions = vtkSmartPointer<vtkPoints>::New();
  this->Positions->SetDataTypeToFloat();
  this->Positions->SetNumberOfPoints(numBodies);
  //
  this->Vertices  = vtkSmartPointer<vtkCellArray>::New();
  vtkIdType *cells = this->Vertices->WritePointer(numBodies, numBodies*2);
  for (vtkIdType i=0; i<numBodies; ++i) {
    cells[i*2]   = 1;
    cells[i*2+1] = i;
  }

  //
  this->GlobalIds = vtkSmartPointer<vtkIdTypeArray>::New();
  this->GlobalIds->SetName("global_id");
  this->GlobalIds->SetNumberOfTuples(numBodies);

 // Storing the points and cells in the output data object.
  output->SetPoints(this->Positions);
  output->SetVerts(this->Vertices); 

  // allocate velocity first as it uses the most memory and on my win32 machine 
  // this helps load really big data without alloc failures.
  if (this->GetPointArrayStatus("Velocity")) 
    this->Velocity = AllocateGadgetDataArray(output,"velocity",3,numBodies);
  else 
    this->Velocity = NULL;
  if (this->GetPointArrayStatus("Potential")) 
    this->Potential = AllocateGadgetDataArray(output,"potential",1,numBodies);
  else 
    this->Potential = NULL;
  if (this->GetPointArrayStatus("Mass"))
    this->Mass = AllocateGadgetDataArray(output,"mass",1,numBodies);
  else 
    this->Mass = NULL;
  if (this->GetPointArrayStatus("Eps")) 
    this->EPS = AllocateGadgetDataArray(output,"eps",1,numBodies);
  else 
    this->EPS = NULL;
  if (this->GetPointArrayStatus("Rho")) 
    this->RHO = AllocateGadgetDataArray(output,"rho",1,numBodies);
  else 
    this->RHO = NULL;
  if (this->GetPointArrayStatus("Hsmooth")) 
    this->Hsmooth = AllocateGadgetDataArray(output,"hsmooth",1,numBodies);
  else 
    this->Hsmooth = NULL;
  if (this->GetPointArrayStatus("Temperature"))
    this->Temperature = AllocateGadgetDataArray(output,"temperature",1,numBodies);
  else 
    this->Temperature = NULL;

  if (this->GetPointArrayStatus("Metals"))
    this->Metals = AllocateGadgetDataArray(output,"metals",1,numBodies);
  else 
    this->Metals = NULL;

	
	if (this->GetPointArrayStatus("Age"))
    this->Age = AllocateGadgetDataArray(output,"age",1,numBodies);
  else 
    this->Age = NULL;
	
	if (this->GetPointArrayStatus("Type"))
    this->Type = AllocateGadgetDataArray(output,"type",1,numBodies);
  else 
    this->Type = NULL;
	
  if (this->GetPointArrayStatus("Tform"))
    this->Tform = AllocateGadgetDataArray(output,"tform",1,numBodies);
  else 
    this->Tform = NULL;
}
//----------------------------------------------------------------------------
int vtkGadgetReader::RequestInformation(
	vtkInformation* vtkNotUsed(request),
	vtkInformationVector** vtkNotUsed(inputVector),
	vtkInformationVector* outputVector)
{
	vtkInformation* outInfo = outputVector->GetInformationObject(0);
	// means that the data set can be divided into an arbitrary number of pieces
	outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(),
		-1);

  this->PointDataArraySelection->AddArray("Potential");
  this->PointDataArraySelection->AddArray("Mass");
  this->PointDataArraySelection->AddArray("Eps");
  this->PointDataArraySelection->AddArray("Rho");
  this->PointDataArraySelection->AddArray("Hsmooth");
  this->PointDataArraySelection->AddArray("Temperature");
  this->PointDataArraySelection->AddArray("Metals");
	this->PointDataArraySelection->AddArray("Age");
  this->PointDataArraySelection->AddArray("Tform");
	this->PointDataArraySelection->AddArray("Type");
  this->PointDataArraySelection->AddArray("Velocity");

	return 1;
}
/*
* Reads a file, optionally only the marked particles from the file, 
* in the following order:
* 1. Open Gadget binary
* 2. Read Gadget header (tells us the number of particles of each type we are 
*    dealing with)
* NOTE: steps 3, 5 are currently not parallel
* 3. Read mark file indices from marked particle file, if there is one
* 4. Read either marked particles only or all particles
* 5. If an attribute file is additionally specified, reads this additional
* 	 attribute into a data array, reading only those marked if necessary.
*/
//----------------------------------------------------------------------------
int vtkGadgetReader::RequestData(vtkInformation*,
	vtkInformationVector**,vtkInformationVector* outputVector)
{
  FILE   *icfile;
  int     ipart;
  //
	// Make sure we have a file to read.
  //
  if(!this->FileName)
	  {
    vtkErrorMacro("A FileName must be specified.");
    return 0;
    }
  
  // Get output information
	vtkInformation* outInfo = outputVector->GetInformationObject(0);
  // get this->UpdatePiece information
  this->UpdatePiece = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
	this->UpdateNumPieces =outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());

  // get the output polydata
  vtkPolyData *output = \
    vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkSmartPointer<vtkPolyData> GadgetReadInitialOutput = \
  vtkSmartPointer<vtkPolyData>::New();  
	int mympirank=vtkMultiProcessController::GetGlobalController()->GetLocalProcessId();


	// TODO: Open the Gadget standard file and abort if there is an error.
  
  if((icfile = fopen(this->FileName,"rb")) == NULL) 
    {
    // TODO: maybe there are multiple GADGET files
      vtkErrorMacro("maybe there are multiple input GADGET files but support for this has not yet been implemented");
    }
  else 
    {

      gadget.no_gadget_files  = 1;
      gadget.i_gadget_file    = 0;
      
      /* allocate temporary storage for no. of particles arrays */
      gadget.np[0]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
      gadget.np[1]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
      gadget.np[2]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
      gadget.np[3]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
      gadget.np[4]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
      gadget.np[5]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
      
      this->ReadGadget(icfile,output);
      
      fclose(icfile);

    
    }

 	return 1;
}


// actually does the heavy lifting in terms of reading in from the file
int vtkGadgetReader::ReadGadget(FILE *icfile, vtkPolyData *output) {
  //--------------------------------
	
  long unsigned  ipart;
  
  double         tot_mass[6];
  
  int            i,j,k;
  int            no_part;
  int            massflag;
  char           DATA[MAXSTRING];
  float          fdummy[3];
  double         ddummy[3];
  double         x_fac, v_fac, m_fac;
  long           pid, ldummy;
  int            idummy;
  
  /*================= read in GADGET IO header =================*/
  if(this->Format == 1)
  {
    GADGET_SKIP;
    
    fread(DATA,sizeof(char),4,icfile);
    DATA[4] = '\0';
    GADGET_SKIP;
    GADGET_SKIP;
    fprintf(stderr,"reading... %s\n",DATA);
  }
  
  GADGET_SKIP;
  
  ReadInt(icfile,&(gadget.header.np[0]),this->Swap);
  ReadInt(icfile,&(gadget.header.np[1]),this->Swap);
  ReadInt(icfile,&(gadget.header.np[2]),this->Swap);
  ReadInt(icfile,&(gadget.header.np[3]),this->Swap);    /* number of particles in current file */
  ReadInt(icfile,&(gadget.header.np[4]),this->Swap);
  ReadInt(icfile,&(gadget.header.np[5]),this->Swap);
  ReadDouble(icfile,&(gadget.header.massarr[0]),this->Swap);
  ReadDouble(icfile,&(gadget.header.massarr[1]),this->Swap);
  ReadDouble(icfile,&(gadget.header.massarr[2]),this->Swap);
  ReadDouble(icfile,&(gadget.header.massarr[3]),this->Swap);
  ReadDouble(icfile,&(gadget.header.massarr[4]),this->Swap);
  ReadDouble(icfile,&(gadget.header.massarr[5]),this->Swap);
  ReadDouble(icfile,&(gadget.header.expansion),this->Swap);
  ReadDouble(icfile,&(gadget.header.redshift),this->Swap);
  ReadInt(icfile,&(gadget.header.flagsfr),this->Swap);
  ReadInt(icfile,&(gadget.header.flagfeedback),this->Swap);
  ReadInt(icfile,&(gadget.header.nall[0]),this->Swap);
  ReadInt(icfile,&(gadget.header.nall[1]),this->Swap);
  ReadInt(icfile,&(gadget.header.nall[2]),this->Swap);  /* total number of particles in simulation */
  ReadInt(icfile,&(gadget.header.nall[3]),this->Swap);
  ReadInt(icfile,&(gadget.header.nall[4]),this->Swap);
  ReadInt(icfile,&(gadget.header.nall[5]),this->Swap);
  ReadInt(icfile,&(gadget.header.flagcooling),this->Swap);
  ReadInt(icfile,&(gadget.header.NumFiles),this->Swap);
  ReadDouble(icfile,&(gadget.header.BoxSize),this->Swap);
  ReadDouble(icfile,&(gadget.header.Omega0),this->Swap);
  ReadDouble(icfile,&(gadget.header.OmegaLambda),this->Swap);
  ReadDouble(icfile,&(gadget.header.HubbleParam),this->Swap);
  ReadChars(icfile,&(gadget.header.unused[0]),SIZEOFGADGETHEADER- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8);
  
  GADGET_SKIP;
  /*================= read in GADGET IO header =================*/

  
  
  // here's where the ParaView specific code comes in
	
	int GADGET_SIZE = 1234; // TODO: change to read from file of course
	// Allocate the arrays
  // TODO Size from the gadget header
	this->AllocateAllGadgetVariableArrays(GADGET_SIZE, output);
	// Loop through and add to PV arrays
	double pos[3];
	double vel[3];
	double den[3] = {0.0,0.0,0.0};
	double pmass;
	
	// TODO: ignoring age and metals for now.
	for(unsigned i=0;i< GADGET_SIZE;i++) {
		pos[0]=0;
		pos[1]=0;
		pos[2]=0;
		vel[0]=0;
		vel[1]=0;
		vel[2]=0;
		// all particle types have this
		this->Positions->SetPoint(i, pos);
		if (this->Velocity)  this->Velocity->SetTuple(i, vel);
		if (this->Mass)      this->Mass->SetTuple1(i, 0.0);
		if (this->Type)      this->Type->SetTuple1(i, 0.0);
    
    if (this->Metals)      this->Metals->SetTuple1(i, 0.0);
    if (this->Age)      this->Age->SetTuple1(i, 0.0);
    
		// These currently are unused, should be removed
		if (this->Potential) this->Potential->SetTuple1(i,0.0);		
	  if (this->RHO)         this->RHO->SetTuple(i, den);
		if (this->Temperature) this->Temperature->SetTuple1(i, 0.0);
		if (this->Hsmooth)     this->Hsmooth->SetTuple1(i, 0.0);
		if (this->EPS)    this->EPS->SetTuple1(i, 0.0);
    
	}
	
	
	// Done, vis o'clock
	// can free memory allocated in vectors above
	
  vtkDebugMacro("Reading all points from file " << this->FileName);
  // Read Successfully
  vtkDebugMacro("Read " << output->GetPoints()->GetNumberOfPoints() \
                << " points.");
  // release memory smartpointers - just to play safe.
  this->Vertices    = NULL;
  this->GlobalIds   = NULL;
  this->Positions   = NULL;
  this->Potential   = NULL;
  this->Mass        = NULL;
  this->EPS         = NULL;
  this->RHO         = NULL;
  this->Hsmooth     = NULL;
  this->Temperature = NULL;
  this->Metals      = NULL;
	this->Age         = NULL;
  this->Tform       = NULL;
	this->Type        = NULL;
  this->Velocity    = NULL;
  //

  return 0;
}

//----------------------------------------------------------------------------
// Below : Boiler plate code to handle selection of point arrays
//----------------------------------------------------------------------------
const char* vtkGadgetReader::GetPointArrayName(int index) {
  return this->PointDataArraySelection->GetArrayName(index);
}
//----------------------------------------------------------------------------
int vtkGadgetReader::GetPointArrayStatus(const char* name)
{
  return this->PointDataArraySelection->ArrayIsEnabled(name);
}
//----------------------------------------------------------------------------
void vtkGadgetReader::SetPointArrayStatus(const char* name, int status)
{
  if (status!=this->GetPointArrayStatus(name)) {
    if (status) {
      this->PointDataArraySelection->EnableArray(name);
    }
    else {
      this->PointDataArraySelection->DisableArray(name);
    }
    this->Modified();
  }
}
//----------------------------------------------------------------------------
void vtkGadgetReader::Enable(const char* name)
{
  this->SetPointArrayStatus(name, 1);
}
//----------------------------------------------------------------------------
void vtkGadgetReader::Disable(const char* name)
{
  this->SetPointArrayStatus(name, 0);
}
//----------------------------------------------------------------------------
void vtkGadgetReader::EnableAll()
{
  this->PointDataArraySelection->EnableAllArrays();
}
//----------------------------------------------------------------------------
void vtkGadgetReader::DisableAll()
{
  this->PointDataArraySelection->DisableAllArrays();
}
//----------------------------------------------------------------------------
int vtkGadgetReader::GetNumberOfPointArrays()
{
  return this->PointDataArraySelection->GetNumberOfArrays();
}
