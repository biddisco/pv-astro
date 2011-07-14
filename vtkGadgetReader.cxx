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
vtkCxxRevisionMacro(vtkGadgetReader, "$Revision: 1.0 $");
vtkStandardNewMacro(vtkGadgetReader);



/*=============================================================================
 *                                TODO: move this out of file, linking errors on 64 bit OS X so doing this
 *=============================================================================*/

/*=============================================================================
 *                              COMMON DEFINES
 *=============================================================================*/

double       GADGET_LUNIT;
double       GADGET_MUNIT;
int          LGADGET;
int          DGADGET;
unsigned int blklen;

#define MAXSTRING          2048 
#define GADGET_SKIP        ReadUInt(icfile,&blklen,this->Swap);
#define SIZEOFGADGETHEADER 256
#define MZERO             (1e-10)
#define X                  0
#define Y                  1
#define Z                  2


enum GadgetParticleTypes 
{
	GADGET_GAS=0,
  GADGET_HALO=1,
  GADGET_DISK=2,
  GADGET_BULGE=3,
  GADGET_STARS=4,
  GADGET_BNDRY=5
};

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


/*=============================================================================
 *                                Reading functions
 *=============================================================================*/

/*
 Read a string of n characters
 */
int ReadString(FILE *fptr,char *s,int n)
{
  int i,c;
  
  if(sizeof(char) != 1)
  {
    fprintf(stderr,"ReadString: sizeof(char)=%ld and not 1\n",sizeof(char));
    exit(0);
  }
  s[0] = '\0';
  for (i=0;i<n;i++) {
    c = fgetc(fptr);
    if (c == EOF)
      return(false);
    s[i] = c;
    s[i+1] = '\0';
  }
  return(true);
}

/*
 Read an array of n characters
 NOTE: the difference to ReadString() is that we do not '\0'-terminate the array
 */
int ReadChars(FILE *fptr,char *s,int n)
{
  int i,c;
  
  if(sizeof(char) != 1)
  {
    fprintf(stderr,"ReadChars: sizeof(char)=%ld and not 1\n",sizeof(char));
    exit(0);
  }
  
  s[0] = '\0';
  for (i=0;i<n;i++) {
    c = fgetc(fptr);
    if (c == EOF)
      return(false);
    s[i] = c;
  }
  return(true);
}

/*
 Read a possibly byte swapped integer
 */
int ReadInt(FILE *fptr,int *n,int swap)
{
  unsigned char *cptr,tmp;
  
  if(sizeof(int) != 4)
  {
    fprintf(stderr,"ReadInt: sizeof(int)=%ld and not 4\n",sizeof(int));
    exit(0);
  }
  
  if (fread(n,4,1,fptr) != 1)
    return(false);
  if (swap) {
    cptr = (unsigned char *)n;
    tmp     = cptr[0];
    cptr[0] = cptr[3];
    cptr[3] = tmp;
    tmp     = cptr[1];
    cptr[1] = cptr[2];
    cptr[2] = tmp;
  }
  return(true);
}

/*
 Read a possibly byte swapped unsigned integer
 */
int ReadUInt(FILE *fptr,unsigned int *n,int swap)
{
  unsigned char *cptr,tmp;
  
  if(sizeof(int) != 4)
  {
    fprintf(stderr,"ReadInt: sizeof(int)=%ld and not 4\n",sizeof(int));
    exit(0);
  }
  
  if (fread(n,4,1,fptr) != 1)
    return(false);
  if (swap) {
    cptr = (unsigned char *)n;
    tmp     = cptr[0];
    cptr[0] = cptr[3];
    cptr[3] = tmp;
    tmp     = cptr[1];
    cptr[1] = cptr[2];
    cptr[2] = tmp;
  }
  return(true);
}

/*
 Read a possibly byte swapped long integer
 */
int ReadLong(FILE *fptr,long *n,int swap)
{
  unsigned char *cptr,tmp;
  
  if(sizeof(long) == 4)
  {
    if (fread(n,4,1,fptr) != 1)
      return(false);
    if (swap) {
      cptr = (unsigned char *)n;
      tmp     = cptr[0];
      cptr[0] = cptr[3];
      cptr[3] = tmp;
      tmp     = cptr[1];
      cptr[1] = cptr[2];
      cptr[2] = tmp;
    }
  }
  else if(sizeof(long) == 8)
  {
    if (fread(n,8,1,fptr) != 1)
      return(false);
    if (swap) {
      cptr = (unsigned char *)n;
      tmp     = cptr[0];
      cptr[0] = cptr[7];
      cptr[7] = tmp;
      tmp     = cptr[1];
      cptr[1] = cptr[6];
      cptr[6] = tmp;
      tmp     = cptr[2];
      cptr[2] = cptr[5];
      cptr[5] = tmp;
      tmp     = cptr[3];
      cptr[3] = cptr[4];
      cptr[4] = tmp;
    }
  }
  else
  {
    fprintf(stderr,"ReadLong: something wrong...cannot read long\n");
    exit(0);
  }
  
  
  
  return(true);
}

/*
 Read a possibly byte swapped long long integer
 */
int ReadLongLong(FILE *fptr,long long *n,int swap)
{
  unsigned char *cptr,tmp;
  
  if (fread(n,8,1,fptr) != 1)
    return(false);
  if (swap) {
    cptr = (unsigned char *)n;
    tmp     = cptr[0];
    cptr[0] = cptr[7];
    cptr[7] = tmp;
    tmp     = cptr[1];
    cptr[1] = cptr[6];
    cptr[6] = tmp;
    tmp     = cptr[2];
    cptr[2] = cptr[5];
    cptr[5] = tmp;
    tmp     = cptr[3];
    cptr[3] = cptr[4];
    cptr[4] = tmp;
  }
  return(true);
}

/*
 Read a possibly byte swapped double precision number
 Assume IEEE
 */
int ReadDouble(FILE *fptr,double *n,int swap)
{
  unsigned char *cptr,tmp;
  
  if(sizeof(double) != 8)
  {
    fprintf(stderr,"ReadDouble: sizeof(double)=%ld and not 8\n",sizeof(double));
    exit(0);
  }
  
  if (fread(n,8,1,fptr) != 1)
    return(false);
  if (swap) {
    cptr = (unsigned char *)n;
    tmp     = cptr[0];
    cptr[0] = cptr[7];
    cptr[7] = tmp;
    tmp     = cptr[1];
    cptr[1] = cptr[6];
    cptr[6] = tmp;
    tmp     = cptr[2];
    cptr[2] = cptr[5];
    cptr[5] = tmp;
    tmp     = cptr[3];
    cptr[3] = cptr[4];
    cptr[4] = tmp;
  }
  
  return(true);
}

/*
 Read a possibly byte swapped floating point number
 Assume IEEE format
 */
int ReadFloat(FILE *fptr,float *n, int swap)
{
  unsigned char *cptr,tmp;
  
  if(sizeof(float) != 4)
  {
    fprintf(stderr,"ReadFloat: sizeof(float)=%ld and not 4\n",sizeof(float));
    exit(0);
  }
  
  if (fread(n,4,1,fptr) != 1)
    return(false);
  if (swap) 
  {
    cptr = (unsigned char *)n;
    tmp     = cptr[0];
    cptr[0] = cptr[3];
    cptr[3] = tmp;
    tmp     = cptr[1];
    cptr[1] = cptr[2];
    cptr[2] = tmp;
  }
  return(true);
}

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
  this->FilePrefix        = 0;
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
  this->Mass        = NULL;
  this->Energy         = NULL;
  this->Velocity    = NULL;
  this->Controller = NULL;
  this->Controller=vtkMultiProcessController::GetGlobalController();
  
}

//----------------------------------------------------------------------------
vtkGadgetReader::~vtkGadgetReader()
{
  this->SetFileName(0);
  this->SetFilePrefix(0);
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
  if(this->Positions == NULL) {
    // We are creating arrays for the first time
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
    this->GlobalIds->SetNumberOfComponents(1);
    this->GlobalIds->SetNumberOfTuples(numBodies);
    this->GlobalIds->FillComponent(0, 0.0);
    
    output->GetPointData()->AddArray(this->GlobalIds);

   // Storing the points and cells in the output data object.
    output->SetPoints(this->Positions);
    output->SetVerts(this->Vertices); 

    // allocate velocity first as it uses the most memory and on my win32 machine 
    // this helps load really big data without alloc failures.
    if (this->GetPointArrayStatus("Velocity")) 
      this->Velocity = AllocateGadgetDataArray(output,"velocity",3,numBodies);
    else 
      this->Velocity = NULL;
    if (this->GetPointArrayStatus("Mass"))
      this->Mass = AllocateGadgetDataArray(output,"mass",1,numBodies);
    else 
      this->Mass = NULL;
    if (this->GetPointArrayStatus("Energy")) 
      this->Energy = AllocateGadgetDataArray(output,"Energy",1,numBodies);
    else 
      this->Energy = NULL;
    
    if (this->GetPointArrayStatus("Type"))
      this->Type = AllocateGadgetDataArray(output,"type",1,numBodies);
    else 
      this->Type = NULL;
  }
  else{
    // We are extending existing arrays
    vtkIdType currentSize = this->Positions->GetNumberOfPoints();
    vtkIdType newSize=currentSize+numBodies;
    // Arrays relevant for positions:
    this->Positions->SetNumberOfPoints(newSize);
    this->Vertices=vtkSmartPointer<vtkCellArray>::New();
    vtkIdType *cells = this->Vertices->WritePointer(newSize, newSize*2);
    for (vtkIdType i=0; i<newSize; ++i) {
      cells[i*2]   = 1;
      cells[i*2+1] = i;
    }
    output->SetVerts(this->Vertices);
    this->GlobalIds->Resize(newSize);
    vtkErrorMacro("new global ids array size is" << this->GlobalIds->GetSize() << "\n");

    if(this->Velocity!=NULL) {  
      this->Velocity->Resize(newSize);  
      //vtkErrorMacro("new velocity array size is" << this->Velocity->GetSize());

    }
    if(this->Mass!=NULL){
      this->Mass->Resize(newSize);
      //vtkErrorMacro("new mass array size is" << this->Mass->GetSize());

    }
    if(this->Energy!=NULL){
      this->Energy->Resize(newSize);
      //vtkErrorMacro("new energy array size is" << this->Energy->GetSize());

    }
    if(this->Type!=NULL){
      this->Type->Resize(newSize);
      //vtkErrorMacro("new type array size is" << this->Type->GetSize());

    }


    
  }
	
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
  

  this->PointDataArraySelection->AddArray("Mass");
  this->PointDataArraySelection->AddArray("Velocity");
  this->PointDataArraySelection->AddArray("Energy");
  this->PointDataArraySelection->AddArray("Type");

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
  char    gadget_file[MAXSTRING];
  int     no_gadget_files, i_gadget_file;

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
  
  if(this->FilePrefix != NULL && *this->FilePrefix != '\0') {
    vtkErrorMacro("supposedly more than one gadget file" << this->FilePrefix);
    // Supposedly more than one gadget file
    /* maybe there are multiple GADGET files ... count them! */
    no_gadget_files = 0;
    i_gadget_file   = 0;
    sprintf(gadget_file,"%s.%d",this->FilePrefix,i_gadget_file);
    while((icfile = fopen(gadget_file,"rb")) != NULL)
    {
      no_gadget_files++;
      i_gadget_file++;
      sprintf(gadget_file,"%s.%d",this->FilePrefix,i_gadget_file);
    }
    vtkErrorMacro("number of gadget files" << no_gadget_files);
    
    if(no_gadget_files > 1)
    {
      fprintf(stderr,"\nreading GADGET data from %d files:\n",no_gadget_files);
      gadget.no_gadget_files  = no_gadget_files;
      
      /* allocate temporary storage for no. of particles arrays */
      gadget.np[0]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
      gadget.np[1]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
      gadget.np[2]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
      gadget.np[3]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
      gadget.np[4]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
      gadget.np[5]     = (long *) calloc(gadget.no_gadget_files, sizeof(long *));
      
      /* read multi-GADGET files one by one */
      for(i_gadget_file=0; i_gadget_file<no_gadget_files; i_gadget_file++)
      {
        vtkErrorMacro("reading file" << i_gadget_file);
        sprintf(gadget_file,"%s.%d",this->FilePrefix,i_gadget_file);
        fprintf(stderr,"\n===================================================================\n");
        fprintf(stderr,"=> reading %s\n\n",gadget_file);
        icfile = fopen(gadget_file,"rb");
        
        /* tell read_gadget() which file we are using at the moment */
        gadget.i_gadget_file = i_gadget_file;
        
        /* read files... */
        this->ReadGadget(icfile,output);
        fclose(icfile);
      } 
      
      /* free temporary storage again */
      free(gadget.np[0]);
      free(gadget.np[1]);
      free(gadget.np[2]);
      free(gadget.np[3]);
      free(gadget.np[4]);
      free(gadget.np[5]);
    }
    else
    {
      /* there are no multi-GADGET files */
      fprintf(stderr,"\n\ninput: could not open file with IC's  %s\n",this->FilePrefix);
      exit(0);
    }

    
    
    
  }
  else {
    ////
    // Only one gadget file
    /////
    if((icfile = fopen(this->FileName,"rb")) != NULL) 
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
    else 
      {
        vtkErrorMacro("problem opening" << this->FileName);
      }
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
  vtkErrorMacro("read in gadget io header, omega0 is " << gadget.header.Omega0 << "\n");


  vtkErrorMacro("read in gadget io header NPs are " << gadget.header.np[0] << " " << gadget.header.np[1] <<  " " << gadget.header.np[2] << " " << gadget.header.np[3] << " " <<  gadget.header.np[4] << " " << gadget.header.np[5] << "\n");

  
  // here's where the ParaView specific code comes in
	
		// Allocate the arrays
  
  /* keep track of no. of particles in each GADGET file */
  gadget.np[0][gadget.i_gadget_file] = gadget.header.np[0];
  gadget.np[1][gadget.i_gadget_file] = gadget.header.np[1];
  gadget.np[2][gadget.i_gadget_file] = gadget.header.np[2];
  gadget.np[3][gadget.i_gadget_file] = gadget.header.np[3];
  gadget.np[4][gadget.i_gadget_file] = gadget.header.np[4];
  gadget.np[5][gadget.i_gadget_file] = gadget.header.np[5];
  
  /* conversion factors to Mpc/h, km/sec, Msun/h */
  x_fac  = this->LUnit;
  v_fac  = sqrt(gadget.header.expansion);
  m_fac  = this->MUnit;
  
  /* count total no. of particles in current file (and set massflag) */
  massflag    = 0;
  no_part     = 0;
  gadget.nall = 0;
  for(i=0;i<6;i++) 
  {
    no_part     += gadget.header.np[i];
    gadget.nall += gadget.header.nall[i];
    if(gadget.header.massarr[i] < MZERO && gadget.header.np[i] > 0)
      massflag=1;  
  }  
  
  /* be verbose */
  fprintf(stderr,"expansion factor: %lf\n",             gadget.header.expansion);
  fprintf(stderr,"redshift:         %lf\n",             gadget.header.redshift);
  fprintf(stderr,"boxsize:          %lf (%lf Mpc/h)\n", gadget.header.BoxSize,gadget.header.BoxSize*GADGET_LUNIT);
  fprintf(stderr,"omega0:           %lf\n",             gadget.header.Omega0);
  fprintf(stderr,"lambda0:          %lf\n",             gadget.header.OmegaLambda);
  fprintf(stderr,"HubbleParam:      %lf\n\n",           gadget.header.HubbleParam);
  
  fprintf(stderr,"gas:    np[0]=%9d\t nall[0]=%9d\t massarr[0]=%g\n",gadget.header.np[0],gadget.header.nall[0],gadget.header.massarr[0]); 
  fprintf(stderr,"halo:   np[1]=%9d\t nall[1]=%9d\t massarr[1]=%g\n",gadget.header.np[1],gadget.header.nall[1],gadget.header.massarr[1]); 
  fprintf(stderr,"disk:   np[2]=%9d\t nall[2]=%9d\t massarr[2]=%g\n",gadget.header.np[2],gadget.header.nall[2],gadget.header.massarr[2]); 
  fprintf(stderr,"bulge:  np[3]=%9d\t nall[3]=%9d\t massarr[3]=%g\n",gadget.header.np[3],gadget.header.nall[3],gadget.header.massarr[3]); 
  fprintf(stderr,"stars:  np[4]=%9d\t nall[4]=%9d\t massarr[4]=%g\n",gadget.header.np[4],gadget.header.nall[4],gadget.header.massarr[4]); 
  fprintf(stderr,"bndry:  np[5]=%9d\t nall[5]=%9d\t massarr[5]=%g\n",gadget.header.np[5],gadget.header.nall[5],gadget.header.massarr[5]); 
  
  fprintf(stderr,"\n-> reading %d particles from  GADGET file #%d/%d...\n\n", no_part, gadget.i_gadget_file+1, gadget.no_gadget_files);
  

  // TODO Size from the gadget header
	this->AllocateAllGadgetVariableArrays(no_part, output);
	// Loop through and add to PV arrays
	double pos[3];
	double vel[3];
	double den[3] = {0.0,0.0,0.0};
	double pmass;
	
  // Now reading in 
  if(this->Format == 1)
  {
    GADGET_SKIP;
    
    fread(DATA,sizeof(char),4,icfile);
    DATA[4] = '\0';
    GADGET_SKIP;
    
    GADGET_SKIP;
    fprintf(stderr,"reading %s",DATA);
  }
  else
  {
    fprintf(stderr,"reading ");
  }  
  GADGET_SKIP;

	
  
  // TODO: ignoring age and metals for now.
  for(int i=0;i<no_part;i++)
  {    
    /* read */
    if(DGADGET)
    {
      ReadDouble(icfile,&(ddummy[0]),this->Swap);
      ReadDouble(icfile,&(ddummy[1]),this->Swap);
      ReadDouble(icfile,&(ddummy[2]),this->Swap);
    }
    else
    {
      ReadFloat(icfile,&(fdummy[0]),this->Swap);
      ReadFloat(icfile,&(fdummy[1]),this->Swap);
      ReadFloat(icfile,&(fdummy[2]),this->Swap);
      ddummy[0] = fdummy[0];
      ddummy[1] = fdummy[1];
      ddummy[2] = fdummy[2];
    }
    
    /* get proper position in Part[] array */
    pid = this->GetPid(i);
    
    /* storage and conversion to comoving physical units */
    pos[0] = ddummy[0] * x_fac;
    pos[1] = ddummy[1] * x_fac;
    pos[2] = ddummy[2] * x_fac;  
		vel[0]=0;
		vel[1]=0;
		vel[2]=0;
		// all particle types have this
		this->Positions->SetPoint(pid, pos);
		if (this->Type)      this->Type->SetTuple1(pid, 0.0);
    
    if (this->Energy)      this->Energy->SetTuple1(pid, 0.0);
    
	}
  GADGET_SKIP;
  
  /*================= done read in GADGET particles =================*/
  
  // TODO: remove, just positions for now
  return 1;
  
  /*================= read in GADGET velocities =================*/
  if(this->Format == 1)
  {
    GADGET_SKIP;  
    
    fread(DATA,sizeof(char),4,icfile);
    DATA[4] = '\0';
    GADGET_SKIP;
    
    GADGET_SKIP;
    fprintf(stderr,"reading %s",DATA);
  }
  else
  {
    fprintf(stderr,"reading ");
  }
  
  GADGET_SKIP;
  fprintf(stderr,"(%8.2g MB) ... ",blklen/1024./1024.);
  
  for(i=0;i<no_part;i++)
  {
    /* read */
    if(DGADGET)
    {
      ReadDouble(icfile,&(ddummy[0]),this->Swap);
      ReadDouble(icfile,&(ddummy[1]),this->Swap);
      ReadDouble(icfile,&(ddummy[2]),this->Swap);
    }
    else
    {
      ReadFloat(icfile,&(fdummy[0]),this->Swap);
      ReadFloat(icfile,&(fdummy[1]),this->Swap);
      ReadFloat(icfile,&(fdummy[2]),this->Swap);
      ddummy[0] = fdummy[0];
      ddummy[1] = fdummy[1];
      ddummy[2] = fdummy[2];
    }
    
    /* get proper position in Part[] array */
    pid = this->GetPid(i);
    
    /* storage and conversion to comoving physical units */
    vel[0] = ddummy[0] * v_fac;
    vel[1] = ddummy[1] * v_fac;
    vel[2] = ddummy[2] * v_fac; 
    if (this->Velocity)  this->Velocity->SetTuple(pid, vel);

  }
  
  GADGET_SKIP;
  /*================= done read in GADGET velocities =================*/
  

  /*================= read in GADGET id's =================*/
  if(this->Format == 2)
  {
    GADGET_SKIP;
    
    fread(DATA,sizeof(char),4,icfile);
    DATA[4] = '\0';
    GADGET_SKIP;
    
    GADGET_SKIP;
    fprintf(stderr,"reading %s",DATA);
  }
  else
  {
    fprintf(stderr,"reading ");
  }
  
  GADGET_SKIP;
  
  for(i=0;i<no_part;i++)
  {
    /* get proper position in Part[] array */
    pid = this->GetPid(i);
    
    if(LGADGET)
    {
      ReadLong(icfile,&ldummy,this->Swap);
      this->GlobalIds->SetTuple1(pid,(vtkIdType)ldummy); 
    }
    else
    {
      ReadInt(icfile,&idummy,this->Swap);
      this->GlobalIds->SetTuple1(pid, (vtkIdType)ldummy);
    }
    
  }  
  
  GADGET_SKIP;
  /*================= done read in GADGET id's =================*/
  

  
  k = 0;
  /* massflag == 1 indicates that massarr[i] = 0 and hence need to read in particle masses */
  if(massflag==1) 
  {
    /*================= read in GADGET individual particle masses =================*/
    if(this->Format == 2)
    {
      GADGET_SKIP; 
      fread(DATA,sizeof(char),4,icfile);
      DATA[4] = '\0';
      GADGET_SKIP;
      
      GADGET_SKIP;
      fprintf(stderr,"reading %s",DATA);
    }
    else
    {
      fprintf(stderr,"reading ");
    }
    
    GADGET_SKIP;
    fprintf(stderr,"(%8.2g MB) ... ",blklen/1024./1024.);
    
    for(i=0;i<6;i++)
    {
      tot_mass[i] = 0.;
      if (gadget.header.np[i] > 0 && gadget.header.massarr[i] < MZERO  ) 
      {
        
        fprintf(stderr,"  %d    ",i);
        
        for(j=0; j<gadget.header.np[i]; j++)
        {
          /* read */
          if(DGADGET)
          {
            ReadDouble(icfile,&(ddummy[0]),this->Swap);
          }
          else
          {
            ReadFloat(icfile,&(fdummy[0]),this->Swap);
            ddummy[0] = fdummy[0];
          }
          
          /* get proper position in Part[] array */
          pid = this->GetPid(k);
          
          /* store */
          if (this->Mass)      this->Mass->SetTuple1(pid, ddummy[0]);

          tot_mass[i]    += ddummy[0];
          k++;
        }
      }
      else
      {
        /* simply copy appropriate massarr[i] to particles */
        for(j=0; j<gadget.header.np[i]; j++) 
        {
          /* get proper position in Part[] array */
          pid = this->GetPid(k);
          
          /* store */
          if (this->Mass)      this->Mass->SetTuple1(pid, gadget.header.massarr[i]);
          k++;
        }
        tot_mass[i] = gadget.header.np[i]*gadget.header.massarr[i];
      }
    }
    
    GADGET_SKIP;
    fprintf(stderr,"(%8.2g MB) done.\n",blklen/1024./1024.);
    /*================= read in GADGET individual particle masses =================*/
  } 
  
  /* simply copy appropriate massarr[i] to particles */
  else 
  {
    k=0;
    for(i=0;i<6;i++)
    {
      for(j=0;j<gadget.header.np[i];j++) 
      {
        /* get proper position in Part[] array */
        pid = this->GetPid(k);
        
        /* store */
        if (this->Mass)  this->Mass->SetTuple1(pid, gadget.header.massarr[i]);
        k++;
      }
      tot_mass[i] = gadget.header.np[i]*gadget.header.massarr[i];
    }
  }
  
  /* convert masses to Msun/h */
  k=0;
  for(i=0;i<6;i++)
  {
    for(j=0;j<gadget.header.np[i];j++) 
    {
      /* get proper position in Part[] array */
      pid = this->GetPid(k);
      if (this->Mass)      {
          double dmass = this->Mass->GetTuple1(pid);
          this->Mass->SetTuple1(pid,dmass*this->MUnit);
      }   
      k++;
    }
  }
  
  /*================= read in GADGET gas particle energies =================*/
  if(gadget.header.np[0] > 0) 
  {      
    if(this->Format == 2)
    {
      GADGET_SKIP;
      
      fread(DATA,sizeof(char),4,icfile);
      DATA[4] = '\0';
      GADGET_SKIP;
      
      GADGET_SKIP;
      fprintf(stderr,"reading %s",DATA);
    }
    else
    {
      fprintf(stderr,"reading ");
    }
    
    GADGET_SKIP; 
    fprintf(stderr,"(%8.2g MB) ... ",blklen/1024./1024.);
    
    for(i=0; i<gadget.header.np[0]; i++)
    {
      /* store */
      if(DGADGET)
      {
        ReadDouble(icfile,&(ddummy[0]),this->Swap);
      }
      else
      {
        ReadFloat(icfile,&(fdummy[0]),this->Swap);
        ddummy[0] = fdummy[0];
      }
      
      /* get proper position in Part[] array */
      pid = this->GetPid(i);
      /* store additional gas particle property */
      // TODO: add back in correct!
       //if (this->Energy) this->Energy->SetTuple1(pid,ddummy[0]);         
      if (this->Energy) this->Energy->SetTuple1(pid,(double)pid);         
    }
    
    GADGET_SKIP;
    fprintf(stderr,"(%8.2g MB) done.\n",blklen/1024./1024.);
  } 
  /*================= read in GADGET gas particle energies =================*/
  


  
  

	
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
  this->Mass        = NULL;
  this->Energy         = NULL;
	this->Type        = NULL;
  this->Velocity    = NULL;
  //

  return 0;
}




/*=============================================================================
 *                        get proper position in Part[] array
 *=============================================================================*/
long vtkGadgetReader::GetPid(int i)
{
  long pid;
  long itype, ifile;
  
  pid = 0;
  for(ifile=0; ifile<gadget.i_gadget_file; ifile++)
    for(itype=0; itype<6; itype++)
      pid += gadget.np[itype][ifile];
  pid += i;
  
  return(pid);
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
