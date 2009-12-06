/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkTipsyReader.cxx,v $
=========================================================================*/
#include "vtkTipsyReader.h"
#include "AstroVizHelpersLib/AstroVizHelpers.h"
#include "vtkUnstructuredGrid.h"
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
#include <cmath>
#include <assert.h>

vtkCxxRevisionMacro(vtkTipsyReader, "$Revision: 1.0 $");
vtkStandardNewMacro(vtkTipsyReader);

//----------------------------------------------------------------------------
vtkTipsyReader::vtkTipsyReader()
{
  this->MarkFileName = 0; // this file is optional
  this->FileName = 0;
	this->ReadPositionsOnly = 0;
	this->DistributeDataOn = 1;
  this->SetNumberOfInputPorts(0); 
}

//----------------------------------------------------------------------------
vtkTipsyReader::~vtkTipsyReader()
{
  this->SetFileName(0);
  this->SetMarkFileName(0);
	this->SetReadPositionsOnly(0);
	this->SetDistributeDataOn(0);
}

//----------------------------------------------------------------------------
void vtkTipsyReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "FileName: "
     << (this->FileName ? this->FileName : "(none)") << "\n"
		 << indent << "MarkFileName: "
		 << (this->MarkFileName ? this->MarkFileName : "(none)") << "\n";
}

//----------------------------------------------------------------------------
TipsyHeader vtkTipsyReader::ReadTipsyHeader(ifTipsy& tipsyInfile)
{
	// Initializing for reading  
	TipsyHeader       h; 
	// Reading in the header
  tipsyInfile >> h;
	return h;
}

//----------------------------------------------------------------------------
vtkstd::vector<unsigned long> vtkTipsyReader::ReadMarkedParticleIndices(
	TipsyHeader& tipsyHeader,ifTipsy& tipsyInfile)
{
	ifstream markInFile(this->MarkFileName);
	vtkstd::vector<unsigned long> markedParticleIndices;
	if(!markInFile)
 		{
 		vtkErrorMacro("Error opening marked particle file: " 
									<< this->MarkFileName 
									<< " please specify a valid mark file or none at all.\
									 For now reading all particles.");
 		}
	else
		{
		unsigned long mfIndex,mfBodies,mfGas,mfStar,mfDark,numBodies;
		// first line of the mark file is of a different format:
		// intNumBodies intNumGas intNumStars
		if(markInFile >> mfBodies >> mfGas >> mfStar)
			{
	 		mfDark=mfBodies-mfGas-mfStar;
			if(mfBodies!=tipsyHeader.h_nBodies || mfDark!=tipsyHeader.h_nDark \
				|| mfGas!=tipsyHeader.h_nSph || mfStar!=tipsyHeader.h_nStar)
	 			{
	 			vtkErrorMacro("Error opening marked particle file, wrong format,\
	 										number of particles do not match Tipsy file: " 
											<< this->MarkFileName 
											<< " please specify a valid mark file or none at all.\
											For now reading all particles.");
	 			}
			else
				{
				// The rest of the file is is format: intMarkIndex\n
				// Read in the rest of the file file, noting the marked particles
				while(markInFile >> mfIndex)
					{
					// Insert the next marked particle index into the vector
					// subtracting one as the indices in marked particle file
					// begin at 1, not 0
					markedParticleIndices.push_back(mfIndex-1);
					}
				// closing file
				markInFile.close();
				// read file successfully
				vtkDebugMacro("Read " << numBodies<< " marked point indices.");
				}	
	 		}
		}
		//empty if none were read, otherwise size gives number of marked particles
		return markedParticleIndices;
}

//----------------------------------------------------------------------------
void vtkTipsyReader::ReadAllParticles(TipsyHeader& tipsyHeader,
	ifTipsy& tipsyInfile,int piece,int numPieces,vtkPolyData* output)
{
	unsigned long pieceSize = floor(tipsyHeader.h_nBodies*1./numPieces);
	unsigned long beginIndex = piece*pieceSize;
	unsigned long endIndex = (piece == numPieces - 1) ? \
	 	tipsyHeader.h_nBodies : (piece+1)*pieceSize;
	// Allocates vtk scalars and vector arrays to hold particle data, 
	this->AllocateAllTipsyVariableArrays(endIndex-beginIndex,output);
	for(unsigned long i=beginIndex; i < endIndex; i++)  
  	{ 
		this->ReadParticle(i,tipsyHeader,tipsyInfile,output);
  	}
}

//----------------------------------------------------------------------------
void vtkTipsyReader::ReadMarkedParticles(
	vtkstd::vector<unsigned long>& markedParticleIndices,
	TipsyHeader& tipsyHeader,
	ifTipsy& tipsyInfile,
	vtkPolyData* output)
{
	// Allocates vtk scalars and vector arrays to hold particle data, 
	// As marked file was read, only allocates numBodies which 
	// now equals the number of marked particles
	this->AllocateAllTipsyVariableArrays(markedParticleIndices.size(),output);
	unsigned long nextMarkedParticleIndex=0;
	for(vtkstd::vector<unsigned long>::iterator it = \
	 	markedParticleIndices.begin(); it != markedParticleIndices.end(); ++it)		
		{
 		nextMarkedParticleIndex=*it;
		// reading in the particle
		this->ReadParticle(nextMarkedParticleIndex,tipsyHeader,
			tipsyInfile,output);
		}
}

//----------------------------------------------------------------------------
tipsypos::section_type vtkTipsyReader::SeekToIndex(unsigned long index,
	TipsyHeader& tipsyHeader, ifTipsy& tipsyInfile)
{
	if(index < tipsyHeader.h_nSph)
		{
		// we are seeking a gas particle
		tipsyInfile.seekg(tipsypos(tipsypos::gas,index));	
		return tipsypos::gas;
		}
	else if(index < tipsyHeader.h_nDark+tipsyHeader.h_nSph)
		{
		// we are seeking a dark particle, which begin at index 0 after
		// gas particles
		index-=tipsyHeader.h_nSph;
		tipsyInfile.seekg(tipsypos(tipsypos::dark,index));
		return tipsypos::dark;
		}
	else if(index < tipsyHeader.h_nBodies)
		{
		// we are seeking a star particle, which begin at index zero after gas
		// and dark particles
		index-=(tipsyHeader.h_nSph+tipsyHeader.h_nDark);
		tipsyInfile.seekg(tipsypos(tipsypos::star,index));	
		return tipsypos::star;
		}
	else
		{
		vtkErrorMacro("An index is greater than the number of particles in the 	file, unable to read");
		return tipsypos::invalid;
		}
}
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
//i,tipsyHeader,tipsyInfile,output)
vtkIdType vtkTipsyReader::ReadParticle(unsigned long index, 
	TipsyHeader& tipsyHeader, ifTipsy& tipsyInfile, vtkPolyData* output) 
{
  // allocating variables for reading
  vtkIdType id;
  TipsyGasParticle  g;
  TipsyDarkParticle d;
  TipsyStarParticle s;
	tipsypos::section_type particleSection=this->SeekToIndex(index,
		tipsyHeader,tipsyInfile);
	switch(particleSection)
		{
   case tipsypos::gas:
     tipsyInfile >> g;
		 id=this->ReadGasParticle(output,g);
     break;
   case tipsypos::dark:
     tipsyInfile >> d;
		 id=this->ReadDarkParticle(output,d);
     break;
   case tipsypos::star:
     tipsyInfile >> s;
		 id=this->ReadStarParticle(output,s);
     break;
   default:
     assert(0);
     break;
    }
	SetIdTypeValue(output,"global id",id,index);
  return id;
}


//----------------------------------------------------------------------------
vtkIdType vtkTipsyReader::ReadBaseParticle(vtkPolyData* output,
	TipsyBaseParticle& b) 
{
	vtkIdType id = SetPointValue(output,b.pos);
  if(!this->ReadPositionsOnly)
		{
		SetDataValue(output,"velocity",id,b.vel);
	  SetDataValue(output,"mass",id,&b.mass);
		SetDataValue(output,"potential",id,&b.phi);
		}
	return id;
}

//----------------------------------------------------------------------------
vtkIdType vtkTipsyReader::ReadGasParticle(vtkPolyData* output,
 	TipsyGasParticle& g)
{
 	vtkIdType id=ReadBaseParticle(output,g);
  if(!this->ReadPositionsOnly)
		{
	  SetDataValue(output,"rho",id,&g.rho);
	  SetDataValue(output,"temperature",id,&g.temp);
	  SetDataValue(output,"hsmooth",id,&g.hsmooth);
	  SetDataValue(output,"metals",id,&g.metals);	
		}
	return id;
}

//----------------------------------------------------------------------------
vtkIdType vtkTipsyReader::ReadStarParticle(vtkPolyData* output,
 	TipsyStarParticle& s)
{
	vtkIdType id=ReadBaseParticle(output,s);
	if(!this->ReadPositionsOnly)
		{
	  SetDataValue(output,"eps",id,&s.eps);
	  SetDataValue(output,"metals",id,&s.metals);
	  SetDataValue(output,"tform",id,&s.tform);
		}
	return id;
}
//----------------------------------------------------------------------------	
vtkIdType vtkTipsyReader::ReadDarkParticle(vtkPolyData* output,
	TipsyDarkParticle& d)
{
 	vtkIdType id=this->ReadBaseParticle(output,d);
  if(!this->ReadPositionsOnly)
		{
	 	SetDataValue(output,"eps",id,&d.eps);
		}
	return id;
}
		
//----------------------------------------------------------------------------
void vtkTipsyReader::AllocateAllTipsyVariableArrays(unsigned long numBodies,
	vtkPolyData* output)
{
  // Allocate objects to hold points and vertex cells. 
	// Storing the points and cells in the output data object.
  output->SetPoints(vtkSmartPointer<vtkPoints>::New());
  output->SetVerts(vtkSmartPointer<vtkCellArray>::New()); 
  AllocateIdTypeDataArray(output,"global id",1,numBodies);
	// Only allocate the other arrays if we are to read this data in
	// for global ids
  if(!this->ReadPositionsOnly)
		{
		// the default scalars to be displayed
		// Only allocate if we are to read thes in
	  AllocateDataArray(output,"potential",1,numBodies);
		// the rest of the scalars
	  AllocateDataArray(output,"mass",1,numBodies);
	  AllocateDataArray(output,"eps",1,numBodies);
	  AllocateDataArray(output,"rho",1,numBodies);
	  AllocateDataArray(output,"hsmooth",1,numBodies);
	  AllocateDataArray(output,"temperature",1,numBodies);
	  AllocateDataArray(output,"metals",1,numBodies);
	  AllocateDataArray(output,"tform",1,numBodies);
		// the default vectors to be displayed
	  AllocateDataArray(output,"velocity",3,numBodies);
		}
}
//----------------------------------------------------------------------------
int vtkTipsyReader::RequestInformation(
	vtkInformation* vtkNotUsed(request),
	vtkInformationVector** vtkNotUsed(inputVector),
	vtkInformationVector* outputVector)
{
	vtkInformation* outInfo = outputVector->GetInformationObject(0);
	// means that the data set can be divided into an arbitrary number of pieces
	outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(),
		-1);
	return 1;
}
/*
* Reads a file, optionally only the marked particles from the file, 
* in the following order:
* 1. Open Tipsy binary
* 2. Read Tipsy header (tells us the number of particles of each type we are 
*    dealing with)
* NOTE: steps 3, 5 are currently not parallel
* 3. Read mark file indices from marked particle file, if there is one
* 4. Read either marked particles only or all particles
* 5. If an attribute file is additionally specified, reads this additional
* 	 attribute into a data array, reading only those marked if necessary.
*/
//----------------------------------------------------------------------------
int vtkTipsyReader::RequestData(vtkInformation*,
	vtkInformationVector**,vtkInformationVector* outputVector)
{
	// Make sure we have a file to read.
  if(!this->FileName)
	  {
    vtkErrorMacro("A FileName must be specified.");
    return 0;
    }
	// Open the tipsy standard file and abort if there is an error.
	ifTipsy tipsyInfile;
  tipsyInfile.open(this->FileName,"standard");
  if (!tipsyInfile.is_open()) 
		{
	  vtkErrorMacro("Error opening file " << this->FileName);
	  return 0;	
    }
	//All helper functions will need access to this
	vtkSmartPointer<vtkPolyData> tipsyReadInitialOutput = \
	 	vtkSmartPointer<vtkPolyData>::New();
	tipsyReadInitialOutput->Initialize();
	vtkUnstructuredGrid* output = vtkUnstructuredGrid::GetData(
		outputVector);
	// This tells us which portion of the file to read, relevant if in parallel
	vtkInformation* outInfo = outputVector->GetInformationObject(0);
	int piece = outInfo->Get(
		vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
	int numPieces =outInfo->Get(
		vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
  // Read the header from the input
	TipsyHeader tipsyHeader=this->ReadTipsyHeader(tipsyInfile);
	// Next considering whether to read in a mark file, 
	// and if so whether that reading was a success 
	vtkstd::vector<unsigned long> markedParticleIndices;
	if(strcmp(this->MarkFileName,"")!=0)
		{
		// Reading only marked particles
		// Make sure we are not running in parallel, this filter does not work in 
		// parallel
		if(RunInParallel(vtkMultiProcessController::GetGlobalController()))
			{
			vtkErrorMacro("Reading from a mark file is not supported in parallel.");
			return 0;
			}
		vtkDebugMacro("Reading marked point indices from file:" 
			<< this->MarkFileName);
		markedParticleIndices=this->ReadMarkedParticleIndices(tipsyHeader,
			tipsyInfile);
		}
  // Read every particle and add their position to be displayed, 
	// as well as relevant scalars
	if(markedParticleIndices.empty())
		{
		// no marked particle file or there was an error reading the mark file, 
		// so reading all particles
		vtkDebugMacro("Reading all points from file " << this->FileName);
		this->ReadAllParticles(tipsyHeader,tipsyInfile,piece,numPieces,
			tipsyReadInitialOutput);
		}
	else 
		{
		//reading only marked particles
		assert(numPieces==1);
		vtkDebugMacro("Reading only the marked points in file: " \
				<< this->MarkFileName << " from file " << this->FileName);
		this->ReadMarkedParticles(markedParticleIndices,
			tipsyHeader,tipsyInfile,
			tipsyReadInitialOutput);	
		}
  // Close the tipsy in file.
	tipsyInfile.close();
	// If we need to, run D3 on the tipsyReadInitialOutput
	// producing one level of ghost cells
	if(this->GetDistributeDataOn() && \
	 	RunInParallel(vtkMultiProcessController::GetGlobalController()))
		{
		vtkSmartPointer<vtkDistributedDataFilter> d3 = \
		    vtkSmartPointer<vtkDistributedDataFilter>::New();
		d3->AddInput(tipsyReadInitialOutput);
		d3->UpdateInformation();
		// TODO: work with ghost levels again
		/*
		// adds one ghostlevel
		vtkStreamingDemandDrivenPipeline* exec = \
		 	static_cast<vtkStreamingDemandDrivenPipeline*>(d3->GetExecutive()); 	
		exec->SetUpdateExtent(exec->GetOutputInformation(0), piece, numPieces, 1); 
		*/
		d3->Update();
		// Changing output to output of d3
	 	output->ShallowCopy(d3->GetOutput()); 
		// TODO: work with ghost levels again
		/* 
		output->GetInformation()->Set(
		 	vtkDataObject::DATA_NUMBER_OF_GHOST_LEVELS(), 1);
		*/
		}
	else
		{
		output->ShallowCopy(tipsyReadInitialOutput);
		}
	// Read Successfully
	vtkDebugMacro("Read " << output->GetPoints()->GetNumberOfPoints() \
		<< " points.");
 	return 1;
}
