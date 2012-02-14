/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkAddAdditionalAttribute.cxx,v $
=========================================================================*/
#include "vtkAddAdditionalAttribute.h"
#include "AstroVizHelpersLib/AstroVizHelpers.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStringArray.h"
#include "vtkSphereSource.h"
#include "vtkCenterOfMassFilter.h"
#include "vtkCellData.h"
#include "vtkPoints.h"
#include "vtkLine.h"
#include "vtkUnsignedCharArray.h"
#include "vtkSortDataArray.h"
#include "vtkMultiProcessController.h"
#include <vtkstd/string>
#include "vtkMath.h"
#include <stdio.h>
#include <stdlib.h>


vtkCxxRevisionMacro(vtkAddAdditionalAttribute, "$Revision: 1.72 $");
vtkStandardNewMacro(vtkAddAdditionalAttribute);

//----------------------------------------------------------------------------
vtkAddAdditionalAttribute::vtkAddAdditionalAttribute()
{
	this->SetInputArrayToProcess(
    0,
    0,
    0,
    vtkDataObject::FIELD_ASSOCIATION_POINTS_THEN_CELLS,
    vtkDataSetAttributes::SCALARS);
	this->AttributeFile = 0;
	this->AttributeName = 0; 
}

//----------------------------------------------------------------------------
vtkAddAdditionalAttribute::~vtkAddAdditionalAttribute()
{
  this->SetAttributeFile(0);
  this->SetAttributeName(0);
}

//----------------------------------------------------------------------------
void vtkAddAdditionalAttribute::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
	os << indent << "AttributeFile: "
	 << (this->AttributeFile ? this->AttributeFile : "(none)") << "\n";

}

//----------------------------------------------------------------------------
int vtkAddAdditionalAttribute::FillInputPortInformation(int, 
	vtkInformation* info)
{
  // This filter uses the vtkDataSet cell traversal methods so it
  // suppors any data set type as input.
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  return 1;
}


//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
int vtkAddAdditionalAttribute::RequestData(vtkInformation*,
	vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  // get input and output data
  vtkPointSet* input = vtkPointSet::GetData(inputVector[0]);

  vtkPointSet* output = vtkPointSet::GetData(outputVector);
	output->Initialize();
	output->ShallowCopy(input);
	// Make sure we are not running in parallel, this filter does not work in 
	// parallel
	// Gradually starting to make this work in parallel; removing this for now
	/*
	if(RunInParallel(vtkMultiProcessController::GetGlobalController()))
		{
		vtkErrorMacro("This filter is not supported in parallel.");
		return 0;
		}
	*/
	// Make sure we have a file to read.
  if(!this->AttributeFile)
	  {
    vtkErrorMacro("An attribute file must be specified.");
    return 0;
    }


  if(strcmp(this->AttributeName,"")==0)	
  {
    vtkErrorMacro("Please specify an attribute name.");
    return 0;
  }
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  // get this->UpdatePiece information
  this->UpdatePiece = \
      outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  this->UpdateNumPieces =			\
    outInfo->Get(
		 vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());


  // this algorithm only works if we first sort the globalIdArray
  // in increasing order of ids. then only the first
  // call to SeekInAsciiAttribute has the possibility to involve a long seek.
  if(this->AttributeFileFormatType==FORMAT_SKID_ASCII)
  {

     // Get name of data array containing global id
     vtkDataArray* globalIdArray = this->GetInputArrayToProcess(0, inputVector);
    if(!globalIdArray)
     {
       vtkErrorMacro("Failed to locate global id array");
       return 0;
     }
    // open file
    ifstream attributeInFile(this->AttributeFile);
    if(strcmp(this->AttributeFile,"")==0||!attributeInFile)
 		{
      vtkErrorMacro("Error opening attribute file: " << this->AttributeFile);
      return 0;
     }

    vtkSortDataArray::Sort(globalIdArray);
    if(globalIdArray->GetNumberOfTuples() == \
      output->GetPoints()->GetNumberOfPoints())
      {
      // read additional attribute for all particles
      AllocateDataArray(output,this->AttributeName,1,
        output->GetPoints()->GetNumberOfPoints());		
      double attributeData;
      //always skip the header, which is the total number of bodies
      attributeInFile >> attributeData;
      unsigned long formerGlobalId = -1;
      for(unsigned long localId=0; 
        localId < globalIdArray->GetNumberOfTuples(); 
        localId++)
        {
        vtkIdType globalId = globalIdArray->GetComponent(localId,0);
        // seeking to next data id
        attributeData = SeekInAsciiAttributeFile(attributeInFile,
          globalId-formerGlobalId);
        formerGlobalId=globalId;
        // place attribute data in output
        SetDataValue(output,this->AttributeName,localId,
          &attributeData);			
        }
      }
    // closing file
    attributeInFile.close();

  }
  else if(this->AttributeFileFormatType==FORMAT_HOP_DENSITY_BIN)
  {
     // This part not yet supported in Parallel!
    // TODO: check file opens correctly
     FILE *infile = fopen(this->AttributeFile, "r");
     int numberParticles[1];
     int error = fread(numberParticles, sizeof(int),1, infile);
     float attributeData[1];
    // read additional attribute for all particles
     vtkDebugMacro("number hop particles: "<< numberParticles[0]);

     unsigned long pieceSize = floor(numberParticles[0]*1.0/this->UpdateNumPieces);
     unsigned long beginIndex = this->UpdatePiece*pieceSize;
     unsigned long endIndex = (this->UpdatePiece == this->UpdateNumPieces - 1) ?\
       numberParticles[0] : (this->UpdatePiece+1)*pieceSize;
     unsigned long numberElts = endIndex-beginIndex;
     vtkDebugMacro("updatepiece: " << this->UpdatePiece  
		   << " pieceSize: " << pieceSize
		   << " numberElts: " << numberElts
		   << " outputNumPoints: " << output->GetPoints()->GetNumberOfPoints()
		   << " beginIndex: " << beginIndex
		   << " endIndex: " << endIndex << "\n");



    if(numberParticles[0] == \
       output->GetPoints()->GetNumberOfPoints()) 
      {
	AllocateDataArray(output,this->AttributeName,1,
			  output->GetPoints()->GetNumberOfPoints());		

	fseek (infile , sizeof(float)*beginIndex , SEEK_SET);
	for(unsigned long idx=0; idx<numberElts; idx++)
	  {
	    error = fread(attributeData, sizeof(float),1, infile);
	    SetDataValue(output,this->AttributeName,idx,
			 &attributeData[0]);			
	    
	  }
	return 1;
      }
    
    
    else {
      vtkErrorMacro("number of points in input must be equal to number of points in HOP file " << numberParticles[0]);
    }
    
  }
  else if(this->AttributeFileFormatType==FORMAT_HOP_MARKFILE_BIN)
  {
    /// TODO TODO: TEST THE PARALLEL IMPLEMENTATION
    // TODO: check file opens correctly
     FILE *infile = fopen(this->AttributeFile, "rb");
     int numberParticles[1];
     int error = fread(numberParticles, sizeof(int),1, infile);
     int attributeData[1];
     
     float attributeDataF[1];
     vtkDebugMacro("number hop particles: "<< numberParticles[0]);

     unsigned long pieceSize = floor(numberParticles[0]*1.0/this->UpdateNumPieces);
     unsigned long beginIndex = this->UpdatePiece*pieceSize;
     unsigned long endIndex = (this->UpdatePiece == this->UpdateNumPieces - 1) ?\
       numberParticles[0] : (this->UpdatePiece+1)*pieceSize;
     unsigned long numberElts = endIndex-beginIndex;

     // here's where we do an fseek
    // read additional attribute for all particles
     vtkDebugMacro("updatepiece: " << this->UpdatePiece  
		   << " pieceSize: " << pieceSize
		   << " numberElts: " << numberElts
		   << " outputNumPoints: " << output->GetPoints()->GetNumberOfPoints()
		   << " beginIndex: " << beginIndex
		   << " endIndex: " << endIndex << "\n");

     if(numberElts ==					\
	output->GetPoints()->GetNumberOfPoints()) 
       {

	AllocateDataArray(output,this->AttributeName,1,numberElts);		
	fseek(infile, sizeof(int)*beginIndex , SEEK_SET);
	for(unsigned long idx=0; idx<numberElts; idx++)
	  {
	    error = fread(attributeData, sizeof(int),1, infile);
	    attributeDataF[0]=(float)attributeData[0];
	    SetDataValue(output,this->AttributeName,idx,
			 &attributeDataF[0]);			
	    
	  }
	return 1;
       }
    

     else {
      vtkErrorMacro("number of points in input must be equal to number of points in HOP file " << numberParticles[0]);
      }
    
  }

  else {
    return 0;
  }
}
