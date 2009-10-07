/*=========================================================================
Modified from vtkSimplePointsReader and from Doug Potter's Tipsylib
 -christine
=========================================================================*/
// .NAME vtkTipsyReader - Read points from a Tipsy standard binary file
// .SECTION Description
// here is the desciprtion
#ifndef __vtkTipsyReader_h
#define __vtkTipsyReader_h
#include "vtkPolyDataAlgorithm.h" // needed as this class extends vtkPolyDataAlgorithm
#include "tipsylib/ftipsy.hpp" // needed for functions which take Tipsy particles as arguments
#include "vtkSmartPointer.h" // needed for the functions to initialize arrays
#include "vtkPolyData.h" // needed as most helper functions modify output which is vtkPolyData
#include <queue> // needed for FIFO queue used to store marked particles
using std::queue;
class VTK_IO_EXPORT vtkTipsyReader : public vtkPolyDataAlgorithm
{
public:
  static vtkTipsyReader* New();
  vtkTypeRevisionMacro(vtkTipsyReader,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  // Description:
  // Set/Get the name of the file from which to read points.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);
  // Description:
  // Set/Get the name of the file from which to read the marked points.
  vtkSetStringMacro(MarkFileName);
  vtkGetStringMacro(MarkFileName);
  // Description:
  // Set/Get the name of the file from which to get additional attributes
  vtkSetStringMacro(AttributeFileName);
  vtkGetStringMacro(AttributeFileName);

protected:
  vtkTipsyReader();
  ~vtkTipsyReader();
  char* FileName;
	char* MarkFileName;
	char* AttributeFileName;
  int RequestData(vtkInformation*,
                  vtkInformationVector**,
                  vtkInformationVector*);
private:
  vtkTipsyReader(const vtkTipsyReader&);  // Not implemented.
  void operator=(const vtkTipsyReader&);  // Not implemented.
	/* Help functions for reading */
	// Description:
	// Reads the tipsy header. 
	TipsyHeader ReadTipsyHeader(ifTipsy& tipsyInfile);
	// Description:
	// Reads all particles from the tipsy file
	void ReadAllParticles(TipsyHeader& tipsyHeader,ifTipsy& tipsyInfile,vtkPolyData* output);
	// Description:
	// reads in a particle (either gas, dark or star as appropriate) from the tipsy in file of this class
	vtkIdType ReadParticle(ifTipsy& tipsyInFile,vtkPolyData* output);
	// Description:
	// reads in a particle (either gas, dark or star as appropriate) from 
	// the tipsy in file of this class, also reads in files from an
	// attribute array specified by the user
	vtkIdType ReadParticle(ifstream& attributeInFile,ifTipsy& tipsyInfile,vtkPolyData* output);
	// Description:
	// reads variables common to all particles
	vtkIdType ReadBaseParticle(vtkPolyData* output, TipsyBaseParticle& b);
	// The BTX, ETX comments bracket the portion of the code which should not be
	// attempted to wrap for use by python, specifically the code which uses
	// C++ templates as this code is unable to be wrapped. DO NOT REMOVE.
	//BTX
	// Description:
	// Reads only Marked particles from the tipsy file. Must be called after function ReadMarkedParticleIndices.
	void ReadMarkedParticles(queue<int> markedParticleIndices,TipsyHeader& tipsyHeader,ifTipsy& tipsyInfile,vtkPolyData* output);
	// Description:
	// reads in an array of the indices of marked particles from a file, returns a queue of marked particles
	// which is empty if reading was unsucessful.
	queue<int> ReadMarkedParticleIndices(TipsyHeader& tipsyHeader,ifTipsy& tipsyInfile);
	//ETX
	/* Helper functions for storing data in output vector*/
	// Description:
	// allocates all vtk arrays for Tipsy variables and places them in the output vector
	void AllocateAllTipsyVariableArrays(TipsyHeader& tipsyHeader,vtkPolyData* output);
};
#endif
