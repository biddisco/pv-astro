/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkGadgetReader.h,v $

  Copyright (c) Christine Corbett Moran
  All rights reserved.
     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkGadgetReader - Read points from a Gadget standard binary file
// .SECTION Description
// Read points from a Gadget standard binary file. Fully parallel. Has ability
// to read in additional attributes from an ascii file, and to only load in
// marked particles but both these functions are serial only.
#ifndef __vtkGadgetReader_h
#define __vtkGadgetReader_h

#include "vtkPolyDataAlgorithm.h" // superclass

#include "vtkSmartPointer.h"
#include "tipsylib/ftipsy.hpp" // functions take Gadget particle objects
#include <vtkstd/vector>




class vtkPolyData;
class vtkCharArray;
class vtkIdTypeArray;
class vtkDoubleArray;
class vtkPoints;
class vtkCellArray;
class vtkDataArraySelection;
class   vtkMultiProcessController;
class VTK_EXPORT vtkGadgetReader : public vtkPolyDataAlgorithm
{
public:
  static vtkGadgetReader* New();
  vtkTypeMacro(vtkGadgetReader,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  // Description:
  // Set/Get the name of the file from which to read points.
	vtkSetStringMacro(FileName);
 	vtkGetStringMacro(FileName);

  
  vtkSetStringMacro(FilePrefix);
 	vtkGetStringMacro(FilePrefix);

	
	// Description:
  // Set/Get the LUnit
	vtkSetMacro(LUnit,double);
 	vtkGetMacro(LUnit,double);

	
  // Description:
  // Set/Get the MUnit
	vtkSetMacro(MUnit,double);
 	vtkGetMacro(MUnit,double);


	// Description:
  // Set/Get the Format
	vtkSetMacro(Format,bool);
 	vtkGetMacro(Format,bool);

  // Description:
  // Set/Get the Swap
	vtkSetMacro(Swap,bool);
 	vtkGetMacro(Swap,bool);
  
  // Description:
  // Set/Get the LGadget
	vtkSetMacro(LGadget,bool);
 	vtkGetMacro(LGadget,bool);

  
  // Description:
  // Set/Get the LGadget
	vtkSetMacro(DGadget,bool);
 	vtkGetMacro(DGadget,bool);

	// Description:
  // An H5Part file may contain multiple arrays
  // a GUI (eg Paraview) can provide a mechanism for selecting which data arrays
  // are to be read from the file. The PointArray variables and members can
  // be used to query the names and number of arrays available
  // and set the status (on/off) for each array, thereby controlling which
  // should be read from the file. Paraview queries these point arrays after
  // the (update) information part of the pipeline has been updated, and before the
  // (update) data part is updated.
  int         GetNumberOfPointArrays();
  const char* GetPointArrayName(int index);
  int         GetPointArrayStatus(const char* name);
  void        SetPointArrayStatus(const char* name, int status);
  void        DisableAll();
  void        EnableAll();
  void        Disable(const char* name);
  void        Enable(const char* name);
  //
  int         GetNumberOfPointArrayStatusArrays() { return GetNumberOfPointArrays(); }
  const char* GetPointArrayStatusArrayName(int index) { return GetPointArrayName(index); }
  int         GetPointArrayStatusArrayStatus(const char* name) { return GetPointArrayStatus(name); }
  void        SetPointArrayStatusArrayStatus(const char* name, int status) { SetPointArrayStatus(name, status); }

// The BTX, ETX comments bracket the portion of the code which should not be
// attempted to wrap for use by python, specifically the code which uses
// C++ templates as this code is unable to be wrapped. DO NOT REMOVE. 
//BTX
protected:
  vtkGadgetReader();
  ~vtkGadgetReader();
	char* FileName;
  char* FilePrefix;

	double LUnit;
  double MUnit;
	bool Format;
  bool Swap;
  bool LGadget;
  bool DGadget;
	int RequestInformation(vtkInformation*,	vtkInformationVector**,
		vtkInformationVector*);

  int RequestData(vtkInformation*,vtkInformationVector**,
    vtkInformationVector*);

  vtkIdType                       ParticleIndex;
  vtkSmartPointer<vtkIdTypeArray> GlobalIds;
  vtkSmartPointer<vtkPoints>      Positions;
  vtkSmartPointer<vtkCellArray>   Vertices;

  vtkSmartPointer<vtkDoubleArray>   Mass;
	vtkSmartPointer<vtkDoubleArray>		 Type;
	vtkSmartPointer<vtkDoubleArray>		 Energy;
  vtkSmartPointer<vtkDoubleArray>   Velocity;

	
  //
  vtkMultiProcessController *Controller;
  int           UpdatePiece;
  int           UpdateNumPieces;

  // To allow paraview gui to enable/disable scalar reading
  vtkDataArraySelection* PointDataArraySelection;

private:
  vtkGadgetReader(const vtkGadgetReader&);  // Not implemented.
  void operator=(const vtkGadgetReader&);  // Not implemented.
	/* Helper functions for storing data in output vector*/
	// Description:
	// allocates all vtk arrays for Tipsy variables and places them 
	// in the output vector
	void AllocateAllGadgetVariableArrays(vtkIdType numBodies,
																			vtkPolyData* output);
  int ReadGadgetAllFiles();
  int ReadGadget(FILE *icfile);
  long GetPid(int i);
//ETX

};
#endif
