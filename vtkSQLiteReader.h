/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkSQLiteReader.h,v $

  Copyright (c) Rafael Kueng
  All rights reserved.
     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkSQLiteReader - Read points from a SQLite database
// .SECTION Description
// Reads points from a SQLite DB and displays them on screen


#ifndef __vtkSQLiteReader_h
#define __vtkSQLiteReader_h

#include "vtkPolyDataAlgorithm.h" // superclass
#include "sqlitelib/sqlite3.h" // sqlite headerfile

#include "vtkObjectFactory.h"
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include <vtkstd/vector>

// probably there better ways than list...
//TODO check this...
//#include <list>
//using namespace std;

class vtkPolyData;
class vtkCharArray;
class vtkIdTypeArray;
class vtkFloatArray;
class vtkPoints;
class vtkCellArray;
class vtkDataArraySelection;

//BTX
//struct halo {
//	int id;
//	vtkPoints coordinates;
//};
//ETX

class VTK_EXPORT vtkSQLiteReader : public vtkPolyDataAlgorithm
{
public:
	static vtkSQLiteReader* New();
	vtkTypeRevisionMacro(vtkSQLiteReader,vtkPolyDataAlgorithm);
	// Set/Get the name of the file from which to read points.
	vtkSetStringMacro(FileName);
 	vtkGetStringMacro(FileName);

	vtkSetMacro(DisplaySnapshot,int);
	vtkGetMacro(DisplaySnapshot,int);

//BTX
protected:
	vtkSQLiteReader();
	~vtkSQLiteReader();
	char* FileName;

	int RequestInformation(vtkInformation*,	vtkInformationVector**,
		vtkInformationVector*);

	int RequestData(vtkInformation*,vtkInformationVector**,
		vtkInformationVector*);

	

// functions

// structs
	struct snapinfo {
		int snapshotnr;
		double redshift;
		double time;
	} snapinfo;

// variables
	int numSnaps;
	int DisplaySnapshot;
	std::vector<vtkSmartPointer<vtkPolyData>> data;

	/* not used
	// for halos
	vtkIdType							ParticleIndex;
	vtkSmartPointer<vtkIdTypeArray>		ParticleId;
	vtkSmartPointer<vtkPoints>			Position;
	vtkSmartPointer<vtkCellArray>		Cells;
	vtkSmartPointer<vtkFloatArray>		Velocity;
	vtkSmartPointer<vtkFloatArray>		nParticles;

	vtkSmartPointer<vtkFloatArray>		mVir;
	vtkSmartPointer<vtkFloatArray>		rVir;
	vtkSmartPointer<vtkFloatArray>		RHO;

	// for snapshots
	vtkSmartPointer<vtkIdTypeArray>		SnapId;
	vtkSmartPointer<vtkFloatArray>		Redshift;
	vtkSmartPointer<vtkFloatArray>		Time;
	vtkSmartPointer<vtkDataArray>		Halos;

	// for tracks
	vtkSmartPointer<vtkIdTypeArray>		TrackId;
	vtkSmartPointer<vtkPolyData>		Lines;
*/


// constants


private:
	vtkSQLiteReader(const vtkSQLiteReader&);  // Not implemented.
	void operator=(const vtkSQLiteReader&);  // Not implemented.

	//functions
	int openDB(char*);
	int vtkSQLiteReader::readSnapshots(
		std::vector<vtkSmartPointer<vtkPolyData>> *);
	int vtkSQLiteReader::ReadHeader(vtkInformationVector*);
	int RequestDataDemo(vtkInformationVector*);

	// helpers
	vtkStdString vtkSQLiteReader::Int2Str(int);
	int vtkSQLiteReader::SQLQuery(vtkStdString, sqlite3_stmt*);

	//variables
	sqlite3 * db;
	bool dataIsRead;

	//constants


//ETX
};

#endif