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

class vtkPolyData;
class vtkCharArray;
class vtkIdTypeArray;
class vtkFloatArray;
class vtkPoints;
class vtkCellArray;
class vtkDataArraySelection;

class VTK_EXPORT vtkSQLiteReader : public vtkPolyDataAlgorithm
{
public:
	static vtkSQLiteReader* New();
	vtkTypeRevisionMacro(vtkSQLiteReader,vtkPolyDataAlgorithm);
	// Set/Get the name of the file from which to read points.
	vtkSetStringMacro(FileName);
 	vtkGetStringMacro(FileName);

	vtkSetMacro(HighlightSnapshot,int);
	vtkGetMacro(HighlightSnapshot,int);

	vtkSetMacro(HighlightTrack,int);
	vtkGetMacro(HighlightTrack,int);

//BTX
protected:
	vtkSQLiteReader();
	~vtkSQLiteReader();

	int RequestInformation(vtkInformation*,	vtkInformationVector**,
		vtkInformationVector*);

	int RequestData(vtkInformation*,vtkInformationVector**,
		vtkInformationVector*);

	
// old stuff - can be deleted (with caution!!)
	// functions

	// structs
	struct velocity {
		double vx;
		double vy;
		double vz;
	};

	struct snapshot {
		vtkSmartPointer<vtkPoints> coord;
		vtkSmartPointer<vtkCellArray> cells;
		std::vector<velocity> velo;
		// std::vector<int> npart;
		// std::vector<int> Mvir;
		// std::vector<int> Rvir;
		// ...
	};

	struct snapinfo {
		int snap_id;
		double redshift;
		double time;
		int npart;
	};

	struct trackPoint {
		int snap_id;
		int qid;
	};

	struct track {
		std::vector<trackPoint> point;
		int noOfPoints;
		vtkSmartPointer<vtkPolyLine> line;
	};

	// variables
	int numSnaps;
	int numTracks;
	int totNumPoints;
	int DisplaySnapshot;
	std::vector<vtkSmartPointer<vtkPolyData>> data;
	
	std::vector<snapshot> data2;
	std::vector<snapinfo> snapinfoVector;
	std::vector<track> trackVector;

// --- v3 ------
	//functions

	//structs
	struct SnapshotInfo3 {
		vtkIdType Offset; //stores the id where this snapshot starts
		int lenght; // stores the amount of halos in this snapshot
		double redshift;
		double time;
		int npart;
	};

	//variables
	vtkSmartPointer<vtkPoints>			Position;
	vtkSmartPointer<vtkFloatArray>		Velocity;
	vtkSmartPointer<vtkCellArray>		Cells;
	vtkSmartPointer<vtkCellArray>		Tracks;
	vtkSmartPointer<vtkIdTypeArray>		TrackId;
	vtkSmartPointer<vtkIdTypeArray>		GId;
	vtkSmartPointer<vtkIdTypeArray>		SnapId;
	vtkSmartPointer<vtkFloatArray>		RVir;
	vtkSmartPointer<vtkUnsignedCharArray> colors;

	int nParticles3;
	int nTracks3;

	std::vector<SnapshotInfo3> SnapInfo3;

	char* FileName;

	//gui variables
	int HighlightSnapshot;
	int HighlightTrack;


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

private:
	vtkSQLiteReader(const vtkSQLiteReader&);  // Not implemented.
	void operator=(const vtkSQLiteReader&);  // Not implemented.


// used in v3
	//variables
	sqlite3 * db;
	bool dataIsRead;

	//functions
	int vtkSQLiteReader::ReadHeader(vtkInformationVector*); // reads the database header information

	int vtkSQLiteReader::readSnapshots3(); // reads the snapshots
	int vtkSQLiteReader::readSnapshotInfo3(); 
	int vtkSQLiteReader::readTracks3();
	int vtkSQLiteReader::generateColors();

	int openDB(char*);
	vtkStdString vtkSQLiteReader::Int2Str(int);

// old stuff - can be deleted with caution
	int vtkSQLiteReader::readSnapshots(
		std::vector<vtkSmartPointer<vtkPolyData>> *);
	int vtkSQLiteReader::readSnapshots2();
	int vtkSQLiteReader::ReadTracks();
	int vtkSQLiteReader::GenerateTracks();
	int vtkSQLiteReader::CollectLines(vtkSmartPointer<vtkCellArray>*);
	int vtkSQLiteReader::GenerateOutput(vtkPolyData *);

	int RequestDataDemo(vtkInformationVector*);
	int vtkSQLiteReader::ReadSnapshotInfo();

	int vtkSQLiteReader::SQLQuery(vtkStdString, sqlite3_stmt*);


//ETX
};

#endif