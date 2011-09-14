/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkSQLiteReader.cxx,v $

  Copyright (c) Rafael Kueng
  
=========================================================================*/
#include "vtkTimerLog.h"
#include "vtkSQLiteReader.h"
#include "vtkObjectFactory.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkSmartPointer.h"
#include "vtkFloatArray.h"
#include "vtkLine.h"
#include "vtkPolyLine.h"
#include "vtkArrayData.h"
#include "vtkLookupTable.h"
#include "vtkPolyDataMapper.h"
#include "vtkSMProxy.h"
#include "vtkSMProxyProperty.h"
#include "vtkKdTree.h"
#include "vtkTable.h"
#include "vtkUnicodeString.h"

#include <vector>
#include <sstream>
#include <string>


vtkCxxRevisionMacro(vtkSQLiteReader, "$Revision: 1.0.1 $");
vtkStandardNewMacro(vtkSQLiteReader);


//----------------------------------------------------------------------------
// Constructor
vtkSQLiteReader::vtkSQLiteReader()
{
	this->FileName          = 0;		// init filename
	this->SetNumberOfInputPorts(0);   // set no of input files (0 is just fine)
	this->SetNumberOfOutputPorts(3);

	this->db	= NULL;

	
	AllData = vtkSmartPointer<vtkPolyData>::New();
	SelectedData = vtkSmartPointer<vtkPolyData>::New();
	EmptyData = vtkSmartPointer<vtkPolyData>::New();

	TrackData = vtkSmartPointer<vtkPolyData>::New();
	SnapshotData = vtkSmartPointer<vtkPolyData>::New();


	//collect free gui variables into struct, and init
	Gui.DisplayColliding = &this->DisplayColliding;
	Gui.DisplayMerging = &this->DisplayMerging;
	Gui.DisplayInverted = &this->DisplayInverted;

	Gui.LowerLimit = &this->LowerLimit;
	Gui.UpperLimit = &this->UpperLimit;

	*Gui.DisplayColliding = false;
	*Gui.DisplayMerging = false;
	*Gui.DisplayInverted = false;

	*Gui.LowerLimit = 0.001;
	*Gui.UpperLimit = 0.0;


	//init

	// init dataInfo (maybee just call this->reset!??)
	dataInfo.InitComplete = false;
	dataInfo.dataIsRead = false;
	dataInfo.nPoints = 0;
	dataInfo.nSnapshots = 0;
	dataInfo.nTracks = 0;
	dataInfo.nSnapshots = 0;
	dataInfo.Hubble = 0;

	this->collisionCalc.MergingTolerance = *this->Gui.LowerLimit;
	this->collisionCalc.CollisionTolerance = *this->Gui.UpperLimit;
	this->collisionCalc.CalcIsDone = false;

}

//----------------------------------------------------------------------------
// Deconstructor
vtkSQLiteReader::~vtkSQLiteReader()
{
	sqlite3_close(this->db);
}

//----------------------------------------------------------------------------
// Print self
void vtkSQLiteReader::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
	//TODO be more talkative...
}


//----------------------------------------------------------------------------
// Output port
int vtkSQLiteReader::FillOutputPortInformation(int port,	vtkInformation* info)
{
	if(port == 0)
	{info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");}
	
	else if(port == 1)
	{info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");}
	
	else if(port == 2)
	{info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");}

	return 1;
}

/*----------------------------------------------------------------------------
reads some head data in, is called directly when loading the file (database), before apply is clicked (presents the infos for the window)
this currently does the following:
- checks the database format (right tables, right columns, important data not empty)
- checks what data is available and allocs arrays
*/
int vtkSQLiteReader::RequestInformation(
		vtkInformation* vtkNotUsed(request),
		vtkInformationVector** vtkNotUsed(inputVector),
		vtkInformationVector* outputVector)
{
	/* Stuff for doing it in parallel, leave it for the moment...
	vtkInformation* outInfo = outputVector->GetInformationObject(0);
	// means that the data set can be divided into an arbitrary number of pieces
	outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(),-1);
	*/

	if(!this->dataInfo.InitComplete)
	{
		vtkSmartPointer<vtkTimerLog> inittimer = vtkSmartPointer<vtkTimerLog>::New();
		inittimer->StartTimer();

		// opening database
		if(!this->openDB(this->FileName)){return 0;}

		//init data
		vtkSmartPointer<vtkPoints> points;
		points = vtkSmartPointer<vtkPoints>::New();
		this->AllData->SetPoints(points);
		points = vtkSmartPointer<vtkPoints>::New();
		this->TrackData->SetPoints(points);
		points = vtkSmartPointer<vtkPoints>::New();
		this->SnapshotData->SetPoints(points);

		// add hepler arrays
		this->dataInfo.pSnapIdArray = this->CreateIdTypeArray(this->AllData,"SnapId");
		this->dataInfo.pGIdArray = this->CreateIdTypeArray(this->AllData,"GId");
		this->dataInfo.pTrackIdArray = this->CreateIdTypeArray(this->AllData,"TrackId");
		this->dataInfo.dataArrayOffset = 3;

		// read in database header
		if(!this->ReadHeader()){return 0;}

		//add aditional arrays
		this->CreateArray(this->AllData, "Cv"); //halo concentration
		this->CreateArray(this->AllData, "Redshift");
		this->CreateArray(this->AllData, "a"); //scale factor = 1/(1+redshift)


		this->SelectedData->DeepCopy(this->AllData);
		this->EmptyData->DeepCopy(this->AllData);
		
		this->dataInfo.InitComplete = true;

		//init timers
		this->timing.tot = 0;
		this->timing.read = 0;
		this->timing.calc = 0;
		this->timing.refresh = 0;
		this->timing.generateSelect = 0;
		this->timing.executeSelect = 0;

		inittimer->StopTimer();
		//vtkErrorMacro("Initialisation took: " << inittimer->GetElapsedTime() << " s");
		this->timing.init = inittimer->GetElapsedTime();
	}
	return 1;
}

/*----------------------------------------------------------------------------
is called after clicking apply, reads the actual, selected data
*/
int vtkSQLiteReader::RequestData(vtkInformation*,
								 vtkInformationVector**,
								 vtkInformationVector* outputVector)
{
	vtkSmartPointer<vtkTimerLog> tottimer = vtkSmartPointer<vtkTimerLog>::New();
	tottimer->StartTimer();

	vtkPolyData * out = vtkPolyData::GetData(outputVector->GetInformationObject(0));
	vtkPolyData * snapshotout = vtkPolyData::GetData(outputVector->GetInformationObject(1));
	vtkPolyData * trackout = vtkPolyData::GetData(outputVector->GetInformationObject(2));

	if (!this->dataInfo.dataIsRead)
	// only read if it's not been read before
	{
		vtkSmartPointer<vtkTimerLog> readtimer = vtkSmartPointer<vtkTimerLog>::New();
		readtimer->StartTimer();

		readSnapshotInfo();
		readSnapshots();
		readTracks();
		calculatePointData();

		this->dataInfo.dataIsRead = true;

		readtimer->StopTimer();
		//vtkErrorMacro(" reading data took: " << readtimer->GetElapsedTime() << " s");
		this->timing.read = readtimer->GetElapsedTime();
	}

	if (this->collisionCalc.MergingTolerance != *this->Gui.LowerLimit ||
		this->collisionCalc.CollisionTolerance != *this->Gui.UpperLimit)
	{
		this->collisionCalc.MergingTolerance = *this->Gui.LowerLimit;
		this->collisionCalc.CollisionTolerance = *this->Gui.UpperLimit;
		this->collisionCalc.CalcIsDone = false;
	}

	if((!this->collisionCalc.CalcIsDone) && (*Gui.DisplayColliding || *Gui.DisplayMerging))
	{
		vtkSmartPointer<vtkTimerLog> calctimer = vtkSmartPointer<vtkTimerLog>::New();
		calctimer->StartTimer();

		this->collisionCalc.MergingTolerance = *this->Gui.LowerLimit;
		this->collisionCalc.CollisionTolerance = *this->Gui.UpperLimit;
		findCollisions(&this->collisionCalc);
		this->collisionCalc.CalcIsDone = true;

		calctimer->StopTimer();
		//vtkErrorMacro(" calculations took: " << calctimer->GetElapsedTime() << " s");
		this->timing.calc = calctimer->GetElapsedTime();
	}

	// refresh the ouput
	vtkSmartPointer<vtkTimerLog> refreshtimer = vtkSmartPointer<vtkTimerLog>::New();
	refreshtimer->StartTimer();

	if(!*Gui.DisplayColliding && !*Gui.DisplayMerging && !*Gui.DisplayInverted)
	//nothing selected = display all points
	{
		out->DeepCopy(this->AllData);
	}
	else if (!*Gui.DisplayColliding && !*Gui.DisplayMerging && *Gui.DisplayInverted)
	//dont display anything..
	{
		out->DeepCopy(this->EmptyData);
	}
	else
	{
		vtkSmartPointer<vtkTimerLog> genStimer = vtkSmartPointer<vtkTimerLog>::New();
		genStimer->StartTimer();

		generateSelection(&this->collisionCalc, &this->selection);

		genStimer->StopTimer();
		//vtkErrorMacro(" generating selection: " << genStimer->GetElapsedTime() << " s");
		this->timing.generateSelect = genStimer->GetElapsedTime();

		vtkSmartPointer<vtkTimerLog> seldtimer = vtkSmartPointer<vtkTimerLog>::New();
		seldtimer->StartTimer();

		generateSelectedData(this->SelectedData, &this->selection);

		seldtimer->StopTimer();
		//vtkErrorMacro(" generating selected data: " << seldtimer->GetElapsedTime() << " s");
		this->timing.executeSelect = seldtimer->GetElapsedTime();

		out->DeepCopy(this->SelectedData);
	}
	
	refreshtimer->StopTimer();
	//vtkErrorMacro(" refreshing took: " << refreshtimer->GetElapsedTime() << " s");
	this->timing.refresh = refreshtimer->GetElapsedTime();

	// set additional data output
	snapshotout->DeepCopy(this->SnapshotData);
	trackout->DeepCopy(this->TrackData);

	tottimer->StopTimer();
	//vtkErrorMacro("Total elapsed time: " << tottimer->GetElapsedTime() << " s");
	this->timing.main = tottimer->GetElapsedTime();
	this->timing.tot = this->timing.main + this->timing.init;


		vtkstd::stringstream ss;
		ss<<"\n\nSQLite2Reader was successful!\n\n";
		ss<<"   Headerdata looked like\n";
		ss<<"      Points                : ";
		ss<<this->initData.nPoints<<"\n";
		ss<<"      Tracks                : ";
		ss<<this->initData.nTracks<<"\n";
		ss<<"      Snapshots             : ";
		ss<<this->initData.nSnapshots<<"\n\n";
		ss<<"   actually read in\n";
		ss<<"      Points                : ";
		ss<<this->dataInfo.nPoints<<"\n";
		ss<<"      Tracks                : ";
		ss<<this->dataInfo.nTracks<<"\n";
		ss<<"      Snapshots             : ";
		ss<<this->dataInfo.nSnapshots<<"\n\n";
	if(this->collisionCalc.CalcIsDone){
		ss<<"   calculated\n";
		ss<<"      colliding tracks      : ";
		ss<<this->collisionCalc.CollisionTrackIds->GetNumberOfIds()<<"\n";
		ss<<"      merging tracks        : ";
		ss<<this->collisionCalc.MergingTrackIds->GetNumberOfIds()<<"\n\n";
		ss<<"   selected\n";
		ss<<"      selected points       : ";
		ss<<this->selection.Points->GetNumberOfIds()<<"\n\n";
	}
		ss<<"   Timings\n";
		ss<<"      Init                  : ";
		ss<<this->timing.init<<" s\n";
		ss<<"      Runtime mainmodule    : ";
		ss<<this->timing.main<<" s\n";
		ss<<"        Reading             : ";
		ss<<this->timing.read<<" s\n";
		ss<<"        Calculations        : ";
		ss<<this->timing.calc<<" s\n";
		ss<<"        generating Selection: ";
		ss<<this->timing.generateSelect<<" s\n";
		ss<<"        executing Selection : ";
		ss<<this->timing.executeSelect<<" s\n";
		ss<<"        refreshing Output   : ";
		ss<<this->timing.refresh<<" s\n\n";
		ss<<"      TOTAL TIME TAKEN      : ";
		ss<<this->timing.tot<<" s\n";

	vtkErrorMacro(<<ss.str());


	return 1;
}

/*----------------------------------------------------------------------------
opens the database
	arguments:
		char * filename: path to db
	returns:
		int	errorcode (1 =ok)
	sets:
		sqlite3*	db:		database handle
*/
int vtkSQLiteReader::openDB(char* filename)
{
	if(!filename) //check for existing filename
	{
		vtkErrorMacro("A FileName must be specified.");
		return 0;
	}
	
	if (sqlite3_open(this->FileName, &db) != SQLITE_OK) // open db, returns SQLITE_OK==0 if successful
	{
		vtkErrorMacro("Can't open database: " + *sqlite3_errmsg(db));
		return 0;
	}

	//vtkErrorMacro("opened successfully: " << this->FileName)
	return 1;
}


/*----------------------------------------------------------------------------
Converts a Integer to a string
*/
vtkStdString vtkSQLiteReader::Int2Str(int number)
{
	std::stringstream ss;//create a stringstream
	ss << number;//add number to the stream
	return ss.str();//return a string with the contents of the stream
}

/*----------------------------------------------------------------------------
Reads the header data from db, sets up data structres, checks datastructres
	assumes:
		opened database, db set
	arguments:
		vtkInformationVector Output information
	returns:
		int	errorcode (1 =ok)
	sets:
		numSnap		Number of snapshots

*/
int vtkSQLiteReader::ReadHeader()
{
	// set up local vars
	int counter;
	vtkstd::string name;

	// set up the sql stuff
	vtkStdString	sql_query;	// a querry
	sqlite3_stmt    *res;		// result of sql query
	const char      *tail;		// ???
	int				sql_error;	// return value of sql query
	
	//------------------------------------------------------
	// read the header table
	sql_query = "SELECT * FROM header";
	sql_error = sqlite3_prepare_v2(db,sql_query,1000, &res, &tail);
	if (sql_error != SQLITE_OK){vtkErrorMacro("sqlerror:\nQuerry: "+sql_query+"\nError: "<<sql_error);return 0;}

	if (sqlite3_step(res) == SQLITE_ROW) //there should be only one row..
	{
		/* assume header has always the same structre:
		//TODO chlean this up, autodetection!!
		0 name
		1 format
		2 version
		3 revision
		4 type
		5 Omega0
		6 OmegaLambda
		7 Hubble
		8 ndims
		9 nsnapshots
		10 nsteps

		read (some of) the data in, add more if needed (don't forget to add them in the header, DataInformation struct)
		*/

		// 5 Omega0
		name = sqlite3_column_name(res,5);
		if (name.compare("Omega0")==0)
		{
			// check if set, otherwise use standard value
			double value = sqlite3_column_int(res, 5);
			if (value>0){this->dataInfo.Omega0 = value;}
			else {this->dataInfo.Omega0 = value;} // you could define a standart Omega0 to use here
		}

		// 6 OmegaLambda
		name = sqlite3_column_name(res,6);
		if (name.compare("OmegaLambda")==0)
		{
			// check if set, otherwise use standard value
			double value = sqlite3_column_int(res, 6);
			if (value>0){this->dataInfo.OmegaLambda = value;}
			else {this->dataInfo.OmegaLambda = value;} // you could define a standart OmegaLambda to use here
		}

		// 7 Hubble
		name = sqlite3_column_name(res,7);
		if (name.compare("Hubble")==0)
		{
			// check if set, otherwise use standard value
			double value = sqlite3_column_int(res, 7);
			if (value>0){this->dataInfo.Hubble = value;}
			else {this->dataInfo.Hubble = 70.4;}
		}

		// 9 nsnapshots
		/* not used, instead nSnapshots is gotten by reading snpinfo table
		// TODO use this if available, else use snapinfo
		int num = sqlite3_column_int(res, 9);
		if (num>0){this->dataInfo.nSnapshots = num;}
		else {this->dataInfo.nSnapshots = -1;}
		*/
	}
	else
	{
		vtkErrorMacro("error reading header table");
		return 0;
	}

	//------------------------------------------------------
	// read snapinfo
	sql_query = "SELECT * FROM snapinfo";
	sql_error = sqlite3_prepare_v2(db,sql_query,1000, &res, &tail);
	if (sql_error != SQLITE_OK){vtkErrorMacro("sqlerror:\nQuerry: "+sql_query+"\nError: "<<sql_error);return 0;}

	counter = 0;
	if (sqlite3_step(res) == SQLITE_ROW)
	{
		// read in snapinfo table header
		this->dataInfo.SnapinfoDataColumns.clear();
		this->dataInfo.SnapinfoSnapidColumn = -1;

		for (int i = 0; i<sqlite3_data_count(res); i++)
		{
			name = sqlite3_column_name(res,i);
			const unsigned char * ptr = sqlite3_column_text(res, i);
			if (name.compare("snap_id")==0 && ptr != NULL)
			{
				this->dataInfo.SnapinfoSnapidColumn = i;
			}
			else if (ptr != NULL)
			{
				this->dataInfo.SnapinfoDataColumns.push_back(i);
				this->CreateArray(this->SnapshotData, name.c_str());
			}
		}
		
		//count the num of snapshots
		do {++counter;}	while (sqlite3_step(res) == SQLITE_ROW);

		if(counter!=0)
		{
			this->dataInfo.nSnapshots = counter;
			//this->InitAllArrays(this->SnapshotData,counter);
		}
		else
		{
			vtkErrorMacro("wrong db structre (table 'snapinfo': no data)");
		}
	}
	// something went wrong
	else
	{
		vtkErrorMacro("wrong db structre (table 'snapinfo': querry error)")
		sqlite3_finalize(res);
		return 0;
	}
	sqlite3_finalize(res);

	//------------------------------------------------------
	// read tracks
	sql_query = "SELECT * FROM tracks";
	sql_error = sqlite3_prepare_v2(db,sql_query,1000, &res, &tail);
	if (sql_error != SQLITE_OK){vtkErrorMacro("sqlerror:\nQuerry: "+sql_query+"\nError: "<<sql_error);return 0;}

	int id, gid;
	int max=0, min=0, gidmax=0;
	if (sqlite3_step(res) == SQLITE_ROW)
	{
		vtkstd::string id1 = sqlite3_column_name(res,0);
		vtkstd::string id2 = sqlite3_column_name(res,1);
		vtkstd::string id3 = sqlite3_column_name(res,2);
		if (id1.compare("id")==0 &&
			id2.compare("snap_id")==0 &&
			id3.compare("gid")==0)
		{
			// if necessairy, implement here automatic column detection...
			// for now just assume id in col 0, snap_id in col 1, gid in col 2
			this->dataInfo.TracksTrackidColumn = 0;
			this->dataInfo.TracksSnapidColumn = 1;
			this->dataInfo.TracksGidColumn = 2;

			do
			{
				id = sqlite3_column_int(res,0);
				if(id>max){max=id;}
				if(id<min){min=id;}
				gid = sqlite3_column_int(res,2);
				if(gid>gidmax){gidmax=gid;}

			} while (sqlite3_step(res) == SQLITE_ROW);
			this->dataInfo.gidmax = gidmax;
		}
		// something wrong with db structre
		else
		{
			vtkErrorMacro("wrong db structre (table 'tracks': no id / snap_id / gid columns)");
			return 0;
		}
	}
	
	if (max!=0)
	{
		this->dataInfo.nTracks = max-min+1;
	}
	// something went wrong
	else
	{
		vtkErrorMacro("wrong db structre (table 'tracks': no data)");
		sqlite3_finalize(res);
		return 0;
	}
	sqlite3_finalize(res);

	//------------------------------------------------------
	// count the points
	sql_query = "SELECT * FROM stat";
	sql_error = sqlite3_prepare_v2(db,sql_query,1000, &res, &tail);
	if (sql_error != SQLITE_OK){vtkErrorMacro("sqlerror:\nQuerry: "+sql_query+"\nError: "<<sql_error);return 0;}

	counter = 0;

	if (sqlite3_step(res) == SQLITE_ROW)
	{
		this->dataInfo.StatCordinateColumns.clear();
		this->dataInfo.StatCordinateColumns.resize(3);
		for (int i = 0; i<sqlite3_data_count(res); i++)
		{
			name = sqlite3_column_name(res,i);
			const unsigned char * ptr = sqlite3_column_text(res, i);

			if(name.compare("Xc")==0) {this->dataInfo.StatCordinateColumns.at(0) = i;}
			else if(name.compare("Yc")==0) {this->dataInfo.StatCordinateColumns.at(1) = i;}
			else if(name.compare("Zc")==0) {this->dataInfo.StatCordinateColumns.at(2) = i;}
			else if(name.compare("snap_id")==0) {this->dataInfo.StatSnapidColumn = i;}
			else if(name.compare("gid")==0) {this->dataInfo.StatGidColumn = i;}
			/* not used, instead just get simple colums from velo vector...
			else if(name.compare("VXc")==0) {this->dataInfo.StatCordinateColumns.at(0) = i;}
			else if(name.compare("VYc")==0) {this->dataInfo.StatCordinateColumns.at(1) = i;}
			else if(name.compare("VZc")==0) {this->dataInfo.StatCordinateColumns.at(2) = i;}
			*/
			else if (ptr != NULL)
			{
				this->dataInfo.StatDataColumns.push_back(i);
				this->CreateArray(this->AllData, name.c_str());
			}
		}

		//count the halos
		do{counter++;}
		while (sqlite3_step(res) == SQLITE_ROW);

		this->dataInfo.nPoints = counter;
	} 
	else
	{
		vtkErrorMacro("wrong db structre (table 'stat': no data)");
		sqlite3_finalize(res);
		return 0;
	}
	sqlite3_finalize(res);


	//------------------------------------------------------
	// cleaning up
	/*vtkErrorMacro(" ready to read: " <<
		this->dataInfo.nPoints << " Points in "<<
		this->dataInfo.nSnapshots << " Snapshots, on "<<
		this->dataInfo.nTracks << " tracks.");*/
	this->initData.nPoints = this->dataInfo.nPoints;
	this->initData.nSnapshots = this->dataInfo.nSnapshots;
	this->initData.nTracks = this->dataInfo.nTracks;
	return 1;
}

/*----------------------------------------------------------------------------
reads the snapshots in (reads all the points, and the according data, generates the cells)
	assumes:
		db set and openend
	arguments:
		none
	sets:

	returns:
		int	errorcode (1 = ok)
*/
int vtkSQLiteReader::readSnapshots()
{
	// sql stuff
	vtkStdString	sql_query;	// sql query
	sqlite3_stmt    *res;		// result of sql query
	const char      *tail;		// ???
	int				sql_error;	// return value of sql query

	// Prepare the query
	sql_query = "SELECT * FROM stat";// ORDER BY snap_id";
	sql_error = sqlite3_prepare_v2(db, sql_query, 1000, &res, &tail);
	if (sql_error != SQLITE_OK){vtkErrorMacro("sqlerror:\nQuerry: "+sql_query+"\nError: "<<sql_error);return 0;}

	//shotcuts
	int					xC = this->dataInfo.StatCordinateColumns.at(0);
	int					yC = this->dataInfo.StatCordinateColumns.at(1);
	int					zC = this->dataInfo.StatCordinateColumns.at(2);
	int					SnapidC = this->dataInfo.StatSnapidColumn;
	int					GIdC = this->dataInfo.StatGidColumn;
	vtkstd::vector<int>	*dC = &this->dataInfo.StatDataColumns;

	vtkIdTypeArray		*SnapidA = this->dataInfo.pSnapIdArray;// pointer to snapid array
	vtkIdTypeArray		*GIdA = this->dataInfo.pGIdArray;// pointer to GId array
	int					offsetA = this->dataInfo.dataArrayOffset; // wich is the first of the dataarrays.

	vtkstd::vector< vtkSmartPointer<vtkIdList> > * snapI = &this->SnapInfo2;
	vtkPointData		*pData = this->AllData->GetPointData();
	vtkPoints			*pnts = this->AllData->GetPoints();

	int					gidmax = this->dataInfo.gidmax;

	// prepare variables
	snapI->resize(this->dataInfo.nSnapshots);
	for (int i = 0; i<this->dataInfo.nSnapshots; ++i)
	{
		//snapI->at(i).PointIds.resize(gidmax+1,-1);
		snapI->at(i) = vtkSmartPointer<vtkIdList>::New();
		snapI->at(i)->SetNumberOfIds(gidmax+1);
		for (vtkIdType j = 0; j<gidmax+1;j++){snapI->at(i)->SetId(j,-1);}
	}
	pnts->SetNumberOfPoints(this->dataInfo.nPoints);
	this->InitAllArrays(this->AllData,this->dataInfo.nPoints);


	// loop variables
	int counter = 0;
	int i, snapid, gid;
	double x,y,z;


	while (sqlite3_step(res) == SQLITE_ROW)
	{
		x = sqlite3_column_double(res, xC);
		y = sqlite3_column_double(res, yC);
		z = sqlite3_column_double(res, zC);

		if (x==0.0 && y==0.0 && z==0.0){continue;}

		pnts->InsertPoint(counter,x,y,z);

		snapid = sqlite3_column_int(res, SnapidC);
		gid = sqlite3_column_int(res, GIdC);

		SnapidA->InsertTuple1(counter,snapid);
		GIdA->InsertTuple1(counter,gid);

		//snapI->at(snapid).PointIds.at(gid) = counter;
		snapI->at(snapid)->SetId(gid, counter);

		for (i=0;i<dC->size();++i)
		{
			pData->GetArray(i+offsetA)->InsertTuple1(counter,sqlite3_column_double(res, dC->at(i)));
		}
		++counter;
	}

	pnts->SetNumberOfPoints(counter);
	for (i=0; i<pData->GetNumberOfArrays(); ++i)
	{
		pData->GetArray(i)->SetNumberOfTuples(counter);
	}
	this->dataInfo.nPoints = counter;

	// Create the vertices (one point per vertex, for easy display)
	vtkIdType N = counter;
	vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New();
	vtkIdType *pverts = verts->WritePointer(N, N*2);
    
	for (vtkIdType i=0; i<N; ++i)
	{
		pverts[i*2]   = 1;
		pverts[i*2+1] = i;
    }
	this->AllData->SetVerts(verts);

	this->dataInfo.dataIsRead = true;
	return 1;
}
/*----------------------------------------------------------------------------
reads the tracks, generates the lines
	assumes:
		db set and openend
	arguments:
		none
	sets:
		this->Tracks
		this->TrackId
	returns:
		int	errorcode (1 = ok)
*/
int vtkSQLiteReader::readTracks(){

	//set up sql stuff
	vtkStdString	sql_query;	// sql query
	sqlite3_stmt    *res;		// result of sql query
	const char      *tail;		// ???
	int				sql_error;	// return value of sql query

	//prepare sql querry
	sql_query = "SELECT * FROM tracks";
	sql_error = sqlite3_prepare_v2(db, sql_query, 1000, &res, &tail);
	if (sql_error != SQLITE_OK){vtkErrorMacro("sqlerror:\nQuerry: "+sql_query+"\nError: "<<sql_error);return 0;}

	//shotcuts
	int TrackidC = this->dataInfo.TracksTrackidColumn;
	int SnapidC = this->dataInfo.TracksSnapidColumn;
	int GidC = this->dataInfo.TracksGidColumn;

	vtkIdTypeArray * TracksA = this->dataInfo.pTrackIdArray;
	vtkstd::vector< vtkSmartPointer<vtkIdList> > * ti = &this->TrackInfo2;
	vtkstd::vector< vtkSmartPointer<vtkIdList> > * si = &this->SnapInfo2;

	//init
	vtkSmartPointer<vtkCellArray> Tracks = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkPolyLine> nextLine;
	//vtkstd::vector<int> nPointsOnTrack;

	ti->clear();
	ti->resize(this->dataInfo.nTracks);
	//nPointsOnTrack.clear();
	//nPointsOnTrack.resize(this->dataInfo.nTracks,0);

	for (int i = 0; i<this->dataInfo.nTracks; ++i)
	{
		//ti->at(i).nPoints=0;
		//ti->at(i).PointIds.clear();
		//ti->at(i).PointIds.resize(this->dataInfo.nSnapshots,-1);
		ti->at(i) = vtkSmartPointer<vtkIdList>::New();
		ti->at(i)->SetNumberOfIds(this->dataInfo.nSnapshots);
		for (vtkIdType j = 0; j<this->dataInfo.nSnapshots;++j){ti->at(i)->SetId(j,-1);}
	}

	TracksA->SetNumberOfTuples(this->dataInfo.nPoints);
	TracksA->FillComponent(0,-1);

	// loop vars
	vtkIdType trackid, snapid, gid, newid;
	//vtkIdList * PointsOnTrack;
	//vtkIdList * Track;

	//read in the data
	while (sqlite3_step(res) == SQLITE_ROW)
	{
		trackid = sqlite3_column_int(res, TrackidC);
		snapid = sqlite3_column_int(res, SnapidC);
		gid = sqlite3_column_int(res, GidC);

		//newid = si->at(snapid).PointIds.at(gid);
		newid = si->at(snapid)->GetId(gid);

		if(newid>-1)
		{
			//ti->at(trackid).PointIds.at(snapid) = newid;
			//ti->at(trackid).nPoints++;
			ti->at(trackid)->SetId(snapid,newid);
			//++nPointsOnTrack.at(trackid);

			TracksA->SetTuple1(newid,trackid);
		}
	}

	//delete the fields with no values (this is necessairy, because the data read in from tracks
	// isnt necessairy ordered.. otherwise this could be done much quicker by chaning the above loop
	// starting from an empty vector and just insert the values...
	for (int i = 0; i<ti->size(); ++i)
	{
		ti->at(i)->DeleteId(-1);
	}


	vtkPoints * tmpPoints = this->TrackData->GetPoints();
	tmpPoints->SetNumberOfPoints(this->dataInfo.nTracks);

	vtkSmartPointer<vtkIntArray> nPointsOnTrackDataArray = vtkSmartPointer<vtkIntArray>::New();
	nPointsOnTrackDataArray->SetNumberOfComponents(1);
	nPointsOnTrackDataArray->SetNumberOfTuples(this->dataInfo.nTracks);
	nPointsOnTrackDataArray->SetName("nPointsOnTrack");

	// you can add here more data if you want...

	for (int i = 0; i<ti->size(); ++i)
	{
		Tracks->InsertNextCell(ti->at(i));

		//Fill Data Arrays here..
		tmpPoints->InsertPoint(i,0,0,i);
		//nPointsOnTrackDataArray->InsertTuple1(i,Track->nPoints);
		nPointsOnTrackDataArray->InsertTuple1(i,ti->at(i)->GetNumberOfIds());
	}

	this->AllData->SetLines(Tracks);

	this->TrackData->GetPointData()->AddArray(nPointsOnTrackDataArray);

	//vtkErrorMacro("track count: " << this->dataInfo.nTracks);
	return 1;

}
/*----------------------------------------------------------------------------
Reads in addidtional infos for the snapshots (redshift, time, npart)
	assumes:
		opened database, db set
	sets:
		this->snapinfo	information about snapshots
	arguments:
		none
	returns:
		int	errorcode (1 =ok)
*/

int vtkSQLiteReader::readSnapshotInfo()
{
	// set up the sql stuff
	vtkStdString	sql_query;	// sql query
	sqlite3_stmt    *res;		// result of sql query
	const char      *tail;		// ???
	int				sql_error;	// return value of sql query

	// Prepare the query
	sql_query = "SELECT * FROM snapinfo ORDER BY snap_id, redshift, time";
	sql_error = sqlite3_prepare_v2(db,sql_query,1000, &res, &tail);
	if (sql_error!=SQLITE_OK){vtkErrorMacro("sqlerror:\nQuerry: "+sql_query+"\nError: "<<sql_error);return 0;}
	
	// init vars
	int snapId;
	int counter = 0;
	int SnapinfoSnapidColumn = this->dataInfo.SnapinfoSnapidColumn;
	bool SnapidExists = true;
	if (SnapinfoSnapidColumn == -1) {SnapidExists = false;}
	vtkPointData * pData = this->SnapshotData->GetPointData();

	this->InitAllArrays(this->SnapshotData,this->dataInfo.nSnapshots);

	// fill points
	vtkSmartPointer<vtkPoints> pnt = vtkSmartPointer<vtkPoints>::New();
	pnt->SetNumberOfPoints(this->dataInfo.nSnapshots);
	for (int i = 0; i<this->dataInfo.nSnapshots;i++)
	{
		pnt->InsertPoint(i,0,0,i);
	}
	this->SnapshotData->SetPoints(pnt);

	// main loop
	while (sqlite3_step(res) == SQLITE_ROW)
	{
		if (SnapidExists){snapId = sqlite3_column_int(res, SnapinfoSnapidColumn);}
		else {snapId = counter;}

		for (int i = 0; i<this->dataInfo.SnapinfoDataColumns.size(); i++)
		{
			pData->GetArray(i)->InsertTuple1(snapId,sqlite3_column_double(res, this->dataInfo.SnapinfoDataColumns.at(i)));
		}
		++counter;
	}

	if (this->dataInfo.nSnapshots != counter) {vtkErrorMacro("reading snapshotInfo: something not good..");}
	//this->dataInfo.nSnapshots = counter;
	//this->ResizeAllArrays(this->SnapshotData,this->dataInfo.nSnapshots);

	return 1;
}


/*----------------------------------------------------------------------------
tries to get an estimate for wich tolerance parameter yields to how many colliding tracks.
prints the endresult with vtkErrorMacro

	assumes:
		
	sets:
		int	nCollisions: Number of collisions
	arguments:
	returns:
		int errorcode
*/
int vtkSQLiteReader::calcTolerance()
{
	return 1;
}


/*----------------------------------------------------------------------------
calculates the distance (squared) between the points with gid1 and 2

	assumes:
		
	sets:
		
	arguments:
		ids of two points
	returns:
		the distance
*/
double vtkSQLiteReader::getDistance2(int gid1, int gid2){

	double x1[3];
	double x2[3];

	this->allData.Position->GetPoint(gid1, x1);
	this->allData.Position->GetPoint(gid2, x2);

	return (x1[0]-x2[0])*(x1[0]-x2[0])+
		(x1[1]-x2[1])*(x1[1]-x2[1])+
		(x1[2]-x2[2])*(x1[2]-x2[2]);
}

/*----------------------------------------------------------------------------
Finds the collision and merging Points

	assumes:
		
	sets:
		int	nCollisions: Number of collisions
	arguments:
		pCollCalc		pointer to struct with settings and to fill the data in
	returns:
		int errorcode
*/
int vtkSQLiteReader::findCollisions(CollisionCalculation* pCollCalc){

	//double candidate[3];
	//double dist;
	//double CollTol = pCollCalc->CollisionTolerance;
	//double MergTol = pCollCalc->MergingTolerance;

	////int offset, length;
	//int PointId, SnapId, TrackId, useid;
	vtkIdType SnapId;

	//int nMerging = 0;
	//int nColliding = 0;

	vtkSmartPointer<vtkPoints> points;
	//vtkSmartPointer<vtkKdTree> kDTree;
	vtkSmartPointer<vtkIdList> IdsOfSnap = vtkSmartPointer<vtkIdList>::New();
	//vtkSmartPointer<vtkIdList> IdList;

	/* dont use this, use id list instead..
	// vectores to store the points/tracks / snaps of interesst (indexed by id)
	std::vector<vtkIdType> Points; //points of interest
	std::vector<vtkIdType> Tracks; // tracks of interest
	std::vector<vtkIdType> Snapshots; // 

	// stores info about tracks and points
	//TODO check whtas faster, vector or id list...
	//  -1 unchecked
	//	*0 checked, no collision
	//	1 collision on this point / track
	//	2 merging event (a smaller halo merges in thisone, track has such an event
	//	*3 is the minor of two mergin halos (points only)
	Points.resize(this->dataInfo.nPoints,-1);
	Tracks.resize(this->dataInfo.nTracks,-1);
	Snapshots.resize(this->dataInfo.nSnapshots,-1);
	*/
	vtkSmartPointer<vtkIdList> CollisionPoints = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> CollisionTracks = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> MergingPoints = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> MergingTracks = vtkSmartPointer<vtkIdList>::New();


	for (SnapId=0; SnapId < this->dataInfo.nSnapshots; SnapId++)
	{
		// filter out snapshots with 0 points and possibly other nasty stuff..
		IdsOfSnap->DeepCopy(this->SnapInfo2.at(SnapId));
		IdsOfSnap->DeleteId(-1);
		if (IdsOfSnap->GetNumberOfIds() <= 1){continue;} // dont do this for only 0 or 1 point

		// select points
		points = vtkSmartPointer<vtkPoints>::New();
		this->AllData->GetPoints()->GetPoints(IdsOfSnap,points);

		// build kd tree
		vtkSmartPointer<vtkKdTree> kDTree = vtkSmartPointer<vtkKdTree>::New();
		kDTree->BuildLocatorFromPoints(points);

		// check for best variant:
		//  - FindClosestNPoints
		//  - BuildMapForDuplicatePoints 

		vtkIdTypeArray * result;
		vtkIdType pointid;

		// FIND COLLIDING POINTS
		result = kDTree->BuildMapForDuplicatePoints(pCollCalc->CollisionTolerance);

		for (vtkIdType id = 0; id<result->GetNumberOfTuples(); ++id)
		{
			if (id != result->GetTuple1(id))
			{
				//SPEEDUP: maybe just add ids and check later for uniqueness
				CollisionPoints->InsertUniqueId(IdsOfSnap->GetId(id));
				CollisionPoints->InsertUniqueId(IdsOfSnap->GetId(result->GetTuple1(id)));
			}
		}
		
		/*vtkErrorMacro("Printing the result of kd tree:");
		for (int i = 0; i<result->GetNumberOfTuples(); ++i)
		{
			vtkErrorMacro(" " << i << ": " << result->GetTuple1(i));
		}*/

		result->Delete(); // free memory

		// FIND MERGING POINTS
		result = kDTree->BuildMapForDuplicatePoints(pCollCalc->MergingTolerance);

		for (vtkIdType id = 0; id<result->GetNumberOfTuples(); ++id)
		{
			if (id != result->GetTuple1(id))
			{
				pointid = IdsOfSnap->GetId(id);
				MergingPoints->InsertUniqueId(pointid);
				MergingTracks->InsertUniqueId(this->dataInfo.pTrackIdArray->GetTuple1(pointid));

				pointid = IdsOfSnap->GetId(result->GetTuple1(id));
				MergingPoints->InsertUniqueId(pointid);
				MergingTracks->InsertUniqueId(this->dataInfo.pTrackIdArray->GetTuple1(pointid));
			}
		}
		result->Delete(); // free memory

	}
		
	//DEBUG OUTPUT
	/*
	vtkErrorMacro("Printing the merging track ids (before cleanup):");
	for (int i = 0; i<MergingTracks->GetNumberOfIds(); ++i)
	{vtkErrorMacro(" " << i << ": " << MergingTracks->GetId(i));}
	vtkErrorMacro("Printing the merging point ids (before cleanup):");
	for (int i = 0; i<MergingPoints->GetNumberOfIds(); ++i)
	{vtkErrorMacro(" " << i << ": " << MergingPoints->GetId(i));}
	vtkErrorMacro("Printing the colliding point ids (before cleanup):");
	for (int i = 0; i<CollisionPoints->GetNumberOfIds(); ++i)
	{vtkErrorMacro(" " << i << ": " << CollisionPoints->GetId(i));}
	*/



	//OUTPUT COLLIDING POINTS ONLY
	this->IdListComplement(CollisionPoints, MergingPoints);

	// find colliding, but not merging tracks
	vtkIdType trackid;
	for (vtkIdType i = 0; i<CollisionPoints->GetNumberOfIds(); ++i)
	{
		trackid = this->dataInfo.pTrackIdArray->GetTuple1(CollisionPoints->GetId(i));
		if (MergingTracks->IsId(trackid) == -1)
		{
			CollisionTracks->InsertUniqueId(trackid);
		}
	}

	// DEBUG OUTPUT
	/*
	vtkErrorMacro("Printing the merging track ids (after cleanup):");
	for (int i = 0; i<MergingTracks->GetNumberOfIds(); ++i)
	{vtkErrorMacro(" " << i << ": " << MergingTracks->GetId(i));}
	vtkErrorMacro("Printing the Colliding track ids (after cleanup):");
	for (int i = 0; i<CollisionTracks->GetNumberOfIds(); ++i)
	{vtkErrorMacro(" " << i << ": " << CollisionTracks->GetId(i));}
	
	vtkErrorMacro("Printing the merging point ids (after cleanup):");
	for (int i = 0; i<MergingPoints->GetNumberOfIds(); ++i)
	{vtkErrorMacro(" " << i << ": " << MergingPoints->GetId(i));}
	vtkErrorMacro("Printing the colliding point ids (after cleanup):");
	for (int i = 0; i<CollisionPoints->GetNumberOfIds(); ++i)
	{vtkErrorMacro(" " << i << ": " << CollisionPoints->GetId(i));}
	*/
	
	/*
	vtkErrorMacro("No Merging Tracks: " << MergingTracks->GetNumberOfIds());
	vtkErrorMacro("No Collision Tracks: " << CollisionTracks->GetNumberOfIds());
	vtkErrorMacro("No Merging Points: " << MergingPoints->GetNumberOfIds());
	vtkErrorMacro("No Collision Points: " << CollisionPoints->GetNumberOfIds());
	*/

	/*vtkErrorMacro("the selected points ids are: 1 :");
	for (int i = 0; i<this->SnapInfo2.at(SnapId)->GetNumberOfIds(); ++i)
	{
		vtkErrorMacro(" I  " << i << ": id=" << this->SnapInfo2.at(SnapId)->GetId(i));
	}
	vtkErrorMacro("the selected points are::");
	for (int i = 0; i<points->GetNumberOfPoints(); ++i)
	{
		vtkErrorMacro(" II " << i << ": x=" << points->GetData()->GetComponent(i,0));
	}
	vtkErrorMacro("Printing the result of kd tree:");
	for (int i = 0; i<result->GetNumberOfTuples(); ++i)
	{
		vtkErrorMacro(" " << i << ": " << result->GetTuple1(i));
	}*/
		
	pCollCalc->CollisionPointIds = CollisionPoints;
	pCollCalc->CollisionTrackIds = CollisionTracks;
	pCollCalc->MergingPointIds = MergingPoints;
	pCollCalc->MergingTrackIds = MergingTracks;
	pCollCalc->CalcIsDone = true;

	/*
	vtkErrorMacro("Printing the colliding track ids: (findColl fnc)");
	for (int i = 0; i<CollisionTracks->GetNumberOfIds(); ++i)
	{
		vtkErrorMacro(" " << i << ": " << CollisionTracks->GetId(i));
	}
	*/


	return 1;
}
		

/*----------------------------------------------------------------------------

assumes:
		
	sets:
		
	arguments:
		
	returns:
		int errorcode
*/
int vtkSQLiteReader::generateSelection(
										CollisionCalculation * collCalc,
										SelectionStruct* select)
{
	select->Points = vtkSmartPointer<vtkIdList>::New();
	select->Tracks = vtkSmartPointer<vtkIdList>::New();

	vtkSmartPointer<vtkIdList> allTrackIds = vtkSmartPointer<vtkIdList>::New();
	
	// generate a list of all trackids
	//could be done faster, by assuming ids go from 0 to nTracks... just wanna play safe..
	/* cancel this, it's too slow...
	this->IdTypeArray2IdList(allTrackIds, this->dataInfo.pTrackIdArray);
	*/
	allTrackIds->SetNumberOfIds(this->dataInfo.nTracks);
	for (vtkIdType i = 0; i < this->dataInfo.nTracks; ++i)
	{
		allTrackIds->SetId(i,i);
	}
	
	
	vtkErrorMacro("Printing collisiontracks:");
	for (int i = 0; i<collCalc->CollisionTrackIds->GetNumberOfIds(); ++i)
	{vtkErrorMacro(" " << i << ": " << collCalc->CollisionTrackIds->GetId(i)<<"\n");}

	vtkErrorMacro("Printing the selected track ids (be4 union):");
	for (int i = 0; i<select->Tracks->GetNumberOfIds(); ++i)
	{vtkErrorMacro(" " << i << ": " << select->Tracks->GetId(i)<<"\n");}
	

	if(*this->Gui.DisplayColliding)
	{
		select->Tracks = this->IdListUnionUnique(select->Tracks, collCalc->CollisionTrackIds);
	}

	vtkErrorMacro("Printing the selected track ids (after union w colliding):");
	for (int i = 0; i<select->Tracks->GetNumberOfIds(); ++i)
	{vtkErrorMacro(" " << i << ": " << select->Tracks->GetId(i)<<"\n");}


	if(*this->Gui.DisplayMerging)
	{
		select->Tracks = this->IdListUnionUnique(select->Tracks, collCalc->MergingTrackIds);
	}

	vtkErrorMacro("Printing the selected track ids (after union w merging):");
	for (int i = 0; i<select->Tracks->GetNumberOfIds(); ++i)
	{vtkErrorMacro(" " << i << ": " << select->Tracks->GetId(i)<<"\n");}

	if(*this->Gui.DisplayInverted)
	{
		select->Tracks = this->IdListComplement(allTrackIds, select->Tracks);
	}

	vtkErrorMacro("Printing the selected track ids (after invert):");
	for (int i = 0; i<select->Tracks->GetNumberOfIds(); ++i)
	{vtkErrorMacro(" " << i << ": " << select->Tracks->GetId(i));}


	// collect the points, corresponding to the selected tracks
	for (vtkIdType i = 0 ; i < select->Tracks->GetNumberOfIds(); ++i)
	{
		vtkIdList *tmplist = this->TrackInfo2.at(select->Tracks->GetId(i));
		for (vtkIdType j = 0;
			j<tmplist->GetNumberOfIds();
			++j)
		{
			//SPEEDCHECK01
			//select->Points->InsertUniqueId(tmplist->GetId(j));
			select->Points->InsertNextId(tmplist->GetId(j));
		}
	}

	/*vtkErrorMacro("Printing the selected points:");
	for (int i = 0; i<select->Points->GetNumberOfIds(); ++i)
	{vtkErrorMacro(" " << i << ": " << select->Points->GetId(i));}*/


	return 1;
}

/*----------------------------------------------------------------------------

assumes:
		
	sets:
		
	arguments:
		
	returns:
		int errorcode
*/
int vtkSQLiteReader::generateSelectedData(
		vtkSmartPointer<vtkPolyData> output,
		SelectionStruct* selection)
{
	// generate points
	this->AllData->GetPoints()->GetPoints(
		selection->Points,
		output->GetPoints());
	
	// generate cells
	vtkIdType N = selection->Points->GetNumberOfIds();
	vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New();
	vtkIdType *pverts = verts->WritePointer(N, N*2);
    
	for (vtkIdType i=0; i<N; ++i)
	{
		pverts[i*2]   = 1;
		pverts[i*2+1] = i;
    }
	output->SetVerts(verts);

	// generate tracks
	vtkSmartPointer<vtkCellArray> Tracks = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkIdList> nextLine;
	//vtkIdType oldtrackid;
	//vtkIdList * oldTrack;
	
	int offset = 0;

	for (vtkIdType track = 0; track<selection->Tracks->GetNumberOfIds(); ++track)
	{
		//this is a speedhack, not sure how stable it is..
		vtkIdType nPnts = this->TrackInfo2.at(selection->Tracks->GetId(track))->GetNumberOfIds();
		Tracks->InsertNextCell(nPnts);
		for (vtkIdType i = offset; i<offset+nPnts; ++i)
		{
			Tracks->InsertCellPoint(i);
		}
		offset += nPnts;

		/* Here's the stable way to do it..
		oldTrack = this->TrackInfo2.at(selection->Tracks->GetId(track));
		nextLine = vtkSmartPointer<vtkIdList>::New();
		nextLine->SetNumberOfIds(oldTrack->GetNumberOfIds());

		for (vtkIdType pnt = 0;
			pnt < oldTrack->GetNumberOfIds();
			pnt++)
		{
			nextLine->SetId(pnt,
				selection->Points->IsId(	//gets newid
				oldTrack->GetId(pnt))); //gets oldid

			// if you wanna be really sure nothing goes wrong
			//if (newid>-1){nextLine->InsertId(pnt,newid);}
			//else {vtkErrorMacro("oldid not found");}
		}

		Tracks->InsertNextCell(nextLine);
		*/
	}

	output->SetLines(Tracks);


	// generate pointdata
	vtkPointData * PDall = this->AllData->GetPointData();
	vtkPointData * PDout = output->GetPointData();
	int nPnts = selection->Points->GetNumberOfIds();

	for (int i = 0; i<PDall->GetNumberOfArrays(); ++i)
	{
		PDout->GetArray(i)->SetNumberOfTuples(nPnts);
		PDall->GetArray(i)->GetTuples(selection->Points, PDout->GetArray(i));
	}

	return 1;
}


/*----------------------------------------------------------------------------
Creates a new array (no initialisation is done!)
assumes:
		
	sets:
		
	arguments:
		
		
	returns:
		
*/
vtkSmartPointer<vtkFloatArray> vtkSQLiteReader::CreateArray(
		vtkDataSet *output,
		const char* arrayName,
		int numComponents)
{
	vtkSmartPointer<vtkFloatArray> dataArray=vtkSmartPointer<vtkFloatArray>::New();
	dataArray->SetNumberOfComponents(numComponents);
	dataArray->SetName(arrayName);
	output->GetPointData()->AddArray(dataArray);
	return dataArray;
}
/*----------------------------------------------------------------------------
Creates a new array (no initialisation is done!)
assumes:
		
	sets:
		
	arguments:
		
		
	returns:
		
*/
vtkSmartPointer<vtkIntArray> vtkSQLiteReader::CreateIntArray(
		vtkDataSet *output,
		const char* arrayName,
		int numComponents)
{
	vtkSmartPointer<vtkIntArray> dataArray=vtkSmartPointer<vtkIntArray>::New();
	dataArray->SetNumberOfComponents(numComponents);
	dataArray->SetName(arrayName);
	output->GetPointData()->AddArray(dataArray);
	return dataArray;
}
/*----------------------------------------------------------------------------
Creates a new array (no initialisation is done!)
assumes:
		
	sets:
		
	arguments:
		
		
	returns:
		
*/
vtkSmartPointer<vtkIdTypeArray> vtkSQLiteReader::CreateIdTypeArray(
		vtkDataSet *output,
		const char* arrayName,
		int numComponents)
{
	vtkSmartPointer<vtkIdTypeArray> dataArray=vtkSmartPointer<vtkIdTypeArray>::New();
	dataArray->SetNumberOfComponents(numComponents);
	dataArray->SetName(arrayName);
	output->GetPointData()->AddArray(dataArray);
	return dataArray;
}
/*----------------------------------------------------------------------------
Initialises all array
assumes:
		
	sets:
		
	arguments:
		
		
	returns:
		int errorcode
*/
int vtkSQLiteReader::InitAllArrays(
		vtkDataSet *output,
		unsigned long numTuples)
{
	vtkSmartPointer<vtkDataArray> dataArray;

	for (int i = 0; i<output->GetPointData()->GetNumberOfArrays(); ++i)
	{
		dataArray = output->GetPointData()->GetArray(i);
		dataArray->SetNumberOfTuples(numTuples);
		for (int j=0; j < dataArray->GetNumberOfComponents(); ++j) {
			dataArray->FillComponent(j, 0.0);
		}
	}
	return 1;
}

/*----------------------------------------------------------------------------
Unions two IdLists, adds list 2 to list 1
( set notation: list1 u list2, list1 or list2)
assumes:
		
	sets:
		
	arguments:
		
		
	returns:
		int errorcode
*/

vtkSmartPointer<vtkIdList> vtkSQLiteReader::IdListUnion(
	vtkSmartPointer<vtkIdList> list1,
	vtkSmartPointer<vtkIdList> list2)
{
	int n1 = list1->GetNumberOfIds();
	int n2 = list1->GetNumberOfIds();
	list1->SetNumberOfIds(n1+n2);
	for (vtkIdType i = 0; i<n2; ++i)
	{
		list1->SetId(n1+i,list2->GetId(i));
	}
	return list1;
}

/*----------------------------------------------------------------------------
Unions two IdLists, adds list 2 to list 1
and makes sure, each id is unique
list1 for uniqueness!
( set notation: list1 u list2, list1 or list2)
assumes:
		
	sets:
		
	arguments:
		
		
	returns:
		int errorcode
*/

vtkSmartPointer<vtkIdList> vtkSQLiteReader::IdListUnionUnique(
	vtkSmartPointer<vtkIdList> list1,
	vtkSmartPointer<vtkIdList> list2)
{
	vtkSmartPointer<vtkIdList> result = vtkSmartPointer<vtkIdList>::New();
	//result->SetNumberOfIds(list1->GetNumberOfIds()+list2->GetNumberOfIds());
	
	for (vtkIdType i = 0; i<list1->GetNumberOfIds(); ++i)
	{
		//vtkErrorMacro(" union: l1: id-"<<i<<" - "<<list1->GetId(i));
		result->InsertUniqueId(list1->GetId(i));
	}
	for (vtkIdType i = 0; i<list2->GetNumberOfIds(); ++i)
	{
		//vtkErrorMacro(" union: l2: id-"<<i<<" - "<<list2->GetId(i));
		result->InsertUniqueId(list2->GetId(i));
	}
	result->Squeeze();
	//list1 = result;
	return result;

}
/*----------------------------------------------------------------------------
Intersects two IdLists, result in list1 (wrapper)
( set notation: list1 n list2, list1 and list2)
assumes:
		
	sets:
		
	arguments:
		
		
	returns:
		int errorcode
*/

vtkSmartPointer<vtkIdList> vtkSQLiteReader::IdListIntersect(
	vtkSmartPointer<vtkIdList> list1,
	vtkSmartPointer<vtkIdList> list2)
{
	list1->IntersectWith(*list2);
	return list1;
}
/*----------------------------------------------------------------------------
Gets the complement of two IdLists, result in list1
or: deletes all ids from list 2 in list 1
( set notation: list1 \ list2, list1 - list2)
assumes:
		
	sets:
		
	arguments:
		
		
	returns:
		int errorcode
*/

vtkSmartPointer<vtkIdList> vtkSQLiteReader::IdListComplement(
	vtkSmartPointer<vtkIdList> list1,
	vtkSmartPointer<vtkIdList> list2)
{
	for (vtkIdType i = 0; i<list2->GetNumberOfIds(); ++i)
	{
		list1->DeleteId(list2->GetId(i));
	}
	list1->Squeeze();
	return list1;
}

/*----------------------------------------------------------------------------
Converts a vtkidTypeArray to a vtkIdList
(adds only unique ids to the list)
assumes:
		
	sets:
		
	arguments:
		
		
	returns:
		int errorcode
*/

int vtkSQLiteReader::IdTypeArray2IdList(vtkIdList* destination, vtkIdTypeArray* source)
{
	for (vtkIdType i = 0; i<source->GetNumberOfTuples(); ++i)
	{
		destination->InsertUniqueId(source->GetTuple1(i));
	}
	return 1;
}


/*----------------------------------------------------------------------------
Adds/calculates aditional point data
assumes:
		
	sets:
		
	arguments:
		
		
	returns:
		
*/
int vtkSQLiteReader::calculatePointData()
{
	vtkDataArray * cvarr = this->AllData->GetPointData()->GetArray("Cv");
	vtkDataArray * rsarr = this->AllData->GetPointData()->GetArray("Redshift");
	vtkDataArray * aarr = this->AllData->GetPointData()->GetArray("a");

	cvarr->SetNumberOfTuples(this->dataInfo.nPoints);
	rsarr->SetNumberOfTuples(this->dataInfo.nPoints);
	aarr->SetNumberOfTuples(this->dataInfo.nPoints);
	
	//cvarr->FillComponent(0,0);

	vtkDataArray * Vmax = this->AllData->GetPointData()->GetArray("Vmax");
	vtkDataArray * Rmax = this->AllData->GetPointData()->GetArray("Rmax");

	vtkDataArray * ssdrsarr = this->SnapshotData->GetPointData()->GetArray("redshift");

	float cv, rs, a;
	int snapid;

	for(vtkIdType i = 0; i<this->dataInfo.nPoints; ++i)
	{
		if (Rmax->GetTuple1(i) == 0)
		{
			cv = 0;
		}
		else
		{
			cv = 2 * ( Vmax->GetTuple1(i) / (this->dataInfo.Hubble * Rmax->GetTuple1(i))) *
				( Vmax->GetTuple1(i) / (this->dataInfo.Hubble * Rmax->GetTuple1(i)));
		}
		cvarr->SetTuple1(i,cv);

		snapid = this->AllData->GetPointData()->GetArray(this->dataInfo.StatSnapidColumn)->GetTuple1(i);
		
		rs = ssdrsarr->GetTuple1(snapid);
		rsarr->SetTuple1(i,rs);
		
		a=1/(1+rs);
		aarr->SetTuple1(i,a);
	}

	return 1;
}