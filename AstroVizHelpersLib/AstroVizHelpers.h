/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: AstroVizHelpers.h,v $

  Copyright (c) Christine Corbett Moran
  All rights reserved.
     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME AstroVizHelpers
// .SECTION Description
// Helper library for many common tasks in AstroViz
#include "vtkIdTypeArray.h" // TODO: needed to include this, but should figure out how to remove
#include <iostream>
#include <sstream>
class vtkPolyData;
class vtkPointSet;
class vtkDataSet;
class vtkTable;
class vtkFieldData;
class vtkPointData;
class vtkDataSetAttributes;
class vtkIdTypeArray;
class vtkIdList;
class vtkPointLocator;
class vtkCell;
class vtkCellArray;
class vtkFloatArray;
class vtkInformationVector;
class vtkMultiProcessController;
class vtkVariant;

// CGS units needed
static const float kpcInCm=3.08568025*pow(10.0,21.0);
static const float pcInCm=3.08568025* pow(10.0,18.0);
static const float kmInCm=pow(10.0,5.0);
static const float GyrInS=3.1536*pow(10.0,16.0);
static const float yrInS=3.1553*pow(10.0,7.0);
static const float msolInG=1.98892*pow(10.0,33.0);



/*
* The following methods take and modify vtkPolyData
*/
// Description:
// sets the point vertices in the output vector, assigning the point 
// a unique id and places one point per cell. Works only with PolyData
vtkIdType SetPointValue(vtkPolyData* output,float pos[]);

// Description:
// creates a sphere of radius R, center center in the output.
void CreateSphere(vtkPolyData* output,double radius,double center[]);

/*
* The following methods take and modify vtkPointSet data
*/

// Description:
// sets the data value in the output vector in array arrayName at 
// position id to data.
void SetDataValue(vtkPointSet* output, const char* arrayName,
	vtkIdType id,float data[]);
// Description:
// sets the data value in the output vector in array arrayName at 
// position id to data.
void SetDataValue(vtkPointSet* output, const char* arrayName,
	vtkIdType id,double data[]);
// Description:
// sets the data value in the output vector in array arrayName at 
// position id to data.
// SetIdTypeValue performs a safe downcast from vtkDataArray to 
// vtkIdTypeValue, unlike the other SetDataValue methods, which have
// no need for downcast. Thus array requested MUST be a vtkIdType array
void SetIdTypeValue(vtkPointSet* output, const char* arrayName,
	const vtkIdType indexId,const vtkIdType globalId);


// Description:
// points the double* data to the value in the output vector in 
// array arrayName at position id would like to be float, but the 
// version of vtk I am working with does not support 
double* GetDataValue(vtkPointSet* output, const char* arrayName,
	vtkIdType id);

// Description:
// create a vtkDataArray of floats with the  name arrayName, number /
// of components. place it in the vtkPointSet
void AllocateDataArray(vtkPointSet* output, const char* arrayName,
	int numComponents, unsigned long numTuples);

// Description:
// create a vtkDataArray with the  name arrayName, number of components 
// numComponents and number of tuples numTuples of type T. place it in
// the vtkTable
void AllocateDataArray(vtkTable* output, const char* arrayName,
	int numComponents, unsigned long numTuples);

// Description:
// create a vtkDataArray with the  name arrayName, number of components 
// numComponents and number of tuples numTuples of type T. place it in
// the output
void AllocateIdTypeDataArray(vtkPointSet* output, const char* arrayName,
	int numComponents, unsigned long numTuples);
	

// Description:
// create a vtkDataArray with the  name arrayName, number of components 
// numComponents and number of tuples numTuples of type T
// e.g. AllocateDoubleDataArray("density",1,100) creates a array of 100 
// scalar double densities
void AllocateDoubleDataArray(vtkPointSet* output, const char* arrayName,
	int numComponents, unsigned long numTuples);

// Description:
// create a vtkDataArray with the  name arrayName, number of components 
// numComponents and number of tuples numTuples of type T
// e.g. AllocateDoubleDataArray("density",1,100) creates a array of 100 
// scalar double densities
void AllocateIntDataArray(vtkPointSet* output, const char* arrayName,
	int numComponents, unsigned long numTuples);
// Description:
// Given a VTK array, sets it to have arrayName, numComponents and numTuples
// unlike the above methods, does not create the array, nor add it to the 
// output	
void InitializeDataArray(vtkDataArray* dataArray, const char* arrayName,
	int numComponents, unsigned long numTuples);

// Description:
// returns a pointer to the point's coordinates in output which corresponds 
// to this id
double* GetPoint(vtkPointSet* output,vtkIdType id);

// Description:
// takes in a double array of size three representing a point
// and converts it to an array of the same size but in float precision
float* DoublePointToFloat(double point[]);

// Description:
// Creates a new information vector that is a deep copy of the old one
// Note, returned vector of pointers to vtkInformationVectors must be 
// deleted by caller. See DeleteDeepCopyInput below.
vtkInformationVector** DeepCopyInputVector(
	vtkInformationVector** inputVector,unsigned long inputVectorSize); 

// Description: 
// Deletes an array of pointers to vtkInformationVectors
void DeleteDeepCopyInputVector(vtkInformationVector** inputVector, 
	unsigned long inputVectorSize);

// Description:
// converts a double to a vtkstd::string
inline vtkstd::string ToString(double x)
{
  std::ostringstream o;
	o << x;
  return o.str();
}
// Description:
// returns true if this process should be run in parallel
// (i.e. we have  non-null controller and more than one process to work with)
bool RunInParallel(vtkMultiProcessController* controller);
enum PointsInRadiusMPIData
{
	TOTAL_MASS_IN_SPHERE,
	TOTAL_NUMBER_IN_SPHERE,
	MAX_R,
	NUMBER_PARTICLE_OFFSET
};
// Description:
// Uses the Illinois root finding method to find the root of the function
// func. The root must lie between r and s. Root is returned when it is found 
// within the accuracy xacc, yacc. pnIter indicates how many iterations the
// algorithm took to converge. Returns -1 if there are problems finding the 
// root, thus can only be used to find positive roots.
double IllinoisRootFinder(double (*func)(double,void *),void *ctx,double r,double s,double xacc,double yacc,int *pnIter);

// Description:
// Given a point and an input data set, computes the maximum distance from the
// point to the data set boundaries.
// This is the maximum distance from the point to the 8 corners of the 
// bounding cube of this dataset
/*    	  
*										             (xmin,ymax,zmax)
*														       / |\
*														      /  | \
*														     /   |  \
*														   //    |   \\                							^  ^
*														  /      |     \              						z| /y
*														 /       \(xmin,ymax,zmin)    						|/
*			    (xmin,ymin,zmax)  |       / \      | (xmax,ymax,zmax)       \
*														|\     /   \    /|             						x\
*														| \   /  .p \\ / |               						v
*														|  \ /       //  |
*														|  /\       /  \ |
*														| /  \\    /    \|
*				  (xmin,ymin,zmin)	|/     \  /      | (xmax,ymax,zmin)
*														/   (xmax,ymin,zmax)  
*														 \       |      /
*														  \      |     /
*														   \\    |   //
*														     \   |  /
*														      \  | /
*														       \  /
*										            (xmax,ymin,zmin)
*
*/
double ComputeMaxR(vtkPointSet* input,double point[]);

// Description:
// Runs ComputeMaxR on each process, then syncronizes results taking the
// over all global maximum. This result will be slightly different that that
// if run serially, as the bounding boxes may be tighter on subprocesses
// than an overall bounding box. But if anything, the parallel version should
// be an even better approximation to the maximum radius (as we use the
// bounding box, not the individual points therein to compute) than the
// serial version, and is always guaranteed to be greater than or equal
// to the true maximum (if one used the points within the data set)
double ComputeMaxRadiusInParallel(
	vtkMultiProcessController* controller,vtkPointSet* input,double point[]);
	
// Description:
// This assumes locatorStruct is a VirialRadiusInfo struct, which contains
// a locator for a given vtkdataset, a center from which to calculate
// the volume 
// Given a radius, a center, calculates the density of within the sphere
// of radius r around the center and subtracts the critical density of 
// the user's choice.
// 
double OverDensityInSphere(double r, void* locatorStruct);

// Description:
// This assumes locatorStruct is a VirialRadiusInfo struct, which contains
// a locator for a given vtkdataset, a center from which to calculate
// the volume 
// Given a radius, a center, calculates the number density of within the
// sphere of radius r around the center and subtracts the critical density of 
// the user's choice.
double OverNumberInSphere(double r, void* locatorStruct);

// Descriptions:
// If this is called by each process in the domain of the controller, finds
// all points within a given radius on this process
vtkIdList* FindPointsWithinRadius(double r, double* center,
	vtkPointLocator* locatorOfThisProcess);

// Description:
// The VirialRadiusInfo struct is an containing:
// .locator which is a vtkPointLocator
// .center  which is a double[3]
// .criticalDensity which is a double
// .virialRadius
struct VirialRadiusInfo 
{
	vtkPointLocator* locator;
	vtkMultiProcessController* controller;
	double center[3];
	double criticalValue;
	double virialRadius;
	double softening;
	vtkstd::string massArrayName;
};


// Description:
// for processors p=0 to n, computes
//     particle offset = sum(particles handled by pi<p)
unsigned long ComputeParticleOffset(vtkMultiProcessController* controller,vtkPointSet* input);

// Description:
// Computes the virial radius >=0 base upon the user defined 
// overdensity and center. Returns -1 if there is a problem. 
// Works in parallel if a controller is specified not equal to null and if 
// the number of processors is > 1
VirialRadiusInfo ComputeVirialRadius(
	vtkMultiProcessController* controller, vtkPointLocator* locator,
	vtkstd::string massArrayName, double softening,double overdensity,
	double maxR,double center[]);

// Description:
// shifts every item in array one to left (the first element is thrown away)
// then sets inserts updateValue in the last, free slot
template <class T> void shiftLeftUpdate(T* array,int size, T updateValue);

// Description:
// From dataSet, initializes a newDataSet with point and cell arrays (empty),
// and corresponding  data arrays (empty) with identical names, and number of 
// components to dataSet. Then copys, if listed in the vtkIdList,
// from the old data set to the new data set the points and their cell data.
vtkPolyData* CopyPointsAndData(vtkPointSet* dataSet, vtkIdList*
 	pointsInRadius);

// Description:
// Given a populated virialradiusinfo struct, returns a dataset corresponding
// to only those points within the virial radius.
// This method only works if input was vtkPolyData...
vtkPolyData* GetDatasetWithinVirialRadius(VirialRadiusInfo virialRadiusInfo);

// Description
// Given an input data set, the bin number, a list of points in the relevant
// bin,  and the output table, computes the average radial velocity in the
// bin, the radial velocity dispersion in the bin, the tangential velocity
// in the bin, the tangentical velocity dispersion in the bin, and the angular
// momentum in the bin and adds it to the output table in the appropriate bin.
//
// Average velocity: 
// \vec v_{ave} = \frac{\sum_{i=1}^{N} \vec v_i}{N}
//
// Velocity Dispersion:
// \vec \sigma_v = \frac{\sum_{i=1}^{N} 
// \sqrt{(\vec v_i)^2-(\vec v_{ave}^2)}}{N}
//
// Specific Angular Momentum:
// \vec j = \frac{\sum_{i=1}^{N} \vec r_i x \vec v_i}{N}
//
void ComputeStatisticsInBin(vtkPolyData* inputDataSet, double center[],
	int binNum, vtkIdList* pointsInBin, vtkTable* output);
	
// Description:
// given 3-vector vectorOne and 3-vector vectorTwo, computes the 
// projection of vectorOnew onto vector two.
// Can be used to calculate e.g.
// radialVelocity = ComputeProjection(v,r);
//  or velocitySquared = ComputeProjection(v,v);
void ComputeProjection(double  vectorOne[],double vectorTwo[], double result[]);

// Description:
// Computes the vector difference between two 3-vectors
// Can be used to calculate e.g. tangential velocity
// tangentialVelocity = PointVectorDifference(v,radialVelocity);
void PointVectorDifference(double vectorOne[], double vectorTwo[], double result[]);

// Description
// Give a 3 vector v and a three vector r computes the specific angular 
// momentum = r x v
void ComputeAngularMomentum(double v[], double r[], double result[]);

// Description:
// Multiplies in place a 3-vector by a constant
void VecMultConstant(double vector[],double constant);

// Description
// Given a vSquaredAve and a vAve calculates the velocity dispersion
// placing it in the output vector velocityDispersion
void ComputeVelocityDispersion(double *vSquaredAve, double *vAve, double *result);
	
// Description:
// helper function to compute radial velocity
void ComputeRadialVelocity(double v[],double r[], double result[]);

// Description:
// helper function to compute tangential velocity
void ComputeTangentialVelocity(double v[],double r[], double result[]);

// Description:
// helper function to compute angular momentum
void ComputeAngularMomentum(double v[], double r[], double result[]);

// Description:
// helper function to compute velocity squared
double ComputeVelocitySquared(double v[]);

// Description:
// helper function to compute radial velocity squared
double ComputeRadialVelocitySquared(double v[],double r[]);

// Description:
// helper function to compute tangential velocity squared
double ComputeTangentialVelocitySquared(double v[],double r[]);

// Description:
// helper function to compute circular velocity
void ComputeCircularVelocity(double *cumulativeMass, double *binRadius, double result[]);

// Description:
// helper function to compute density
void ComputeDensity(double *cumulativeMass, double *binRadius, double result[]);
	
// Description:
// Helper function to compute the midpoint between two points
void ComputeMidpoint(double pointOne[], double pointTwo[], double result[]);

// Description:
// Given an AdditionalAttributeFile of the format
// int numBodies
// float attribute
// float attribute
// ...
// float attribute
// Given an offset seeks to a given attribute in the file from the current
// position in the file. 
double SeekInAsciiAttributeFile(
	ifstream& asciiFile,const int offset);






