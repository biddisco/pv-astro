/*=========================================================================

  Program:   AstroViz plugin for ParaView
  Module:    $RCSfile: vtkProfileFilter.h,v $

  Copyright (c) Christine Corbett Moran
  All rights reserved.
     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkProfileFilter 
// .SECTION Description
// Calculates various physical profiles as a function of radius.
// Fully parallel.
#ifndef __vtkProfileFilter_h
#define __vtkProfileFilter_h
#include "vtkTableAlgorithm.h" // super class
#include "vtkStringArray.h" // some class variables are vtkStringArrays
#include <vtkstd/vector>
#include "vtkSmartPointer.h"
#include "vtkFloatArray.h"
//
#include <functional>

class vtkPointSet;
class vtkMultiProcessController;
class vtkPlane;
class vtkPoints;

//----------------------------------------------------------------------------

enum ProfileMPIData 
{
	DATA_TABLE
};

enum BinUpdateType
{
	ADD, 
	MULTIPLY,
	SET
};

enum ColumnType
{
	AVERAGE,
	TOTAL,
	CUMULATIVE
};


//----------------------------------------------------------------------------
class VTK_EXPORT vtkProfileFilter : public vtkTableAlgorithm
{
public:
  static vtkProfileFilter* New();
  vtkTypeRevisionMacro(vtkProfileFilter, vtkTableAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  // Description:
  // Get/Set the center
  vtkSetVector3Macro(Center,double);
  vtkGetVectorMacro(Center,double,3);
  // Description:
  // Get/Set the number of bins
  vtkSetMacro(BinNumber,int);
  vtkGetMacro(BinNumber,int);
    
  // Description:
  // Get/Set the Profile Axis
  vtkSetVector3Macro(ProfileAxis,double);
  vtkGetVectorMacro(ProfileAxis,double,3);

  // Description:
  // Get/Set the number of bins
  vtkSetMacro(ProfileHeight,double);
  vtkGetMacro(ProfileHeight,double);

  // Description:
  // Specify the point locations used to probe input. Any geometry
  // can be used. New style. Equivalent to SetInputConnection(1, algOutput).
  void SetSourceConnection(vtkAlgorithmOutput* algOutput);
  // Description:
  // By defualt this filter uses the global controller,
  // but this method can be used to set another instead.
  virtual void SetController(vtkMultiProcessController*);

//BTX
protected:
  vtkProfileFilter();
  ~vtkProfileFilter();

  // Description:
  // This is called within ProcessRequest when a request asks the algorithm
  // to do its work. This is the method you should override to do whatever the
  // algorithm is designed to do. This happens during the fourth pass in the
  // pipeline execution process.
  virtual int RequestData(vtkInformation*, 
		vtkInformationVector**, vtkInformationVector*);
  vtkMultiProcessController *Controller;

  // These are indices into the array containers and must match
  // the order which profiles are created in InitializeBins
  enum ProfileIndex {
    BinRadius=0,
    BinRadiusMin,
    NumberInBinTotal,
    NumberInBinCumulative,
    // then there are N=PointData()->GetNumberOfArrays() from here onwards
    PointDataArrays,
  };
 
	// Description:
	// the ProfileElement protected nested class holds all the information
	// necessary to initialize and at runtime compute the value
	// of an additional profile element. Each of these are initialized
	// via the constructor
	// name : unique string name describing this element, 
	// only the base name, an affix will be added for any
	// quantities desired to be computed
	// number elements : the number of elements in each entry
	// funcPtr : the function to use to evaluate an update
	// average : if 1 compute the average of this quantity
	// total : if 1 compute the total of this quantity
	// cumulative: if 1 compute the cumulative value of this quantity
	// postprocess: if 1 only update during post processing
  //
  // This would be better using boost::function and boost:bind
  //
  typedef void   (*DoubleFunction0)(double [], double [], double []);
  typedef double (*DoubleFunction1)(double [], double []);
  typedef double (*DoubleFunction2)(double []);
//  typedef void   (*PostProcessFunc)(vtkVariant, vtkVariant, double []);
  enum FuncType 
  {
    FUNC0_TYPE=0,
    FUNC1_TYPE,
    FUNC2_TYPE,
    NULL_FUNC,
    POSTPROCESS_TYPE,
  };
  //
  class ProfileElement
  {
  public:
    //
    typedef std::vector<ProfileElement>::iterator PIterator;
    //
    DoubleFunction0 func0;
    DoubleFunction1 func1;
    DoubleFunction2 func2;
    DoubleFunction0 funcP; // can be removed now, use func0
		vtkstd::string  BaseName;
    vtkstd::string  DerivedName;
		int             NumberComponents;
		FuncType        FuntionType;
		ColumnType      ProfileColumnType;
	  PIterator       ArgOne; 
		PIterator       ArgTwo;
    vtkSmartPointer<vtkFloatArray> Data;
    float          *DataPointer;

    // Description:
		// quantities to be processed for each element in each bin with the
		// function *functPtr which takes in a radius and a velocity given
		// by double arrays
		ProfileElement(const char *baseName, int numberComponents,
			DoubleFunction0, ColumnType columnType);

		ProfileElement(const char *baseName, int numberComponents,
			DoubleFunction1, ColumnType columnType);

		ProfileElement(const char *baseName, int numberComponents,
			DoubleFunction2, ColumnType columnType);

		ProfileElement(const char *baseName, int numberComponents,
			ColumnType columnType);

		// Description:
		// if post processing is desired, then must specify two arguments, which
		// are profile elements. Their data (all computed) will be retrieved from
		// the output in postprocessing and a final computation will be performed
		// with functionPtr.
		//
		// The last four arguments specify which two columns
		// data should be handed to the postprocessing function, which
		// takes two vtkVariants as arguments and returns a double*,
		// thus requires that they are part of the input (for which
		//  CUMULATIVE,AVERAGE and TOTAL are computed for each array name)
		// or that they are specified as an additional profile element above
		ProfileElement(const char *baseName, int numberComponents,
			DoubleFunction0,
			PIterator argOne, 
			PIterator argTwo);
		~ProfileElement();
    
    // String Bookkeeping
    void        CreateColumnNames();
    const char *GetColumnName() const;

    // Array Bookkeeping
    vtkFloatArray *AllocateBinArray(vtkIdType numTuples);

    //----------------------------------------------------------------------------
    // A functor we use to locate a given profile object in a list
    //----------------------------------------------------------------------------
    struct FindProfile : std::unary_function<vtkProfileFilter::ProfileElement,bool> {
      FindProfile(std::string name, ColumnType columntype) {
        Name = name;
        Columntype  = columntype;
      }
      bool operator()(const vtkProfileFilter::ProfileElement &lhs) {
        if (std::string(lhs.BaseName)==Name && lhs.ProfileColumnType==Columntype) {
          return true;
        }
        return false;
      }
      std::string Name;
      ColumnType  Columntype;
    };
 	};

	double Delta;
  
  // Description:
	// Center around which to compute radial bins
	double Center[3];
  
  // Description:
	// If non-zero specified, cylindrical profiles are computed about this axis
	double ProfileAxis[3];

  // Description:
  // if specified and non zero, items above this height from center are ignored in profile
	// user input, with default
	double ProfileHeight;

  // Description:
	// Number of bins to use
	// user input, with default
	int BinNumber;

  // Description:
	// Spacing between bins, automatically calculated based upon other
	// selections
	double BinSpacing;

	// Description:
	// Max distance from center point to the data set boundaries, or to
	// the virial radius if applicable
	double MaxR;

  // Description:
	// Profiles we generate from the input
	vtkstd::vector<ProfileElement> ProfileQuantities;

  // Description:
	// Quantities to add to the input
	vtkstd::vector<ProfileElement> AdditionalProfileQuantities;
  
	// Description:
  // Calculates the center and the maximum distance from the center
	// based upon the user's input and the boundaries of the dataset. Also
	// calculates the bin spacing based on the center, maxR and the user's 
	// desired number of bins
	void SetBoundsAndBinExtents(vtkPointSet* input, vtkDataSet* source);

	// Description:
	// SetBoundsAndBinExtents must have been called first.
	// 
	// Initializes a row for each bin, with the max radius of that bin
	// 
	// For each dataarray in the input, define a total and an
	// average column in the binned output table.
	//
	// For each additional quantity as specified by the 
	// AdditionalProfileQuantities array, define a average column in the
	// binned output table
	//
	// For each cumulative array as specified by the CumulativeQuantitiesArray
	//
	void InitializeBins(vtkPointSet* input, vtkTable* output);

	// Description:
	// Calculates the bin spacing 
	double CalculateBinSpacing(double maxR,int binNumber);

	// Description:
	// For each point in the input, update the relevant bins and bin columns
	// to reflect this point's data. Finally compute the averages, relevant
	// dispersions, and global statistics.
	void UpdateStatistics(vtkPointSet* input,vtkTable* output);

	// Description:
	// For each quantity initialized in InitializeBins updates the statistics
	// of the correct bin to reflect the data values of the point identified
	// with pointGlobalId. Note: for quantities that are averages, or require
	// post processing this are updated additively as with totals. This is
	// why after all points have updated the bin statistics,
	// BinAveragesAndPostprocessing must be called to do the proper averaging
	// and/or postprocessing on  the accumlated columns.
	void UpdateBinStatistics(vtkPointSet* input, 
    vtkPoints *points, vtkDataArray *velocity,
		vtkIdType pointGlobalId,vtkTable* output);
	
	// Description:
	// returns the bin number in which this point lies.
	int GetBinNumber(double x[]);
	  
	// Description:
	// Updates the data values of attribute specified in attributeName
	// in the bin specified by binNum, either additively or 	
	// multiplicatively as specified by updatetype by dataToAdd. Calls
	// either UpdateArrayBin or UpdateDoubleBin depending on the
	// type of data in the bin.
	
  void UpdateBin(int binNum, BinUpdateType updateType,
    ProfileElement &profile, double *updateData);

  template<typename T>
  void UpdateBin(BinUpdateType updateType,
    ProfileElement &profile, float *tableData, T *newData);

  void UpdateCumulativeBins(int binNum, BinUpdateType updateType, 
    ProfileElement &profile, 
    double* dataToAdd);

	// Description:
	// Based upon the additionalQuantityName, returns a double
	// array representing the computation of this quantity. 
	// Always outputs a 3-vector, for quantities such as number in bin
	// which are scalars, result will be in the first position of the vector
	// and equal to the norm of the vector.
	// Currently all of these are calculated
	// from v and from r, which are 3-vectors taken as inputs. Would have to be 
	// rediefined to be more general if other quantities were desired.
	double* CalculateAdditionalProfileQuantity(
    const char *additionalQuantityName,double v[], double r[]);
	
	// Description:
   // After all points have updated the bin statistics, UpdateBinAverages
	// must be called to do the proper averaging and/or postprocessing on 
	// the accumlated columns.
	void BinAveragesAndPostprocessing(vtkPointSet* input, vtkTable* output);

	// Description:
	// merges tableToMerge into originalTable by adding each row,column in
	// tableToMerge to the corresponding row,column in originalTable
	void MergeTables(vtkPointSet* input, vtkTable* originalTable,
		vtkTable* tableToMerge);

	// Description:
	// helper function for MergeTables, to merge two bins by addition, making
	// the original table contain the updated result
	void MergeBins(int binNum, BinUpdateType updateType,
	 	const char *baseName, ColumnType columnType, vtkTable* originalTable,
		vtkTable* tableToMerge);

	// Description:
	// given a base name, a variable index and a column type
	// (TOTAL, AVERAGE,or CUMULATIVE) returns a string representing
	// this data column's name
  std::string GetColumnNameSlow(const char *baseName, ColumnType columnType);

	// Description:
	// Gets a column's data
	vtkVariant GetData(int binNum, const char *baseName,
		ColumnType columnType, vtkTable* output);

  virtual int FillInputPortInformation (int port, vtkInformation *info);
private:
  vtkProfileFilter(const vtkProfileFilter&); // Not implemented
  void operator=(const vtkProfileFilter&); // Not implemented
//ETX
};

#endif


