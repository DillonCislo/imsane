/* =============================================================================================
 *
 * 	point_set_normals.cpp
 *
 * 	Calculates the oriented unit normal vector field for a disorder 3D point set.
 * 	Depends on the "Point Set Processing" package of CGAL.
 *
 * 	by Dillon Cislo
 * 	02/18/2019
 *
 * 	This is a MEX-file for MATLAB
 *
 * ============================================================================================*/

#include "mex.h" // for MATLAB

#include <vector>
#include <utility>
#include <list>
#include <stdexcept>
#include <iostream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/vcm_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/property_map.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel 	Kernel;
typedef Kernel::Point_3 					Point;
typedef Kernel::Vector_3 					Vector;

typedef std::pair<Point, Vector> 				PointVectorPair;
typedef CGAL::First_of_pair_property_map<PointVectorPair> 	Point_map;
typedef CGAL::Second_of_pair_property_map<PointVectorPair> 	Vector_map;

// Concurrency
#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

///
/// The enumeration of type for different normal estimations
///
enum NORMAL_ESTIMATION {

	// Jet Surface Estimation
	JET_NORMALS = 1,

	// PCA Estimation
	PCA_NORMALS = 2,

	// Voronoi Covariance Measure Estimation
	VCM_NORMALS = 3

};

///
/// Brief main function to call computational functionalities
///
void mexFunction( int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[] ) {

	// -------------------------------------------------------------------------------------
	// INPUT PROCESSING
	// -------------------------------------------------------------------------------------
	
	// Check for proper number of arguments
	if ( nrhs != 3 ) {
		mexErrMsgIdAndTxt( "MATLAB:point_set_normals:nargin",
				"POINT_SET_NORMALS requires three input arguments." );
	} else if ( nlhs != 3 ) {
		mexErrMsgIdAndTxt( "MATLAB:point_set_normals:nargout",
				"POINT_SET_NORMALS requires three output arguments." );
	}

	// The point coordinate list
	double *pts = mxGetPr( prhs[0] );
	std::size_t numPoints = mxGetM( prhs[0] ); // The number of points
	std::size_t dim = mxGetN( prhs[0] ); // The dimensionality of the point list

	// Check the dimensionality of the point list
	if ( dim != 3 ) {
		mexErrMsgIdAndTxt( "MATLAB:point_set_normals:point_dim", 
				"Point coordinates must be 3D." );
	}

	// Estimation procedure processing -----------------------------------------------------
	
	int idx;
	int nn = 0;
	double offset_radius = 0.0;
	double convolution_radius = 0.0;

	// The normal estimation procedure type
	int estimation_procedure;
	if( (idx = mxGetFieldNumber( prhs[1], "estimation_procedure" )) == -1 ) {
		mexErrMsgTxt("No estimation procedure field provided!");
	} else {
		estimation_procedure = (int) *mxGetPr(mxGetFieldByNumber( prhs[1], 0, idx ));
	}

	if ( ( estimation_procedure == 1 ) || ( estimation_procedure == 2 ) ) {

		if( (idx = mxGetFieldNumber( prhs[1], "number_of_neighbors" )) == -1 ) {
			mexErrMsgTxt("Number of point neighbors must be provided!");
		} else {
			nn = (int) *mxGetPr(mxGetFieldByNumber( prhs[1], 0, idx ));
		}

	} else if ( estimation_procedure == 3 ) {

		if( (idx = mxGetFieldNumber( prhs[1], "offset_radius" )) == -1 ) {
			mexErrMsgTxt("No offset radius provided!");
		} else {
			offset_radius = *mxGetPr(mxGetFieldByNumber(prhs[1],0,idx));
		}

		if( (idx = mxGetFieldNumber( prhs[1], "convolution_radius" )) == -1) {
			mexErrMsgTxt("No convolution radius provided!");
		} else {
			convolution_radius = *mxGetPr(mxGetFieldByNumber(prhs[1],0,idx));
		}

	} else {

		mexErrMsgIdAndTxt( "MATLAB:point_set_normals:normal_proc",
				"Invalid normal estimation procedure." );

	}

	// Number of neighbors for normal orientation
	int orient_neighbors = (int) *mxGetPr( prhs[2] );

	// Check number of neighbors
	if ( orient_neighbors < 1 ) {
		mexErrMsgIdAndTxt("MATLAB:point_set_normals:orient_neighbors",
				"Number of neighbors used to orient normals must be positive.");
	}

	// Create Point-Vector Pair vector range -----------------------------------------------
	
	std::vector<PointVectorPair> points;
	points.reserve( numPoints );
	for( int i = 0; i < numPoints; i++ ) {

		Point pp = Point( pts[i], pts[i+numPoints], pts[i+(2*numPoints)] );
		Vector vp = Vector( 0, 0, 0 );

		points.push_back( std::make_pair( pp, vp ) );

	}	

	// -------------------------------------------------------------------------------------
	// ESTIMATE NORMALS
	// -------------------------------------------------------------------------------------
	
	switch( estimation_procedure ) {

		case JET_NORMALS :

			/*
            CGAL::jet_estimate_normals<Concurrency_tag>
				( points, nn,
				  CGAL::parameters::point_map(Point_map()).
				  normal_map(Vector_map()) );
             */
			break;

		case PCA_NORMALS :

			CGAL::pca_estimate_normals<Concurrency_tag>
				( points, nn,
				  CGAL::parameters::point_map(Point_map()).
				  normal_map(Vector_map()) );
			break;

		case VCM_NORMALS :

			CGAL::vcm_estimate_normals( points,
				  offset_radius, convolution_radius,
				  CGAL::parameters::point_map(Point_map()).
				  normal_map(Vector_map()) );

			break;

	}

	// -------------------------------------------------------------------------------------
	// ORIENT NORMAL VECTOR FIELD
	// -------------------------------------------------------------------------------------
	
	std::vector<PointVectorPair>::iterator unoriented_points_begin =
		CGAL::mst_orient_normals( points, orient_neighbors,
				CGAL::parameters::point_map(Point_map()).
				normal_map(Vector_map()) );

	// Erase all points with unoriented normals
	points.erase( unoriented_points_begin, points.end() );

	// -------------------------------------------------------------------------------------
	// OUTPUT PROCESSING
	// -------------------------------------------------------------------------------------
	
	plhs[0] = mxCreateDoubleMatrix( points.size(), 3, mxREAL );
	plhs[1] = mxCreateDoubleMatrix( points.size(), 3, mxREAL );

	double *vn_out = mxGetPr( plhs[0] );
	double *vv_out = mxGetPr( plhs[1] );

	for( int i = 0; i < points.size(); i++ ) {

		vn_out[i] = points[i].second[0];
		vn_out[i+points.size()] = points[i].second[1];
		vn_out[i+(2*points.size())] = points[i].second[2];

		vv_out[i] = points[i].first[0];
		vv_out[i+points.size()] = points[i].first[1];
		vv_out[i+(2*points.size())] = points[i].first[2];

	}

	plhs[2] = mxCreateLogicalMatrix(1,1);
	bool *unoriented = mxGetLogicals( plhs[2] );

	if ( points.size() != numPoints ) {
		*unoriented = true;
	} else {
		*unoriented = false;
	}

	return;

};
