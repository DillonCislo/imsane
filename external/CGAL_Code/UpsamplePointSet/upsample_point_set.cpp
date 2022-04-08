/* =============================================================================================
 *
 *  upsample_point_set.cpp
 *
 *  This function uses an edge aware method to upsample an disordered point set.
 *  Depends on CGAL.
 *
 *  by Dillon Cislo
 *  12/09/2021
 *
 *  This is a MEX-file for MATLAB
 *  
 * ============================================================================================*/

#include "mex.h" // for MATLAB

#include <vector>
#include <fstream>
#include <iostream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/edge_aware_upsample_point_set.h>

typedef CGAL::Simple_cartesian<double>              Kernel;
typedef Kernel::Point_3                             Point;
typedef Kernel::Vector_3                            Vector;
typedef std::pair<Point, Vector> 				            Pwn;
typedef CGAL::First_of_pair_property_map<Pwn>		 	  Point_map;
typedef CGAL::Second_of_pair_property_map<Pwn>		 	Vector_map;

// Concurrency
#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

// Main function
void mexFunction( int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[] ) {

	// -------------------------------------------------------------------------------------
	// INPUT PROCESSING
	// -------------------------------------------------------------------------------------

  // Check for proper number of arguments
  if ( nrhs != 6 ) {
    mexErrMsgIdAndTxt( "MATLAB:upsample_point_set:nargin",
        "UPSAMPLE_POINT_SET requires 6 input arguments" );
  } else if ( nlhs != 2 ) {
    mexErrMsgIdAndTxt("MATLAB:upsample_point_set:nargout",
        "UPSAMPLE_POINT_SET requries 2 output arguments" );
  }

	// The input point cloud list
	double *pts = mxGetPr( prhs[0] );
	std::size_t numPoints = mxGetM( prhs[0] ); // The number of points
	std::size_t dim = mxGetN( prhs[0] ); // The dimensionality of the vertex list

	// Check the dimensionality of the point cloud
	if ( dim != 3 ) {
		mexErrMsgIdAndTxt( "MATLAB:upsample_point_set:point_dim",
				"Point coordinates must be 3D" );
	}

	// The input point normal list
	double *pn = mxGetPr( prhs[1] );
	std::size_t numNormals = mxGetM( prhs[1] );
	std::size_t dimNormals = mxGetN( prhs[1] );

	// Check the dimensionality of the point normals
	if ( dimNormals != 3 ) {
		mexErrMsgIdAndTxt( "MATLAB:upsample_point_set:normal_dim",
				"Normal vectors must be 3D" );
	}

	// Check that the number of elements in each list is consistent
	if ( numPoints != numNormals ) {
		mexErrMsgIdAndTxt( "MATLAB:upsample_point_set:invalid_point_cloud",
				"Point cloud is inconsistently sized." );
	}

  // The number of output points
  int number_of_output_points = (int) *mxGetPr( prhs[2] );

  // Control the sharpness of the result
  double sharpness_angle = *mxGetPr( prhs[3] );

  // Controls sensitivity to edges
  // Higher values will sample more points near the edges
  double edge_sensitivity = *mxGetPr( prhs[4] );

  // Initial size of neighborhood
  double neighbor_radius = *mxGetPr( prhs[5] );

	// Load points and normal vectors into pair structure ----------------------------------
	
	std::vector<Pwn> points;
	points.reserve( numPoints );
	for( int i = 0; i < numPoints; i++ ) {

		Point pp = Point( pts[i], pts[i+numPoints], pts[i+(2*numPoints)] );
		Vector vp = Vector( pn[i], pn[i+numPoints], pn[i+(2*numPoints)] );

		points.push_back( std::make_pair( pp, vp ) );

	}

	// -------------------------------------------------------------------------------------
	// UPSAMPLE POINT SET
	// -------------------------------------------------------------------------------------
  
  CGAL::edge_aware_upsample_point_set<Concurrency_tag>(
      points.begin(), points.end(),
      std::back_inserter( points ), Point_map(), Vector_map(),
      sharpness_angle, edge_sensitivity, neighbor_radius,
      number_of_output_points );

  // ------------------------------------------------------------------------------------
  // OUTPUT PROCESSING
  // ------------------------------------------------------------------------------------
  
  plhs[0] = mxCreateDoubleMatrix( points.size(), dim, mxREAL );
  plhs[1] = mxCreateDoubleMatrix( points.size(), dim, mxREAL );

  double *vv_out = mxGetPr( plhs[0] );
  double *vn_out = mxGetPr( plhs[1] );

  for ( int i = 0; i < points.size(); i++ ) {
    for ( int j = 0; j < dim; j++ ) {

      vv_out[i+(j*points.size())] = points[i].first[j];
      vn_out[i+(j*points.size())] = points[i].second[j];

    }
  }

  return;

};
