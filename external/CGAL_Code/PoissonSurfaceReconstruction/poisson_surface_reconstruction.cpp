/* =============================================================================================
 *
 * 	poisson_surface_reconstruction.cpp
 *
 * 	This function creates a mesh triangulation from a disordered 3D point cloud
 * 	with oriented vertex normals. Depends on CGAL.
 *
 * 	by Dillon Cislo
 * 	02/18/2019
 *
 * 	This is a MEX-file for MATLAB
 *
 * ============================================================================================*/

#include "mex.h" // for MATLAB

#include <vector>
#include <iostream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/poisson_surface_reconstruction.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel 	Kernel;
typedef Kernel::Point_3 					Point;
typedef Kernel::Vector_3 					Vector;

typedef std::pair<Point, Vector> 				Pwn;
typedef CGAL::First_of_pair_property_map<Pwn>		 	Point_map;
typedef CGAL::Second_of_pair_property_map<Pwn>		 	Vector_map;

typedef CGAL::Polyhedron_3<Kernel> 				Polyhedron;

//Concurrency
#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

///
/// Brief main function to call computational functionalities
///
void mexFunction( int nlhs, mxArray *plhs[],
		int nrhs, const mxArray *prhs[] ) {

	// -------------------------------------------------------------------------------------
	// INPUT PROCESSING
	// -------------------------------------------------------------------------------------
	
	// Check for proper number of arguments
	if ( nrhs != 2 ) {
		mexErrMsgIdAndText( "MATLAB:poisson_surface_reconstruction:nargin",
				"POISSON_SURFACE_RECONSTRUCTION requires two input arguments" );
	} else if ( nlhs != 2 ) {
		mexErrMsgIdAndText( "MATLAB:poisson_surface_reconstruction:nargout",
				"POISSON_SURFACE_RECONSTRUCTION requires two output arguments" );
	}

	// The input point cloud list
	double *pts = mxGetPr( prhs[0] );
	std::size_t numPoints = mxGetM( prhs[0] ); // The number of points
	std::size_t dim = mxGetN( prhs[0] ); // The dimensionality of the vertex list

	// Check the dimensionality of the point cloud
	if ( dim != 3 ) {
		mexErrMsgIdAndText( "MATLAB:poisson_surface_reconstruction:point_dim",
				"Point coordinates must be 3D" );
	}

	// The input point normal list
	double *pn = mxGetPr( prhs[1] );
	std::size_t numNormals = mxGetM( prhs[1] );
	std::size_t dimNormals = mxGetN( prhs[1] );

	// Check the dimensionality of the point normals
	if ( dimNormals != 3 ) {
		mexErrMsgIdAndText( "MATLAB:poisson_surface_reconstruction:normal_dim",
				"Normal vectors must be 3D" );
	}

	// Check that the number of elements in each list is consistent
	if ( numPoints != numNormals ) {
		mexErrMsgIdAndText( "MATLAB:poisson_surface_reconstruction:invalid_point_cloud",
				"Point cloud is inconsistently sized." );
	}

	// Load points and normal vectors into pair structure ----------------------------------
	
	std::vector<Pwn> points;
	points.reserve( numPoints );
	for( int i = 0; i < numPoints; i++ ) {

		Point pp = Point( pts[i], pts[i+numPoints], pts[i+(2*numPoints)] );
		Vector vp = Vector( pn[i], pn[i+numPoints], pn[i+(2*numPoints)] );

		points.push_back( std::make_pair( pp, vp ) );

	}

	// -------------------------------------------------------------------------------------
	// CONSTRUCT MESH TRIANGULATION
	// -------------------------------------------------------------------------------------
	
	Polyhedron output_mesh;

	double average_spacing = CGAL::compute_average_spacing<Concurrency_tag>
		( points, 6, CGAL::parameters::point_map(Point_map()) );

	bool mesh_success = CGAL::poisson_surface_reconstruction_delaunay
		( points.begin(), points.end(), Point_map(), Vector_map(),
		  output_mesh, average_spacing );

	if (!mesh_success) {
		mesErrMsgTxt( "Mesh could not be constructed properly" );
	}

	// -------------------------------------------------------------------------------------
	// OUTPUT PROCESSING
