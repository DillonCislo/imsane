/* ==============================================================================================
 *
 * 	surface_geodesic_pairs.cpp
 *
 * 	This function calculates the surface geodesics between sets of point pairs on the surface
 * 	of a 3D triangulation.  The points are assumed to elements of face interiors, 
 * 	i.e. not vertices or elements of edges.
 *
 * 	by Dillon Cislo
 * 	01/30/2019
 *
 * 	This is a MEX-file for MATLAB
 *
 * ============================================================================================*/

#include "mex.h" // for MATLAB

#include <cstdlib>
#include <iostream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/Surface_mesh_shortest_path/function_objects.h>
#include <CGAL/Surface_mesh_shortest_path/barycentric.h>
#include <CGAL/Surface_mesh_shortest_path/internal/misc_functions.h>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>

#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

#include <boost/lexical_cast.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel 		Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3>				Triangle_mesh;

typedef boost::property_map<Triangle_mesh, CGAL::vertex_point_t>::type 	VertexPointMap;
typedef boost::graph_traits<Triangle_mesh> 				Graph_traits;
typedef Graph_traits:: vertex_descriptor 				vertex_descriptor;
typedef Graph_traits::vertex_iterator 					vertex_iterator;
typedef Graph_traits::face_descriptor 					face_descriptor;
typedef Graph_traits::face_iterator 					face_iterator;

typedef CGAL::Surface_mesh_shortest_path_traits<Kernel, Triangle_mesh> 	Traits;
typedef CGAL::Surface_mesh_shortest_path<Traits> 	Surface_mesh_shortest_path;
typedef Surface_mesh_shortest_path::Face_location 			Face_location;
typedef Traits::Point_3 						Point_3;

typedef CGAL::AABB_face_graph_triangle_primitive<Triangle_mesh, VertexPointMap> AABB_primitive;
typedef CGAL::AABB_traits<Kernel, AABB_primitive> 				AABB_traits;
typedef CGAL::AABB_tree<AABB_traits> 						AABB_tree;

///
/// Calculates the surface geodesics between pairs of points on a mesh triangulation
///
std::vector< std::vector<Point_3> > calculate_geodesic_pairs(
		Surface_mesh_shortest_path &shortest_paths,
		const std::vector< std::pair<std::size_t, std::size_t> > bondIDx,
		const std::vector<Face_location> cell_locations ) {

	// Construct the output vector
	std::vector< std::vector<Point_3> > geodesic_pairs;
	geodesic_pairs.reserve( bondIDx.size() );

	// Iterate over the bonds to construct the geodesic paths ------------------------------
	
	// Remove any existing source points
	if ( shortest_paths.number_of_source_points() != 0 ) {
		shortest_paths.remove_all_source_points();
	}

	// Add the source point of the first bond
	shortest_paths.add_source_point( cell_locations[ bondIDx[0].first ] );

	for( int i = 0; i < bondIDx.size(); i++ ) {

		// Update the current source if it is not an element of the current bond
		Face_location currentSourceLocation = *shortest_paths.source_points_begin();

		if( currentSourceLocation != cell_locations[ bondIDx[i].first ] ) {

			shortest_paths.remove_all_source_points();
			shortest_paths.add_source_point( cell_locations[ bondIDx[i].first ] );

		}

		std::vector<Point_3> currentBond;
		shortest_paths.shortest_path_points_to_source_points(
				cell_locations[ bondIDx[i].second ].first,
				cell_locations[ bondIDx[i].second ].second,
				std::back_inserter( currentBond ) );

		geodesic_pairs.push_back( currentBond );

	}

	return geodesic_pairs;

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
	if ( nrhs != 4 ) {
		mexErrMsgIdAndTxt( "MATLAB:surface_geodesic_pairs:nargin",
				"SURFACE_GEODESIC_PAIRS requires five input arguments." );
	} else if ( nlhs != 3 ) {
		mexErrMsgIdAndTxt( "MATLAB:surface_geodesic_pairs:nargout",
				"SURFACE_GEODESIC_PAIRS requires one output argument." );
	}

	// The face connectivity list
	std::size_t *faces = ( std::size_t* ) mxGetData( prhs[0] );
	std::size_t numFaces = mxGetM( prhs[0] );   // The number of faces
	std::size_t sizeFaces = mxGetN( prhs[0] );  // The number of vertices in a single face

	// Check the the mesh is a triangulation
	if ( sizeFaces != 3 ) {
		mexErrMsgIdAndTxt( "MATLAB:surface_geodesic_pairs:face_size",
				"Faces must be elements of a trianuglation." );
	}

	// The vertex coordinate list
	double *vertex = mxGetPr( prhs[1] );
	std::size_t numVertex = mxGetM( prhs[1] ); // The number of vertices
	std::size_t dim = mxGetN( prhs[1] ); 	   // The dimensionality of the vertex list

	// Check the dimensionality of the vertex list
	if ( dim != 3 ) {
		mexErrMsgIdAndTxt( "MATLAB:surface_geodesic_pairs:vertex_dim",
				"Vertex coordinates must be 3D." );
	}

	// The bond ID pair list
	std::size_t *bondID = ( std::size_t* ) mxGetData( prhs[2] );
	std::size_t numBonds = mxGetM( prhs[2] ); 	// The number of unique bonds
	std::size_t sizeBonds = mxGetN( prhs[2] );	// The number of vertices in a bond

	// Check bond size
	if ( sizeBonds != 2 ) {
		mexErrMsgIdAndTxt( "MATLAB:surface_geodesic_pairs:bond_dim",
				"Bonds must defined in terms of two cells only." );
	}

	// The 3D coordinates of the cell centroids
	double *cellCoords = mxGetPr( prhs[3] );
	std::size_t numCells = mxGetM( prhs[3] );  // The number of cells
	std::size_t cellDim = mxGetN( prhs[3] );   // The dimensionality of cell coordinates

	// Check the dimensionality of the cell coordinates
	if ( cellDim != 3 ) {
		mexErrMsgIdAndTxt( "MATLAB:surface_geodesic_pairs:cell_dim",
				"Cell coordinates must be 3D." );
	}

	// Re-format the bond ID pair list -----------------------------------------------------
	std::vector< std::pair<std::size_t, std::size_t> > bondIDx;
	bondIDx.reserve( numBonds );
	for( int i = 0; i < numBonds; i++ ) {

		bondIDx.push_back( std::make_pair( bondID[i], bondID[i+numBonds] ) );

	}

	// Create and populate the polyhedral mesh ---------------------------------------------
	
	// Create vector of 3D point objects
	std::vector<Kernel::Point_3> points;
	points.reserve( numVertex );
	for( int i = 0; i < numVertex; i++ ) {

		points.push_back( Kernel::Point_3( vertex[i],
				  	 	   vertex[i+numVertex],
						   vertex[i+(2*numVertex)] ) );

	}

	// Create vector of polygon objects
	std::vector< std::vector<std::size_t> > polygons;
	polygons.reserve( numFaces );
	for( int i = 0; i < numFaces; i++ ) {

		std::vector<std::size_t> currentPolygon;
		currentPolygon.reserve( sizeFaces );
		for( int j = 0; j < sizeFaces; j++ ) {

			currentPolygon.push_back( faces[i+(j*numFaces)] );

		}

		polygons.push_back( currentPolygon );

	}

	// Populate the mesh
	Triangle_mesh tmesh;
	CGAL::Polygon_mesh_processing::orient_polygon_soup( points, polygons );
	CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh( points, polygons, tmesh );
	
	// Create the Surface_mesh_shortest_path object ----------------------------------------
	Surface_mesh_shortest_path shortest_paths( tmesh );

	// Find the face locations of each cell centroid ---------------------------------------
	
	// Cache the AABB tree object
	AABB_tree tree;
	shortest_paths.build_aabb_tree( tree );

	// Construct a vector of face locations
	std::vector<Face_location> cell_locations;
	cell_locations.reserve( numCells );
	for( int i = 0; i < numCells; i++ ) {

		Point_3 cellPoint = Point_3( cellCoords[i],
					     cellCoords[i+numCells],
				             cellCoords[i+(2*numCells)] );

		cell_locations.push_back( shortest_paths.locate( cellPoint, tree ) );

	}

	// -------------------------------------------------------------------------------------
	// CALCULATE BOND GEODESICS
	// -------------------------------------------------------------------------------------
	
	std::vector< std::vector<Point_3> > geodesic_pairs;
	geodesic_pairs = calculate_geodesic_pairs( shortest_paths, bondIDx, cell_locations );

	// -------------------------------------------------------------------------------------
	// OUTPUT PROCESSING
	// -------------------------------------------------------------------------------------
	
	// Create geodesic_path output cell array ----------------------------------------------
	plhs[0] = mxCreateCellMatrix( numBonds, 1 );

	for( int i = 0; i < numBonds; i++ ) {

		std::vector<Point_3> currentGeodesic = geodesic_pairs[i];
		int numPoints = currentGeodesic.size();

		mxArray *pointSequenceOut = mxCreateDoubleMatrix( numPoints, 3, mxREAL );
		double *pointSequence = mxGetPr( pointSequenceOut );
		for( int j = 0; j < numPoints; j++ ) {

			pointSequence[j] = currentGeodesic[j].x();
			pointSequence[j+numPoints] = currentGeodesic[j].y();
			pointSequence[j+(2*numPoints)] = currentGeodesic[j].z();

		}

		mxSetCell( plhs[0], i, pointSequenceOut );

	}

	// Create barycentric face coordinate output for each cell centroid --------------------
	
	plhs[1] = mxCreateDoubleMatrix( numCells, 1, mxREAL );
	double *cellIndexOut = mxGetPr( plhs[1] );

	plhs[2] = mxCreateDoubleMatrix( numCells, 3, mxREAL );
	double *cellBaryOut = mxGetPr( plhs[2] );

	for( int i = 0; i < numCells; i++ ) {

		Face_location currentCell = cell_locations[i];

		cellIndexOut[i] = (double) currentCell.first;

		cellBaryOut[i] = currentCell.second[0];
		cellBaryOut[i+numCells] = currentCell.second[1];		
		cellBaryOut[i+(2*numCells)] = currentCell.second[2];

	}

	return;

};
