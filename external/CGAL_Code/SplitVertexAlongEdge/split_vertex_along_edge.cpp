/* =============================================================================================
 *
 * 	split_vertex_along_edge.cpp
 *
 * 	Splits a user supplied set of vertices in a 3D mesh triangulation along
 * 	a set of user supplied edges incident to those vertices in such a way
 * 	that the triangulation property is maintained. Depends on CGAL.
 *
 * 	by Dillon Cislo
 * 	03/05/2019
 *
 * 	This is a MEX-File for MATLAB
 *
 * ============================================================================================*/

#include "mex.h" // for MATLAB

#include <vector>
#include <stdexcept>
#include <iostream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

#include <CGAL/boost/graph/Euler_operations.h>
#include <boost/function_output_iterator.hpp>
#include <boost/foreach.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel 	Kernel;
typedef Kernel::Point_3 					Point;
typedef Kernel::Vector_3 					Vector;
typedef CGAL::Surface_mesh<Point> 				Mesh;

typedef boost::graph_traits<Mesh>::face_descriptor 		Face_index;
typedef boost::graph_traits<Mesh>::vertex_descriptor 		Vertex_index;
typedef boost::graph_traits<Mesh>::halfedge_descriptor		Halfedge_index;
typedef boost::graph_traits<Mesh>::edge_descriptor 		Edge_index;

namespace PMP = CGAL::Polygon_mesh_processing;

///
/// Split a single vertex along a given set of edges
///
void split_single_vertex( Mesh &m, Vertex_index div_ID,
		Vertex_index source_h1_ID, Vertex_index source_h2_ID ) {

	// Find the halfedges along which the vertex will be split
	std::pair< Halfedge_index, bool > h1_pair;
	h1_pair = halfedge( source_h1_ID, div_ID, m );

	Halfedge_index h1;
	if ( h1_pair.second ) { h1 = h1_pair.first;
	} else { std::runtime_error( "Edge not found!" ); }

	// Find the vertices for which the corresponding edges will be removed
	std::vector<Vertex_index> rmEdges;
	BOOST_FOREACH( Vertex_index vv, vertices_around_target( h1, m ) ) {

		if ( vv == source_h2_ID ) { break;
		} else if ( vv != source_h1_ID ) {
			rmEdges.push_back( vv );
		}

	}

	// Remove the appropriate edges
	for( int i = 0; i < rmEdges.size(); i++ ) {

		std::pair< Halfedge_index, bool > hh_pair;
		hh_pair = halfedge( rmEdges[i], div_ID, m );

		Halfedge_index hh;
		if ( hh_pair.second ) { hh = hh_pair.first;
		} else { std::runtime_error( "Edge not found for removal!" ); }

		// Join the face
		Halfedge_index hout = CGAL::Euler::join_face( hh, m );

	}

	// Add the new vertex and edges
	Halfedge_index hnew = CGAL::Euler::add_center_vertex( h1, m );

	// Find the new vertex ID
	Vertex_index vnew = target( hnew, m );

	// Create the new vertex point
	double x = 0.0, y = 0.0, z = 0.0;
	int count = 0;
	BOOST_FOREACH( Vertex_index vv, vertices_around_target( hnew, m ) ) {
		
		Point curPoint = m.point( vv );

		x += curPoint.x();
		y += curPoint.y();
		z += curPoint.z();

		++count;

	}

	x = x / ((double) count);
	y = y / ((double) count);
	z = z / ((double) count);

	m.point( vnew ) = Point( x, y, z );

	// Reposition the old point
	x = 0.0; y = 0.0; z = 0.0;
	count = 0;
	BOOST_FOREACH( Vertex_index vv, vertices_around_target( h1, m ) ) {

		Point curPoint = m.point( vv );

		x += curPoint.x();
		y += curPoint.y();
		z += curPoint.z();

		++count;

	}

	x = x / ((double) count);
	y = y / ((double) count);
	z = z / ((double) count);

	m.point( div_ID ) = Point( x, y, z );

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
		mexErrMsgIdAndTxt( "MATLAB:split_vertex_along_edge:nargin",
				"SPLIT_VERTEX_ALONG_EDGE requires three input arguments." );
	} else if ( nlhs != 3 ) {
		mexErrMsgIdAndTxt( "MATLAB:split_vertex_along_edge:nargout",
				"SPLIT_VERTEX_ALONG_EDGE requires three output arguments." );
	}

	double *face = mxGetPr( prhs[0] );	   // The face connectivity list
	std::size_t numFaces = mxGetM( prhs[0] );  // The number of faces
	std::size_t sizeFaces = mxGetN( prhs[0] ); // The number of vertices in a single face

	// Check that the mesh is a triangulation
	if ( sizeFaces != 3 ) {
		mexErrMsgIdAndTxt( "MATLAB:split_vertex_along_edge:face_degree",
				"Faces must be elements of a triangulation." );
	}

	double *vertex = mxGetPr( prhs[1] ); 	   // The vertex coordinate list
	std::size_t numVertex = mxGetM( prhs[1] ); // The number of vertices
	std::size_t dim = mxGetN( prhs[1] ); 	   // The dimensionality of the vertex list

	// Check the dimensionality of the vertex list
	if ( dim != 3 ) {
		mexErrMsgIdAndTxt( "MATLAB:split_vertex_along_edge:vertex_dim",
				"Vertex coordinates must be 3D." );
	}

	double *divIDx_in = mxGetPr( prhs[2] );
	std::size_t numDivide = mxGetM( prhs[2] );
	std::size_t dimDivide = mxGetN( prhs[2] );

	// Check the dimensionality of the division list
	if ( dimDivide != 3 ) {
		mexErrMsgIdAndTxt( "MATLAB:split_vertex_along_edge:div_dim",
				"Division array is improperly sized." );
	}

	// Format division vertex ID data ------------------------------------------------------
	std::vector<Vertex_index> divIDx;
	divIDx.reserve( numDivide );

	std::vector<Vertex_index> source_h1_IDx;
	source_h1_IDx.reserve( numDivide );

	std::vector<Vertex_index> source_h2_IDx;
	source_h1_IDx.reserve( numDivide );

	for( int i = 0; i < numDivide; i++ ) {

		// NOTE: We subtract 1 from the index to account for
		// MATLAB's 1-indexed array structures
		double div_ID = divIDx_in[i]-1.0;
		double sh1_ID = divIDx_in[i+numDivide]-1.0;
		double sh2_ID = divIDx_in[i+(2*numDivide)]-1.0;

		divIDx.push_back( (Vertex_index) div_ID );
		source_h1_IDx.push_back( (Vertex_index) sh1_ID );
		source_h2_IDx.push_back( (Vertex_index) sh2_ID );

	}

	// Create and populate the polyhedral mesh ---------------------------------------------
	
	// Create vector of 3D point objects
	std::vector<Point> points;
	points.reserve( numVertex );
	for( int i = 0; i < numVertex; i++ ) {

		points.push_back( Point( vertex[i],
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

			// NOTE: We subtract 1 from the index to account for
			// MATLAB's 1-indexed array structures
			double index = face[i+(j*numFaces)]-1.0;
			currentPolygon.push_back( (std::size_t) index );

		}

		polygons.push_back( currentPolygon );

	}

	// Populate the mesh
	Mesh mesh;
	PMP::orient_polygon_soup( points, polygons );
	PMP::polygon_soup_to_polygon_mesh( points, polygons, mesh );

	// -------------------------------------------------------------------------------------
	// MESH PROCESSING
	// -------------------------------------------------------------------------------------
	
	for( int i = 0; i < numDivide; i++ ) {

		split_single_vertex( mesh, divIDx[i],
			       	source_h1_IDx[i], source_h2_IDx[i] );

	}

	// Collect any garbage that may have accumulated in the mesh
	if ( mesh.has_garbage() ) { mesh.collect_garbage(); }

	// Confirm that all faces are triangles
	BOOST_FOREACH( Face_index f, faces( mesh ) ) {
		if ( next(next(halfedge(f,mesh),mesh),mesh)
				!= prev(halfedge(f,mesh),mesh) ) {
			std::cerr << "ERROR: Non-triangular faces present in mesh" << std::endl;

		}
	}

	// -------------------------------------------------------------------------------------
	// OUTPUT PROCESSING
	// -------------------------------------------------------------------------------------
	
	std::size_t numFacesFinal = mesh.number_of_faces(); // Final number of faces
	std::size_t numVertexFinal = mesh.number_of_vertices(); // Final number of vertices
	std::size_t numEdgesFinal = mesh.number_of_edges(); // Final number of edges


	plhs[0] = mxCreateDoubleMatrix( numFacesFinal, sizeFaces, mxREAL );
	double *facesOut = mxGetPr( plhs[0] );

	plhs[1] = mxCreateDoubleMatrix( numVertexFinal, 3, mxREAL );
	double *vertexOut = mxGetPr( plhs[1] );

	plhs[2] = mxCreateDoubleMatrix( numEdgesFinal, 2, mxREAL );
	double *edgesOut = mxGetPr( plhs[2] );


	// Collect face connectivity list
	int i = 0;
	BOOST_FOREACH( Face_index f, mesh.faces() ) {

		// Iterate around the current face to find connectivity
		int j = 0;
		BOOST_FOREACH( Vertex_index v, 
				vertices_around_face( halfedge( f, mesh ), mesh ) ) {

			// NOTE: We add 1 to the index to account for
			// MATLAB's 1-indexed array structures
			facesOut[i+(j*numFacesFinal)] = (double) v+1.0;
			j++;

		}

		i++;
	}

	// Collect vertex coordinates
	i = 0;
	BOOST_FOREACH( Vertex_index v, vertices( mesh ) ) {

		// Vertex coordinates
		Point pp = mesh.point( v );

		vertexOut[i] = pp[0];
		vertexOut[i+numVertexFinal] = pp[1];
		vertexOut[i+(2*numVertexFinal)] = pp[2];

		i++;

	}

	// Collect edge connectivity list
	i = 0;
	BOOST_FOREACH( Edge_index ee, edges( mesh ) ) {

		// NOTE: We add 1 to the index to account for 
		// MATLAB's 1-indexed array structures
		edgesOut[i] = (double) ( source( ee, mesh ) + 1.0 );
		edgesOut[i+numEdgesFinal] = (double) ( target( ee, mesh ) + 1.0 );

		i++;

	}

	return;

};

