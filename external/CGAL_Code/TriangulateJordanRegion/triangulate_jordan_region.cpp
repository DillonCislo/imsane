/* =============================================================================================
 *
 *  triangulate_jordan_region.cpp
 *
 *  Triangulates a 2D Jordan region defined by an ordered set of input vertices
 *
 *  by Dillon Cislo
 *  05/02/2019
 *
 *  This is a MEX-file for MATLAB
 *
 * ============================================================================================*/

#define CGAL_MESH_2_OPTIMIZER_VERBOSE

#include "mex.h" // for MATLAB

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

#include <CGAL/lloyd_optimize_mesh_2.h>

#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Delaunay_mesh_vertex_base_2<K>                  Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K>                    Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>          Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds>    CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT>              Criteria;
typedef CGAL::Delaunay_mesher_2<CDT, Criteria>                Mesher;

typedef CDT::Vertex_handle  Vertex_handle;
typedef CDT::Point          Point;

///
/// Brief main function to call computational functionalities
///
void mexFunction( int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[] ) {

  // -------------------------------------------------------------------------------------------
  // INPUT PROCESSING
  // -------------------------------------------------------------------------------------------
  
  // Check for proper number of arguments
  if ( nrhs != 1 ) {
    mexErrMsgIdAndTxt("MATLAB:triangulate_jordan_region:nargin",
        "TRIANGULATE_JORDAN_REGION requires one input argument.");
  } else if ( nlhs != 2 ) {
    mexErrMsgIdAndTxt("MATLAB:triangulate_jordan_region:nargout",
        "TRIANGULATE_JORDAN_REGION requires two output arguments.");
  }

  double *vjr = mxGetPr( prhs[0] ); // The original boundary vertex coordinate list
  int numVJR = (int) mxGetM( prhs[0] ); // The number of vertices
  int dim = (int) mxGetN( prhs[0] ); // The dimensions of the vertex coordinates

  // Check dimensionality of the vertex list
  if ( dim != 3) {
    mexErrMsgIdAndTxt("MATLAB:triangulate_jordan_region:dimension",
        "Vertices must be 2D.");
  }


