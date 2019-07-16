/*=============================================================
 * ricci_flow_extremal_length.cpp
 *
 * A test of employing the RFEL pipeline written by David Gu
 * as a mex function
 * 
 * The calling syntax is
 * 	ricci_flow_extremal_length.cpp
 *
 * This is a MEX-file for MATLAB
 * ==========================================================*/

#include "mex.h" // for MATLAB

#include "../MeshLib/algorithm/Structure/Structure.h"
// for Extremal Lengths
#include "../MeshLib/algorithm/Riemannian/RicciFlow/TangentialRicciExtremalLength.h"
#include "../MeshLib/algorithm/Riemannian/RicciFlow/EuclideanEmbed.h"

using namespace MeshLib;

void _main();

unsigned int CRicciFlowVertex::traits = 0;

/**********************************************************
 *
 * 	Extremal Length
 *
 *********************************************************/

void _tangent_ricci_extremal_length( const char * _input_mesh, const char * _mesh_with_uv ) {
	CRicciFlowVertex::traits = CRicciFlowVertex::traits | TRAIT_UV;

	CRFMesh mesh;
	mesh.read_m( _input_mesh );

	CTangentialRicciFlowExtremalLength<CRicciFlowVertex,CRicciFlowEdge,CRicciFlowFace,CRicciFlowHalfEdge> mapper(&mesh);
	mapper._calculate_metric();

	CRFEmbed embed( &mesh );
	embed._embed();
	mesh.write_m( _mesh_with_uv );

}

// Brief main function to call functionalities

void mexFunction( int nlhs, mxArray *plhs[], 
        int nrhs, const mxArray *prhs[])
{

	// Variable definitions
	int i;
	char *fCommand, *meshFileIn, *meshFileOut;

	// Check for proper number of arguments
	if (nrhs != 3){
		mexErrMsgIdAndTxt("MATLAB:ricci_flow_extremal_length:nargin",
				"RICCI_FLOW_EXTREMAL_LENGTH requires three input arguments.");
	} else if (nlhs >= 1){
		mexErrMsgIdAndTxt("MATLAB:ricci_flow_extremal_length:nargout",
				"RICCI_FLOW_EXTREMAL_LENGTH requires no output arguments.");
	}

	// Make sure that the arguments are strings (char arrays)
	for(i=0; i<nrhs; i++){
		if(!mxIsChar(prhs[i]) || mxIsComplex(prhs[i])){
			mexErrMsgIdAndTxt("MATLAB:ricci_flow_extremal_length:inputNotString",
					"Inputs must be strings.");
		}
	}
	
	// Get the flow command and file names
	fCommand = mxArrayToString(prhs[0]);
	meshFileIn = mxArrayToString(prhs[1]);
	meshFileOut = mxArrayToString(prhs[2]);

	if( strcmp(fCommand, "-tangent_ricci_extremal_length") == 0 ) {
		_tangent_ricci_extremal_length(meshFileIn, meshFileOut);
		return;
	}

	return;
	
}
