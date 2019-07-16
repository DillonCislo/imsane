/*!
*      \file main.cpp
*      \brief Main program for Tangential Ricci Flow for Extremal Length
*	   \author David Gu
*      \date Document 10/07/2012
*
*		The code performs the following tasks
*
*/

#include <iostream>

#include "../MeshLib/algorithm/Structure/Structure.h"
//for Extremal Lengths
#include "../MeshLib/algorithm/Riemannian/RicciFlow/TangentialRicciExtremalLength.h"
#include "../MeshLib/algorithm/Riemannian/RicciFlow/EuclideanEmbed.h"
//for Disk Maps
#include "../MeshLib/algorithm/Riemannian/RicciFlow/TangentialRicciFlow.h"

using namespace MeshLib;

unsigned int CRicciFlowVertex::traits = 0;

/******************************************************************************************************************************
*
*	Extremal Length
*
*******************************************************************************************************************************/

//-tangent_ricci_extremal_length sophie.remesh.m sophie.uv.m
void _tangent_ricci_extremal_length( const char * _input_mesh, const char * _mesh_with_uv )
{
	CRicciFlowVertex::traits = CRicciFlowVertex::traits | TRAIT_UV;

	CRFMesh mesh;
	mesh.read_m( _input_mesh );
	//mesh.read_obj( _input_mesh );

	CTangentialRicciFlowExtremalLength<CRicciFlowVertex,CRicciFlowEdge,CRicciFlowFace,CRicciFlowHalfEdge> mapper(&mesh);
	mapper._calculate_metric();

	CRFEmbed embed( &mesh );
	embed._embed();
	mesh.write_m( _mesh_with_uv );
}

// -tangent_ricci
void _tangent_ricci( const char * _input_mesh, const char * _mesh_with_uv )
{
	CRicciFlowVertex::traits = CRicciFlowVertex::traits | TRAIT_UV;

	CRFMesh mesh;
	mesh.read_m( _input_mesh );

	CTangentialRicciFlow<CRicciFlowVertex,CRicciFlowEdge,CRicciFlowFace,CRicciFlowHalfEdge> mapper(&mesh);
	mapper._calculate_metric();

	std::cout << "Metric has been calculated\n";

	CRFEmbed embed( &mesh );
	embed._embed();
	mesh.write_m( _mesh_with_uv );
}


/*!	\brief main function to call all the functionalities
 * 
 */

int main( int argc, char * argv[] )
{

/*------------------------------------------------------------------------------

	Extremal Length

-------------------------------------------------------------------------------*/
std::cout << "You have entered " << argc
	<< " arguments:" << "\n";

//for (int i = 0; i < argc; ++i){
//	std::cout << argv[i] << "\n";
//}

//-tangent_ricci_extremal_length sophie.remesh.m sophie.uv.m
if( strcmp( argv[1] , "-tangent_ricci_extremal_length") == 0 && argc == 4)
{
	_tangent_ricci_extremal_length( argv[2], argv[3]);
	return 0;
}

if( strcmp( argv[1] , "-tangent_ricci") == 0 && argc == 4)
{
	_tangent_ricci( argv[2], argv[3]);
	return 0;
}

	return 0;
}
 