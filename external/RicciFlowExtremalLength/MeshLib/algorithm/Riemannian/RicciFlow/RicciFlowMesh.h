/*!
*      \file RicciFlowMesh.h
*      \brief Mesh for Ricci Flow
*	   \author David Gu
*      \date Documented 10/16/2010
*
*/

#ifndef  _RICCI_FLOW_MESH_H_
#define  _RICCI_FLOW_MESH_H_

#include <map>
#include <vector>

#include "../../../core/Mesh/BaseMesh.h"
#include "../../../core/Mesh/Vertex.h"
#include "../../../core/Mesh/HalfEdge.h"
#include "../../../core/Mesh/Edge.h"
#include "../../../core/Mesh/Face.h"
#include "../../../core/Mesh/iterators.h"
#include "../../../core/Mesh/boundary.h"
#include "../../../core/Parser/parser.h"
#include "../../../core/Parser/traits_io.h"

namespace MeshLib
{

#define TRAIT NORMAL  1
#define TRAIT_FATHER  2
#define TRAIT_UV      4
#define TRAIT_RGB     8
#define TRAIT_PARENT  16

/*!
*	\brief CRicciFlowVertex class
*
*	Vertex class for Ricci flow
*   Traits:
*   index
*	father
*   log_radius
*   curvature
*   vertex uv
*   target_curvature
*   touched
*/
class CRicciFlowVertex : public CVertex
{

public:

	/*! 
	 *	Each bit in the traits indicates whether the vertex class has the correpsonding trait.
	 *  e.g. if ( traits & TRAIT_UV ), then vertex->huv() needs to be stored int the vertex string.
	 */
	static unsigned int traits;

	/*!
	 *	CRicciFlowVertex constructor
	 */
	CRicciFlowVertex() { m_index = 0; };
	/*!
	 *	CRicciFlowVertex destructor
	 */
	~CRicciFlowVertex() {};

	/*!
	 *	Vertex uv trait
	 */
	CPoint2 & huv() { return m_huv; };
	/*!
	 *	Vertex index trait
	 */
	int &     idx() { return m_index; };
	/*!
	 *  Vertex father trait
	 */
	int &	  father(){ return m_father;};
	/*!
	 *	Vertex log radius
	 */
	double &   u()    { return m_log_radius; };
	/*!
	 *	Vertex curvature trait
	 */
	double &   k()    { return m_curvature; };
	/*!
	 *	Vertex target curvature
	 */
	double & target_k() { return m_target_curvature; };
	/*!
	 *	whether the vertex has been touched
	 */
	bool  & touched()  { return m_touched; };

	/*!
	 *	Read vertex uv from vertex string
	 */
	void _from_string();
	/*!
	 *	Write vertex uv to vertex string
	 */
	void _to_string();
	/*!
	 * Topological valence of the vertex
	 */
	int & valence() { return m_valence; };
	/*!
	 *	Vertex color
	 */
	//CPoint & rgb() { return m_rgb; };

protected:
	/*! Vertex uv */
	CPoint2 m_huv;
	/*! Vertex index */
	int     m_index;
	//father
	int     m_father;
	//log radius
	double  m_log_radius;
	//curvature
	double  m_curvature;
	//target curvature
	double  m_target_curvature;
	//touched
	bool    m_touched;
	//topological valence
	int     m_valence;
	//vertex color
	//CPoint  m_rgb;
};

//read father from string
inline void CRicciFlowVertex::_from_string()
{
  CParser parser( m_string );

  for( std::list<CToken*>::iterator iter = parser.tokens().begin() ; iter != parser.tokens().end(); ++ iter )
  {
	  CToken * token = *iter;
	  if( token->m_key == "father" )
	  {
		 std::string line = strutil::trim( token->m_value, "()");
		 m_father = strutil::parseString<int>( line );	
		 continue;
		
	  }
	  /*
	  if( token->m_key == "rgb" )
	  {
		 std::string line = strutil::trim( token->m_value, "()");
		 line >> m_rgb;
		 continue;
	  }
	  */
	  /*
	  if( token->m_key == "uv" )
	  {
		 std::string line = strutil::trim( token->m_value, "()");
		 line >> m_huv;		
		 continue;
	  }
	  */
  }
}

//converting vertex uv trait to string

inline void CRicciFlowVertex::_to_string()
{
	CParser parser( m_string );
	
	if( traits & TRAIT_UV )
	{
		parser._removeToken( "uv" );
	}
	//if( traits & TRAIT_RGB )
	//{
	//	parser._removeToken( "rgb" );
	//}
	parser._toString( m_string );
	std::stringstream iss;
	
	if( traits & TRAIT_UV )
	{
		iss << "uv=(" << m_huv[0] << " " << m_huv[1] << ") ";
	}
	//if( traits & TRAIT_RGB )
	//{
	//	iss << "rgb=(" << m_rgb[0] << " " << m_rgb[1] << " " << m_rgb[2] << ") ";
	//}
	if( m_string.length() > 0 )
	{
		m_string += " ";
	}
	m_string += iss.str();
}

/*!
*	\brief CRicciFlowEdge class
*
*	Edge class for Ricci flow
*   trait:
*   edge length
*   edge weight \$f \frac{\partial \theta_i}{\partial u_j} = w_k \f$
*   edge inversive distance \f$ \cos \phi \f$, intersection angle is \f$\phi\f$
*/
class CRicciFlowEdge : public  CEdge
  {
  public:
    /*!	CRicciFlowEdge constructor
	 */
	 CRicciFlowEdge() { m_length = 0; m_weight = 0; m_sharp = false; };
    /*!	CRicciFlowEdge destructor
	 */
    ~CRicciFlowEdge(){};
	/*!	edge weight trait
	 */
	double & weight() { return m_weight; };
	/*!	edge length trait
	 */
	double & length() { return m_length; };
	/*! edge inversive distance
	*/
	double & inversive_distance() { return m_inversive_distance; };
	/*! edge sharp
	 */
	bool sharp() { return m_sharp; };
    /*!
	 *	read sharp trait from the edge string
	 */
	void _from_string() 
	{ 
		CEdge::_from_string(); 
		CParser parser( m_string );

		  for( std::list<CToken*>::iterator iter = parser.tokens().begin() ; iter != parser.tokens().end(); ++ iter )
		  {
			  CToken * token = *iter;
			  if( token->m_key == "sharp" )
			  {
				  m_sharp = true;
				  break;
			  }
		  }	
	};
	/*!
	 *	write no traits to edge string
	 */
	void _to_string() { CEdge::_to_string(); };

  protected:
	  /*! edge weight trait */
	double   m_weight;
	  //edge length
	double   m_length;
	//inversive distance
	double   m_inversive_distance;
	/*!
	 *	if edge is sharp
	 */
	bool     m_sharp;
};

/*!
*	\brief CRicciFlowHalfEdge class
*
*	HalfEdge class for Ricci flow
*   trait:
*   corner angle
*/

class CRicciFlowHalfEdge : public  CHalfEdge
{
public:
  CRicciFlowHalfEdge() { m_angle = 0.0; };
  ~CRicciFlowHalfEdge() {};
  /*! corner angle */
  double & angle() { return m_angle; };
  /*! on half edge [v_i,v_j], \f$ \frac{\parital \theat_i}{\partial u_j} */		
  double & dtheta_du() { return m_theta_u; };
  /*! on half edge [v_i,v_j], \f$ \frac{\parital l_{ij}}{\partial u_j} */		
  double & dl_du()     { return m_l_u; };

protected:
  double m_angle;
  CPoint m_s;
  double m_theta_u;
  double m_l_u;

};

/*!
*	\brief CRicciFlowFace class
*
*	Face class for Ricci flow
*   trait:
*	whether the face has been touched
*/

class CRicciFlowFace: public CFace
{
public:
	CRicciFlowFace() { m_touched = false; };
	~CRicciFlowFace() {};
	bool & touched() { return m_touched; };
	CPoint & normal(){ return m_normal; };
protected:
	bool m_touched;
	CPoint m_normal;
};

/*-------------------------------------------------------------------------------------------------------------------------------------

	Ricci flow Mesh Class

--------------------------------------------------------------------------------------------------------------------------------------*/
/*!
 *	\brief CHarmonicMapperMesh class
 *
 *	Mesh class for harmonic mapping purpose
 */
//template<typename V, typename E, typename F, typename H>
template< class V, class E, class F, class H > // Added by Dillon 2017/08/22
class CRicciFlowMesh : public CBaseMesh<V,E,F,H>
{
public:
	typedef MeshLib::CBoundary<V,E,F,H> CBoundary;
	typedef MeshLib::CLoop<V,E,F,H> CLoop;
	typedef MeshLib::MeshVertexIterator<V,E,F,H> MeshVertexIterator;
	typedef MeshLib::MeshEdgeIterator<V,E,F,H> MeshEdgeIterator;
	typedef MeshLib::VertexVertexIterator<V,E,F,H> VertexVertexIterator;
	typedef MeshLib::VertexEdgeIterator<V,E,F,H> VertexEdgeIterator;
	typedef MeshLib::VertexFaceIterator<V,E,F,H> VertexFaceIterator;
	typedef MeshLib::VertexInHalfedgeIterator<V,E,F,H> VertexInHalfedgeIterator;
	typedef MeshLib::VertexOutHalfedgeIterator<V,E,F,H> VertexOutHalfedgeIterator;
	typedef MeshLib::FaceHalfedgeIterator<V,E,F,H> FaceHalfedgeIterator;
	typedef MeshLib::FaceEdgeIterator<V,E,F,H> FaceEdgeIterator;
	typedef MeshLib::MeshFaceIterator<V,E,F,H> MeshFaceIterator;
	typedef MeshLib::FaceVertexIterator<V,E,F,H> FaceVertexIterator;
	typedef MeshLib::MeshHalfEdgeIterator<V,E,F,H> MeshHalfEdgeIterator;
public:
	static unsigned long long m_input_traits;
	static unsigned long long m_output_traits;
};

template< class V, class E, class F, class H >
unsigned long long CRicciFlowMesh<V,E,F,H>::m_input_traits = CBaseMesh<V,E,F,H>::m_input_traits;

template< class V, class E, class F, class H >
unsigned long long CRicciFlowMesh<V,E,F,H>::m_output_traits = CBaseMesh<V,E,F,H>::m_output_traits;

typedef CRicciFlowMesh<CRicciFlowVertex, CRicciFlowEdge, CRicciFlowFace, CRicciFlowHalfEdge> CRFMesh;

template<> unsigned long long CRFMesh::m_input_traits = EDGE_SHARP| VERTEX_FATHER | VERTEX_RGB;
template<> unsigned long long CRFMesh::m_output_traits = VERTEX_UV;
};
#endif  //_RICCI_FLOW_MESH_H_
