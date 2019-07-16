/*! \file TangentialRicciFlow.h
 *  \brief General Euclidean Ricci flow algorithm
 *  \author David Gu
 *  \date   documented on 10/17/2010
 *
 *	Algorithm for general Ricci Flow
 */

#ifndef _TANGENTIAL_RICCI_FLOW_H_
#define _TANGENTIAL_RICCI_FLOW_H_

#include <map>
#include <vector>
// #include <Eigen/Sparse> Modified by Dillon Cislo 2018/08/17
#include "../../../../Eigen/Sparse"

#include "../../../core/Mesh/BaseMesh.h"
#include "../../../core/Mesh/Vertex.h"
#include "../../../core/Mesh/HalfEdge.h"
#include "../../../core/Mesh/Edge.h"
#include "../../../core/Mesh/Face.h"
#include "../../../core/Mesh/iterators.h"
#include "../../../core/Mesh/boundary.h"
#include "../../../core/Parser/parser.h"
#include "./RicciFlowMesh.h"
#include "./BaseRicciFlow.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif

namespace MeshLib
{
/*! \brief Class CTangentialRicciFlow
*
*	Algorithm for computing Ricci flow
*/
//template<typename V, typename E, typename F, typename H>
template< class V, class E, class F, class H > // Added by Dillon 2017/08/23
class CTangentialRicciFlow : public CBaseRicciFlow<V,E,F,H>
  {
  public:
    /*! \brief CTangentialRicciFlow constructor
	 *  \param pMesh the input mesh
	 */
	  CTangentialRicciFlow( CRicciFlowMesh<V,E,F,H> * pMesh );
    /*! \brief CTangentialRicciFlow destructor
	 */
	  ~CTangentialRicciFlow(){};
	/*!	Computing the metric
	 */
	void _calculate_metric();
	 /*!
	 *	Curvature flow, override 
	 */
     bool _flow( double );

  protected:

	  /*!
	   *	Calculate each edge length, has to be defined in the derivated classes
	   */
   void _length( double u1, double u2, E * e );

	/*!
	 *	Cosine law, has to be defined in the derivated classes
	 */
	double _cosine_law( double a, double b, double c ) ;
	 /*!
	  *	Normalization
	  * \param du the du vector
	  * \param n dimension of the du vector
	  */
	void _normalization( Eigen::VectorXd & du, int n );


	/*!
	 *	Calculate the edge weight
	 */
    void _calculate_edge_weight();

	/*!
	 *	Set the target curvature on each vertex
	 */
    virtual void    _set_target_curvature();

  };


//template<typename V, typename E, typename F, typename H>
template< class V, class E, class F, class H > // Added by Dillon 2017/08/23
CTangentialRicciFlow<V,E,F,H>::CTangentialRicciFlow( CRicciFlowMesh<V,E,F,H> * pMesh ):CBaseRicciFlow<V,E,F,H>( pMesh)
{
};

//Compute the edge length
//template<typename V, typename E, typename F, typename H>
template< class V, class E, class F, class H > // Added by Dillon 2017/08/23
void CTangentialRicciFlow<V,E,F,H>::_length( double u1, double u2, E * e )
{
	  e->length() = exp(u1) + exp(u2);
};


//Calculate corner angle
//template<typename V, typename E, typename F, typename H>
template< class V, class E, class F, class H > // Added by Dillon 2017/08/23
double CTangentialRicciFlow<V,E,F,H>::_cosine_law( double a, double b, double c )
{
          double cs =  ( a*a + b * b  - c * c )/( 2.0 * a * b );
          assert( cs <= 1.0 && cs >= -1.0 );
          return acos( cs );    
};


//Calculate edge weight

//template<typename V, typename E, typename F, typename H>
template< class V, class E, class F, class H > // Added by Dillon 2017/08/23
void CTangentialRicciFlow<V,E,F,H>::_calculate_edge_weight()
{
	for( typename CRicciFlowMesh<V,E,F,H>::MeshEdgeIterator eiter( this->m_pMesh ); !eiter.end(); eiter ++ )
  {
      E * e = *eiter;
      e->weight() = 0.0;
  }

	for(  typename CRicciFlowMesh<V,E,F,H>::MeshFaceIterator fiter( this->m_pMesh ) ; !fiter.end(); fiter ++ )
  {
      F * f = *fiter;

      double r[3];
      int i = 0;
	  for( typename CRicciFlowMesh<V,E,F,H>::FaceHalfedgeIterator hiter(  f ); !hiter.end(); ++hiter )
      {
         H * he = *hiter;
		 V * v = this->m_pMesh->halfedgeTarget(he);
         r[i++] = exp( v->u() );
      }

      double w = sqrt(r[0]*r[1]*r[2]/(r[0]+r[1]+r[2]));
      
	  for( typename CRicciFlowMesh<V,E,F,H>::FaceEdgeIterator eiter(f); !eiter.end(); ++ eiter )
      {
          E * e = * eiter;
          e->weight() += w/e->length();
      }
  }
};

//set target curvature

//template<typename V, typename E, typename F, typename H>
template< class V, class E, class F, class H > // Added by Dillon 2017/08/23
void CTangentialRicciFlow<V,E,F,H>::_set_target_curvature()
{
  for( typename CRicciFlowMesh<V,E,F,H>::MeshVertexIterator viter( this->m_pMesh ); !viter.end(); viter ++ )
  {
    V * v = *viter;
    v->target_k() = 0;
  }
 
  //typedef typename CRicciFlowMesh<V,E,F,H>::CLoop CLoop1;
  //typename std::vector<CRicciFlowMesh<V,E,F,H>::CLoop*>& pLs = this->m_boundary.loops();
  //typename std::vector<CLoop1*>& pLs = this->m_boundary.loops();
  typename std::vector<CLoop<V,E,F,H>*>& pLs = this->m_boundary.loops();

	
  int id = 0;

  //for( typename std::vector<CRicciFlowMesh<V,E,F,H>::CLoop*>::iterator liter = pLs.begin(); liter != pLs.end(); liter++)
  for( typename std::vector<CLoop<V,E,F,H>*>::iterator liter = pLs.begin(); liter != pLs.end(); liter++)
  {
	  typename CRicciFlowMesh<V,E,F,H>::CLoop * pL = *liter;
	  typename std::list<H*> & pHes = pL->halfedges();

	  double sum = 0;
	  double inv_sum = 0;
	
	  for( typename std::list<H*>::iterator hiter = pHes.begin(); hiter != pHes.end(); hiter ++ )
	  {
			H  * he = *hiter;
			E  * pe = this->m_pMesh->halfedgeEdge( he );

			sum  += pe->length();
			inv_sum += 1/pe->length();
 	  }

	  for( typename std::list<H*>::iterator hiter = pHes.begin(); hiter != pHes.end(); hiter ++ )
	  {
			H * ce = *hiter;
			V   * pv = this->m_pMesh->halfedgeTarget( ce );
			H * he = this->m_pMesh->vertexMostCcwInHalfEdge( pv );
			H * te =  this->m_pMesh->vertexMostClwOutHalfEdge( pv );

			double L = ( this->m_pMesh->halfedgeEdge(he)->length() + this->m_pMesh->halfedgeEdge( te )->length()  )/2.0;

		 // map all boundaries to circular holes
			if( id == 0 )
			{
				double tk = 2*PI*L/sum;
				pv->target_k() = tk;
			}
			else
			{
				double tk = -2*PI*L/sum;
				pv->target_k() = tk;
			}
	  }
	  id ++;
  }

};

//compute metric

//template<typename V, typename E, typename F, typename H>
template< class V, class E, class F, class H > // Added by Dillon 2017/08/23
void CTangentialRicciFlow<V,E,F,H>::_calculate_metric()
{

  double error = 1e-6;
  //double error = 5e-4;

  this->_calculate_edge_length();

  while( true )
  {
	  _set_target_curvature();
	  this->_Newton( error, 1 );
    //break;
      if( _flow( error ) ) break;
}      


};

//gradient flow method

//template<typename V, typename E, typename F, typename H>
template< class V, class E, class F, class H > // Added by Dillon 2017/08/23
bool CTangentialRicciFlow<V,E,F,H>::_flow( double error_threshold )
{
  int num = this->m_pMesh->numVertices();

  for( int k = 0; k < 64; k ++  )
	  {
 		 this-> _calculate_edge_length();
	      _set_target_curvature();
		  _calculate_edge_weight();

		  this->_calculate_corner_angle();
		  this->_calculate_vertex_curvature();

		  double error =  this->_calculate_curvature_error();
		  printf("Current error is %f\r\n", error );

		  if( error < error_threshold)  return true;
  	  

	  //set b 
		  for( typename CRicciFlowMesh<V,E,F,H>::MeshVertexIterator viter( this->m_pMesh ); !viter.end(); ++ viter )
		  {
		    V * v = *viter;
			double dif = v->target_k() - v->k();
            v->u() += dif * 2e-2;
		  }
    }
    return false;
};

//normalization

//template<typename V, typename E, typename F, typename H>
template< class V, class E, class F, class H > // Added by Dillon 2017/08/23
void CTangentialRicciFlow<V,E,F,H>::_normalization( Eigen::VectorXd & x, int num )
{

	double s = 0;
	for(int i = 0; i < num; i ++ )
	{
		s += x(i);
	}
	s /= num;

	for (int i = 0; i < num; i++)
	{
	 x(i) -= s;
	}
};

}

#endif // _TANGENTIAL_RICCI_FLOW_H_
