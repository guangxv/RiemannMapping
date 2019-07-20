/*****************************************************************************
*
* Program:   selfTriangularMesh
* Module:    This module adds the geometric calculators to TriMesh_ArrayKernel
*
* Function:  delaunayConfirmingFilter is used to check the delaunary of a 
*   		 genus zero mesh. Rectify the non-delaunary by simple flipping. 
*   		 Mean Smoothing method is described in Openmesh guide.
*   		 StringConstants and Laplacian operator refer Gu`s paper
*         
* Date:      $Date: 2012-11-04$
* Version:   $Revision: 1.0 $
*
* Reference: X.Gu, Y. Wang et al., Genus Zero Surface Conformal Mapping and
*   		 Its Application to Brain surface Mapping. IEEE Trans. on Medical
*   		 Imaging, Vol.23,No. 8, 2004.
*   		 Meyer et.al.,"Discrete Differential-Geometry Operators for 
*            Triangulated 2-Manifolds"          
=========================================================================*/
#ifndef _selfTriangularMesh_h
#define _selfTriangularMesh_h

#include <iostream>
#include <vector>
//#include <vnl/vnl_math.h>
// --------------------
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
//#include <OpenMesh/Core/Mesh/Traits.hh>


#define SUCCESS 1
#define FAILURE 0


namespace OpenMesh
{
	class TriangularMesh : public TriMesh_ArrayKernelT<>
	{
	public:
		typedef TriMesh_ArrayKernelT<>				MeshType;
		typedef MeshType::Point						PointType;
		typedef MeshType::VertexIter				VertexIterType;
		typedef MeshType::VertexVertexIter			VerVerIterType;
		typedef MeshType::VertexEdgeIter			VerEdgeIterType;
		typedef MeshType::EdgeIter					EdgeIterType;
		typedef MeshType::HalfedgeIter				HEdgeIterType;
		typedef MeshType::FaceIter					FaceIterType;
		typedef MeshType::FaceVertexIter			FaceVerIterType;
		typedef MeshType::VertexHandle				VertexHandleType;
		typedef MeshType::Scalar					ScalarType;
		typedef Vec3d								VectorType;

		/**********************************************************************/
		/////////////Basic operations on genus zero tetrahedral mesh///////////
		/**********************************************************************/

		// String energy constants
		void request_StringConstants( void );
		void release_StringConstants( void );
		// option 0: k = cot_alfa + cot_belt
		void update_StringConstants( unsigned int option = 0 );

		// Must run request_StringConstants() and update_StringConstants() previously
		double getStringConstants( EdgeIterType e_it );
		double getStringConstants( VerEdgeIterType e_it );

		/********************** Discrete Surface Operators ************************************/
		// We must call request_StringConstants() and update_StringConstants()
		// befor using the functions below !!!!!!!!!!!!

		// Laplacian function
		VectorType laplacian( VertexIterType vertexIter );
		VectorType laplacian( VertexHandleType vertexHandle );

		// Discrete Mean Curvature for Delaunay mesh (not exist obtuse)
		// Meyer et.al.,"Discrete Differential-Geometry Operators for Triangulated 2-Manifolds"
		double vertexAreaElement( VertexIterType vertexIter );
		double vertexAreaElement( VertexHandleType vertexHandle );

		double meanCurvature( VertexIterType vertexIter );
		double meanCurvature( VertexHandleType vertexHandle );

		double gaussianCurvature( VertexIterType vertexIter );
		double gaussianCurvature( VertexHandleType vertexHandle );

	private:
		EPropHandleT< double > STRING_CONSTANTS;


	};


}	// end of the namespace

#include "selfTriangularMesh.cxx"
#endif
