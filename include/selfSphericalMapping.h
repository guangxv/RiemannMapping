/*****************************************************************************
*
* Program:   selfSphericalMapping
* Module:    Auxiliary functions of a sphere mesh 
*
* Function:  Transform Cartesian coordinates and spherical coordinates of a 
*			 sphere mesh.
*   		 Mobius Transformation for three-fixed point of spherical conformal
*			 mapping.
*   		 Search for the nearest points on spherical coordinates.
*         
* Date:      $Date: 2013-01-17$
* Version:   $Revision: 1.0 $
*
* Reference:          
=========================================================================*/
#ifndef _selfSphericalMapping_h
#define _selfSphericalMapping_h

#include <iostream>
#include <vector>
// --------------------
#include <selfTriangularMesh.h>

namespace Self
{
	class SphericalMapping
	{
	public:
		typedef OpenMesh::TriangularMesh		MeshType;
		typedef MeshType::Point					PointType;
		typedef MeshType::VertexIter			VertexIterType;
		typedef MeshType::VertexVertexIter		VerVerIterType;
		typedef MeshType::EdgeIter				EdgeIterType;
		typedef MeshType::HalfedgeIter			HEdgeIterType;
		typedef MeshType::FaceIter				FaceIterType;
		typedef MeshType::FaceVertexIter		FaceVerIterType;
		typedef MeshType::VertexHandle			VertexHandleType;
		typedef MeshType::Scalar				ScalarType;
		typedef OpenMesh::Vec3d					VectorType;

		typedef OpenMesh::Vec2d					SpheriCoordType; //Spherical Coordinates Type
											//[0]: elevation(0~PI); [1]: azimuth(0~2PI)

		typedef std::pair<double,double> CoordinatesType;
		typedef std::pair<CoordinatesType, VertexHandleType >	 IndexDataType;

	public:

		/*************************************************************************/
		SphericalMapping( MeshType& );
		~SphericalMapping();

		/******************** Spherial Coordinate System *************************/
		
		// Tansform between the Spherical Coordinate and the Cartesian Coordinate
		void request_SphericalCoords( void );
		void update_SphericalCoords( ScalarType radius = 1 );
		void release_SphericalCoords( void );

		void update_CartesianCoords( ScalarType radius = 1 );
	    
		void sphericalRotation( unsigned int direction, ScalarType radian);
	    
		SpheriCoordType getSphericalCoords( VertexIterType verIter );
		SpheriCoordType getSphericalCoords( VertexHandleType verhandle );

	
		/************************ Stereographic Projection **************************/

		// Get 2d stereographic projective coordinates from 3d spherical coordinates
		SpheriCoordType getStereoProj( PointType &vertex );
		// Set 3d spherical coordinates to 2d stereographic projective coordinates
		PointType setStereoProj( SpheriCoordType &value );

		/************************ Mobius Transformation ***************************/

		void setFixPoints( SpheriCoordType &zero, SpheriCoordType &one, SpheriCoordType &infinity);
		SpheriCoordType mobiusTransFor3Fixs( SpheriCoordType &inputPoint );

		/************************ Search for the nearest points ***************************/
		// Sort the spherical coordinates and store it to m_sortedVertices
		bool sortSpheriCoord();

		// Search vertices in a special range
		bool searchVertices( std::vector< VertexHandleType > &,
										SpheriCoordType cntrValue,double range);

		// Inside of the range, search the nearest point to the refer value after sort
		// Return the i'th item
		VertexHandleType searchNearestVertex( SpheriCoordType value, double range );

		/*************** Variable ******************************/
	private:
		MeshType m_mesh;

		// first: elevation(0~PI); second: azimuth(0~2PI)
		OpenMesh::VPropHandleT< SpheriCoordType > SphericalCoordinates;

		OpenMesh::VPropHandleT< SpheriCoordType > StereographicProjection;

		//Search the nearest vertices < <spherical coordinates>, vertical index>
		std::vector< IndexDataType >	m_sortedVertices;

		// Coefficients for Mobius transformation
		SpheriCoordType m_z0, m_z1, m_z2;

		// Complex Number Operation			
		SpheriCoordType sqrtOfVector2D( SpheriCoordType & ); // Sqrt of Vector 2D
		SpheriCoordType complexNumMulti( SpheriCoordType &a, SpheriCoordType &b);
		SpheriCoordType complexNumDivi( SpheriCoordType &a, SpheriCoordType &b); // a/b

	};


}	// end of the namespace

#include "selfSphericalMapping.cxx"
#endif
