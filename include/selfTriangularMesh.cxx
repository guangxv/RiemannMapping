/****************************************************************
*
*
*
*
*****************************************************************/

#ifndef _selfTriangularMesh_cxx
#define _selfTriangularMesh_cxx

#include <selfTriangularMesh.h>
#include <math.h>
namespace OpenMesh
{

// String constants k
// String energy constants

void 
TriangularMesh::request_StringConstants( void )
{
	this->add_property( STRING_CONSTANTS );
}


void 
TriangularMesh::update_StringConstants( unsigned int option )
{
	switch(option){

	// k = cot_alfa + cot_belt
	case 0:
		{
			// The odd halfedge is the opposite of the even halfedge
			HEdgeIterType he_it(this->halfedges_begin()), e_end(this->halfedges_end());		
			while( he_it != e_end )
			{
				PointType v1, v2, v3;
				VectorType temp2;
				ScalarType temp1;
				double alfa, belt(0);
				EdgeHandle e_handle(this->edge_handle(he_it));

				v1 = this->point(this->from_vertex_handle(he_it));
				v2 = this->point(this->to_vertex_handle(he_it));
				v3 = this->point(this->opposite_vh(he_it));

				temp1 = (v1-v3) | (v2-v3); // dot product
				temp2 = (v1-v3) % (v2-v3); // cross product
				alfa = temp1 / temp2.norm();
				++he_it;

				if(!this->is_boundary(he_it.handle())) // There is not belt in
				{										// boundary half edge
					v3 = this->point(this->opposite_vh(he_it));
					temp1 = (v1-v3) | (v2-v3);				
					temp2 = (v1-v3) % (v2-v3); // cross product
					belt = temp1 / temp2.norm();
					++he_it;
				}

				this->property( STRING_CONSTANTS, e_handle ) = alfa + belt;
			}
			break;
		}
	default:
		{ break; }
	}
}

// Get string constants for edges
// Must run request_StringConstants() and update_StringConstants() previously

double
TriangularMesh::getStringConstants( EdgeIterType e_it )
{
	double k = this->property( STRING_CONSTANTS, e_it );
	if( k < 0) k = 0;
	return 0.5 * k;
}


double
TriangularMesh::getStringConstants( VerEdgeIterType e_it )
{
	double k = this->property( STRING_CONSTANTS, e_it );
	if( k < 0) k = 0;
	return 0.5 * k;
}


void 
TriangularMesh::release_StringConstants( void )
{
	this->remove_property( STRING_CONSTANTS );
}

// Laplacian Operator
TriangularMesh::VectorType
TriangularMesh::laplacian( VertexIterType v_it )
{
	VectorType tempValue(0,0,0);
	VerVerIterType	vv_it(this->vv_iter( v_it ));
	while( vv_it )
	{
		tempValue += 
			(VectorType)(this->point( vv_it ) - this->point( v_it ));
		++vv_it;	
	}
	return tempValue;
}

TriangularMesh::VectorType
TriangularMesh::laplacian( VertexHandleType vh )
{
	VectorType tempValue(0,0,0);
	VerVerIterType	vv_it(this->vv_iter( vh ));
	while( vv_it )
	{
		tempValue += 
			(VectorType)(this->point( vv_it ) - this->point( vh ));
		++vv_it;	
	}
	return tempValue;
}


// Voronoi area element for vertex
double 
TriangularMesh::vertexAreaElement( VertexIterType v_it )
{
	// Voronoi area
	double area(0);
	VerEdgeIterType	ve_it(this->ve_iter( v_it ));
	while( ve_it )
	{
		double stringCons = this->getStringConstants( ve_it );
		area += stringCons * this->calc_edge_sqr_length(ve_it);
		++ve_it;
	}
	return 0.25 * area;
}
double 
TriangularMesh::vertexAreaElement( VertexHandleType vh )
{
	// Voronoi area
	double area(0);
	VerEdgeIterType	ve_it(this->ve_iter( vh ));
	while( ve_it )
	{
		double stringCons = this->getStringConstants( ve_it );
		area += stringCons * this->calc_edge_sqr_length(ve_it);
		++ve_it;
	}
	return 0.25 * area;
}

// Mean Curvature
double
TriangularMesh::meanCurvature( VertexIterType v_it )
{
	VerEdgeIterType	ve_it(this->ve_iter( v_it ));

	// Voronoi area
	double area = this->vertexAreaElement( v_it );
	// Laplace-Beltrami operator
	VectorType tmpVec, sumVec;
	sumVec[0] = 0; sumVec[1] = 0; sumVec[2] = 0;
	VerVerIterType	vv_it(this->vv_iter( v_it ));
	while( vv_it )
	{
		tmpVec = this->point(vv_it) - this->point(v_it);
		double stringCons = this->property( STRING_CONSTANTS, ve_it );
		tmpVec *= stringCons;
		sumVec += tmpVec;
		++vv_it;
	}
	
	if(area != 0)
	{
		return 0.25 * sumVec.norm() / area;
	}
	else
	{
		std::cerr << "curvature calculation is error!" << std::endl;
		return 0;
	}
}

double
TriangularMesh::meanCurvature( VertexHandleType vh )
{
	VerEdgeIterType	ve_it(this->ve_iter( vh ));

	// Voronoi area
	double area = this->vertexAreaElement( vh );
	// Laplace-Beltrami operator
	VectorType tmpVec, sumVec;
	sumVec[0] = 0; sumVec[1] = 0; sumVec[2] = 0;
	VerVerIterType	vv_it(this->vv_iter( vh ));
	while( vv_it )
	{
		tmpVec = this->point(vv_it) - this->point(vh);
		double stringCons = this->property( STRING_CONSTANTS, ve_it );
		tmpVec *= stringCons;
		sumVec += tmpVec;
		++vv_it;
	}
	
	if(area != 0)
	{
		return 0.25 * sumVec.norm() / area;
	}
	else
	{
		std::cerr << "curvature calculation is error!" << std::endl;
		return 0;
	}
}

// Gaussian Curvature
double
TriangularMesh::gaussianCurvature( VertexIterType v_it )
{
	// Voronoi area
	double area = this->vertexAreaElement( v_it );
	// Angles
	double angle(0);
	PointType p0,p1,p2;
	VerVerIterType	vv_it(this->vv_begin(v_it));
	p0 = this->point(v_it);

	p1 = this->point(vv_it);
	do{
		p2 = this->point(++vv_it);
		VectorType e0, e1;
		e0 = p1 - p0;
		e1 = p2 - p0;
		angle +=acos( e0.normalize() | e1.normalize() );
		p1 = p2;
	}while(vv_it);

	if(area != 0)
	{
		return (2. * M_PI - angle) / area;
	}
	else
	{
		std::cerr << "curvature calculation is error!" << std::endl;
		return 0;
	}
}

double
TriangularMesh::gaussianCurvature( VertexHandleType vh )
{
	// Voronoi area
	double area = this->vertexAreaElement( vh );
	// Angles
	double angle(0);
	PointType p0,p1,p2;
	VerVerIterType	vv_it(this->vv_begin(vh));
	p0 = this->point(vh);

	p1 = this->point(vv_it);
	do{
		p2 = this->point(++vv_it);
		VectorType e0, e1;
		e0 = p1 - p0;
		e1 = p2 - p0;
		angle +=acos( e0.normalize() | e1.normalize() );
		p1 = p2;
	}while(vv_it);

	if(area != 0)
	{
		return (2. * M_PI - angle) / area;
	}
	else
	{
		std::cerr << "curvature calculation is error!" << std::endl;
		return 0;
	}
}

}// end namespace

#endif
