/****************************************************************
*
*
*
*
*****************************************************************/

#ifndef _selfSphericalMapping_cxx
#define _selfSphericalMapping_cxx

#include <selfSphericalMapping.h>
#include <math.h>

namespace Self
{

SphericalMapping::SphericalMapping(  MeshType& inputMesh )
{
	m_mesh = inputMesh;
}

SphericalMapping::~SphericalMapping()
{
}


// Creat Spherical Coordinate System from Cartesian Coordinate
void
SphericalMapping::request_SphericalCoords()
{
	m_mesh.add_property( SphericalCoordinates );
}

// Update Spherical Coordinate System from Cartesian Coordinate
void
SphericalMapping::update_SphericalCoords( ScalarType radius )
{	
	PointType				vertex;
	SpheriCoordType			radian;

	VertexIterType  v_it( m_mesh.vertices_begin() ), v_end( m_mesh.vertices_end() );

	while( v_it != v_end )
	{
		vertex = m_mesh.point( v_it );

		// elevation: 0~PI
		double temp = vertex[2] / radius;
        if( temp >1 ) temp = 1; else if(temp < -1) temp = -1;
        radian[0] = acos( temp );
        
		// azimuth: 0~2PI    
        if( vertex[0] != 0 )
        {
            double x = vertex[0];
            double y = vertex[1];

            radian[1] = atan( y / x);
            if(x < 0)
                radian[1] += M_PI;
            else if(radian[1] < 0)
                radian[1] += M_PI * 2;
        }
        else
        {
            if(vertex[1] == 0)
            {
                radian[1] = 0;
            }
            else
            {
                if(vertex[1] > 0)
                radian[1] = M_PI_2;
                else
                radian[1] = 3 * M_PI_2;
            }
        }
         
        m_mesh.property( SphericalCoordinates,v_it ) = radian;
        ++v_it; 
    }
}

void
SphericalMapping::release_SphericalCoords()
{
	m_mesh.remove_property( SphericalCoordinates );
}

// Creat Cartesian Coordinate from Spherical Coordinate System
void
SphericalMapping::update_CartesianCoords( ScalarType radius )
{
	PointType						vertex;
	SpheriCoordType					radian;

	VertexIterType  v_it( m_mesh.vertices_begin() ), v_end( m_mesh.vertices_end() );

	while( v_it != v_end )
	{
        radian = m_mesh.property( SphericalCoordinates,v_it );
		
		vertex[0] = (ScalarType)(radius * sin(radian[1]) * cos(radian[1]));
		vertex[1] = (ScalarType)(radius * sin(radian[1]) * sin(radian[1]));
		vertex[2] = (ScalarType)(radius * cos(radian[1]));

		m_mesh.set_point( v_it, vertex );
        ++v_it;
	}

}

// Get Spherical Coordinates
SphericalMapping::SpheriCoordType 
SphericalMapping::getSphericalCoords( VertexIterType v_it )
{
	return m_mesh.property( SphericalCoordinates,v_it );
}

SphericalMapping::SpheriCoordType 
SphericalMapping::getSphericalCoords( VertexHandleType verhandle )
{
	return m_mesh.property( SphericalCoordinates,verhandle );
}

// Stereographic Projection
// Get 2d stereographic projective coordinates from 3d spherical coordinates
SphericalMapping::SpheriCoordType 
SphericalMapping::getStereoProj( PointType &vertex )
{
	SpheriCoordType project;

	if( vertex[2] != 1)
	{
		double tmp = 1 / ( 1 - vertex[2] );
		project[0] = vertex[0] * tmp;		// Real part
		project[1] = vertex[1] * tmp;		// Imaginary part
	}
	else
	{
		if( vertex[0] > 0)		project[0] = 1.7e308;
		else					project[0] = -1.7e308;
		if( vertex[1] > 0)		project[1] = 1.7e308;
		else					project[1] = -1.7e308;
	}

	return project;
}

// Set 3d spherical coordinates to 2d stereographic projective coordinates
SphericalMapping::PointType
SphericalMapping::setStereoProj( SpheriCoordType &value )
{
	PointType vertex;
	double tmp = 1 / ( 1+ value[0] * value[0] + value[1] * value[1] );
	vertex[0] = (float)(2 * value[0] * tmp);
	vertex[1] = (float)(2 * value[1] * tmp);
	vertex[2] = (float)(1 - 2 * tmp);
	return vertex;
}

// Mobius Transformation
// Set 3 points, which are mapped to 0, 1, infinity seperately
void 
SphericalMapping::setFixPoints( SpheriCoordType &zero, 
				SpheriCoordType &one, SpheriCoordType &infinity)
{
	m_z0 = zero; m_z1 = one; m_z2 = infinity;
}

// Calculate the mobius after set the 3 fixed points
SphericalMapping::SpheriCoordType 
SphericalMapping::mobiusTransFor3Fixs( SpheriCoordType &inputPoint )
{
	SpheriCoordType tmp_z1,tmp_z2;
	tmp_z1[0] = inputPoint[0] - m_z0[0]; tmp_z1[1] = inputPoint[1] - m_z0[1];
	tmp_z2[0] = m_z1[0] - m_z2[0]; tmp_z2[1] = m_z1[1] - m_z2[1];
	SpheriCoordType tmp0 = complexNumMulti(tmp_z1, tmp_z2);
	tmp_z1[0] = inputPoint[0] - m_z2[0]; tmp_z1[1] = inputPoint[1] - m_z2[1];
	tmp_z2[0] = m_z1[0] - m_z0[0]; tmp_z2[1] = m_z1[1] - m_z0[1];
	SpheriCoordType tmp1 = complexNumMulti(tmp_z1, tmp_z2);
	tmp0 = complexNumDivi(tmp0, tmp1);
	return tmp0;
}

// 2d vector square root
SphericalMapping::SpheriCoordType 
SphericalMapping::sqrtOfVector2D( SpheriCoordType &input )
{
	double sqrtSum = sqrt(input[0] * input[0] + input[1] * input[1]);
	SpheriCoordType result;
	result[0] = sqrt(0.5 * (sqrtSum + input[0]));
	result[1] = sqrt(0.5 * (sqrtSum - input[0]));

	return result;
}

// Sort the vertices
bool 
SphericalMapping::sortSpheriCoord()
{
	MeshType::VertexIter  v_it( m_mesh.vertices_begin() );
	while( v_it != m_mesh.vertices_end() )
	{
		CoordinatesType tmp;
		tmp.first = this->getSphericalCoords( v_it )[0];
		tmp.second = this->getSphericalCoords( v_it )[1];
		m_sortedVertices.push_back( std::make_pair(tmp, v_it.handle()) );

	++v_it;
	}

	// sort the spherical coordinates
	std::sort( m_sortedVertices.begin(), m_sortedVertices.end() );

	return SUCCESS;
}


// Search for 2D nearest points
bool
SphericalMapping::searchVertices( std::vector< VertexHandleType > & filteredValue, 
		SpheriCoordType cntrValue, double recRange = 0.01 )
{
		// Limite the search range to a rectangular from (value*0.99) to (value*1.01)
	//// Limite the first coordinate
	IndexDataType upperLeft, lowerRight; 
	upperLeft.first.first = cntrValue[0] - recRange;
	upperLeft.first.second = cntrValue[1] - recRange;

	std::vector< IndexDataType >::iterator upperLeftIter =
			std::lower_bound( m_sortedVertices.begin(), m_sortedVertices.end(), upperLeft );

	lowerRight.first.first = cntrValue[0] + recRange;
	lowerRight.first.second = cntrValue[1] + recRange;
	std::vector< IndexDataType >::iterator lowerRightIter =
			std::upper_bound( upperLeftIter, m_sortedVertices.end(), lowerRight );

	//// Limite the second coordinate
	std::vector< IndexDataType >::iterator		iter(upperLeftIter);
	while( iter != lowerRightIter )
	{
		if( (iter->first).second >= upperLeft.first.second &&
							(iter->first).second <= lowerRight.first.second )
		{
			filteredValue.push_back(iter->second);
		}

	++iter;
	}

	if(filteredValue.size())
		return 1;
	else
		return 0;
}

// Search for 2D nearest points
SphericalMapping::VertexHandleType 
SphericalMapping::searchNearestVertex( SpheriCoordType value, double recRange = 0.01 )
{
	// Limite the search range to a rectangular from (value*0.99) to (value*1.01)
	//// Limite the first coordinate
	IndexDataType upperLeft, lowerRight; 
	upperLeft.first.first = value[0] - recRange;
	upperLeft.first.second = value[1] - recRange;

	std::vector< IndexDataType >::iterator upperLeftIter =
			std::lower_bound( m_sortedVertices.begin(), m_sortedVertices.end(), upperLeft );

	lowerRight.first.first = value[0] + recRange;
	lowerRight.first.second = value[1] + recRange;
	std::vector< IndexDataType >::iterator lowerRightIter =
			std::upper_bound( upperLeftIter, m_sortedVertices.end(), lowerRight );

	//// Limite the second coordinate
	std::vector< IndexDataType >				filteredValue;
	std::vector< IndexDataType >::iterator		iter(upperLeftIter);
	while( iter != lowerRightIter)
	{
		if( (iter->first).second >= upperLeft.first.second &&
							(iter->first).second <= lowerRight.first.second )
		{
			filteredValue.push_back(*iter);
		}

	++iter;
	}

	//std::cout << "Refer:" << referCoordinates[i].first << " " << referCoordinates[i].second << " " << std::endl;
	//std::cout << "low:" << upperLeft.first.first << " " << upperLeft.first.second << " " << std::endl;
	//std::cout << "up:" << lowerRight.first.first << " " << lowerRight.first.second << " " << std::endl;
	double distanRef(65535);
	VertexHandleType vhTmp;
	std::vector< IndexDataType >::iterator it(upperLeftIter);
	while( it != lowerRightIter )
	{
		double tmp0 = (it->first).first - value[0];
		double tmp1 = (it->first).second - value[1];
		double distanTmp(tmp0 * tmp0 + tmp1 * tmp1);
		if( distanTmp < distanRef)
		{
			distanRef = distanTmp;
			vhTmp = it->second;
		}

	++it;
	}
	return vhTmp;
}

/**********************************************************************/
//////////////////////// ComplexNumber Type class //////////////////////
/**********************************************************************/
SphericalMapping::SpheriCoordType 
SphericalMapping::complexNumMulti( SpheriCoordType &a, SpheriCoordType &b)
{
	SpheriCoordType result;
	result[0] = a[0]*b[0] - a[1]*b[1];
	result[1] = a[0]*b[1] + a[1]*b[0];

	return result;
}

SphericalMapping::SpheriCoordType 
SphericalMapping::complexNumDivi( SpheriCoordType &a, SpheriCoordType &b)
{
	SpheriCoordType result;
	if( b[0]*b[1] != 0)
	{
		double deno = 1 / (b[0] * b[0] + b[1] * b[1]);
		result[0] = a[0]*b[0] + a[1]*b[1];
		result[1] = a[1]*b[0] - a[0]*b[1];
		result[0] *= deno;
		result[1] *= deno;
	}
	else
	{
		if( a[0] > 0)			result[0] = 1.7e308;
		else					result[0] = -1.7e308;
		if( a[1] > 0)			result[1] = 1.7e308;
		else					result[1] = -1.7e308;
	}

	return result;
}


}// end namespace

#endif
