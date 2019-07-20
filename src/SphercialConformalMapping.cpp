/*****************************************************************************
*
* Program:   Spherical Confromal Mapping
* Module:    Read Smoothed surface mesh
*            Calculate Tuette Mapping 
*            Calculate Spherical Conformal Mapping
*            (It is an optional that choose the zero-mass constraint or 3-fixed 
*   		 points constraint for spherical conformal mapping)
*            Write the mapping sphere
*
* Function:  Calculate the Tuette Mapping and Spherical conformal mapping
*   		 for a genus zero surface. We can get the spherical conformal
*   		 mapping from surface mesh or its Tuette map. The processing is
*   		 Smoothed mesh --> Tuette mapping --> spherical conformal mapping.
*   		 If trig ZERO_MASS_CENTER, the conformal mapping is operated using
*   		 zero-mass constraint. Else, use 3-fixed points constraint that
*   		 map three specific points to South pole, one point on centerline,
*   		 North pole of the unit sphere.
*   		 The input mesh should genus zero, i.e. there is not hole or island
*   		 on the surface.
*         
* Date:      $Date: 2017-02-09$
* Version:   $Revision: 1.1 $
*             Add main() arguments
*
* Reference: X.Gu, Y. Wang et al., Genus Zero Surface Conformal Mapping and
*   		 Its Application to Brain surface Mapping. IEEE Trans. on Medical
*   		 Imaging, Vol.23,No. 8, 2004.
*   		 X. Gu,S.T. Yau, Global Conformal Surface Parameterization.
*   		 Eurographics Symposium on Geometry Processing,2003.
          
=========================================================================*/
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <iostream>
#include <vector>

// --------------------
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <selfTriangularMesh.h>
#include <selfSphericalMapping.h>

/***********************************************************************
* main
************************************************************************/

int main(int argc, char * argv[] )
{
	if( argc < 4 )
    {
    	std::cerr << "Arguments Missing. " << std::endl;
    	std::cerr <<
      	"Usage:  IterativeClosestPoint3   fixedPointsFile  movingPointsFile "
      	<< std::endl;
    	return EXIT_FAILURE;
    }

    //
	// Type definition
	//
	typedef OpenMesh::TriangularMesh				MeshType;
	typedef MeshType::Point							PointType;
	typedef MeshType::VectorType					VectorType;
	typedef MeshType::Scalar						ScalarType;

	typedef MeshType::VertexIterType				VertexIterType;
	typedef MeshType::HalfedgeHandle				HalfedgeHandleType;

	typedef Self::SphericalMapping					SphericalMappingFilterType;


	// Define Variable
	MeshType							inputMesh, outputMesh;

	// Global variables
	double energy0(0), deltEnergy;
	const double tuetteThreshold = atof(argv[3]); // (0.0005)
	const double conformalThreshold = atof(argv[4]); // (0.00001)
	double delt_t( 0.2 );
	//-----------------------------------------------------------------
	// Read
	//-----------------------------------------------------------------
 
	// read mesh from stdin
	if ( !OpenMesh::IO::read_mesh(inputMesh, argv[1] ) )
	{
		std::cerr << "Error: Cannot read mesh from Case" << argv[1] << std::endl;
		EXIT_FAILURE;
	}
	std::cout << "Reading Case" << argv[1] << std::endl;

	outputMesh = inputMesh;

	/****************** Variable **************************************/
	// set output mapping mesh
  	MeshType::VertexIter  v_it0, v_end0(inputMesh.vertices_end());
	MeshType::EdgeIter	  e_it0, e_end0(inputMesh.edges_end());
	MeshType::VertexIter  v_it1, v_end1(outputMesh.vertices_end());
	MeshType::EdgeIter	  e_it1, e_end1(outputMesh.edges_end());
	//MeshType::FaceIter	  f_it1, f_end1(outputMesh.faces_end());

/***********************************************************************
* Get the center to 0
************************************************************************/

	//PointType cntrPoint(0,0,0);

	//v_it0 = inputMesh.vertices_begin();
	//while( v_it0 != v_end0 )
	//{
	//	cntrPoint += inputMesh.point(v_it0);
	//++v_it0;
	//}
	//cntrPoint /= inputMesh.n_vertices();

	//v_it0 = inputMesh.vertices_begin();
	//while( v_it0 != v_end0 )
	//{
	//	inputMesh.point(v_it0) -= cntrPoint;
	//++v_it0;
	//}

/*************************************************************************************
* Spherical Tuette Mapping
* Procudure:
*	1. Compute Gauss map; initial Tuette energy E0
*	2. For each vertex, compute absolute derivative Dt
*	3. Update t(v)
*	4. Compute the Tuette energy E
*	5. If |E - E0| < threshold, reture t. Else, assign E to E0 repeat steps 2-5
**************************************************************************************/

	//
	// Initialization
	//Gauss map
	inputMesh.request_face_normals();	//we need face normals previously
	inputMesh.request_vertex_normals();
	inputMesh.update_normals();

	v_it0 = inputMesh.vertices_begin();
	v_it1 = outputMesh.vertices_begin();
	while( v_it0 != v_end0 )
	{
		outputMesh.set_point( v_it1, inputMesh.normal(v_it0) );
	++v_it0; ++v_it1;
	}

	// Release memory
	inputMesh.release_face_normals();
	inputMesh.release_vertex_normals();


	// Half of String energy E0
	e_it1 = outputMesh.edges_begin();
	while( e_it1 != e_end1 )
	{
		energy0 += outputMesh.calc_edge_sqr_length( e_it1 );
	++e_it1;
	}
	std::cout << "Tuette Energy0: " << energy0 << std::endl;

	//
	// AbsoluteDerivative
	// repeat steps 2)-5)
	std::cout << "Tuette Mapping..." << std::endl;
	do
	{
		// AbsoluteDerivative
		v_it1 = outputMesh.vertices_begin();
		while( v_it1 != v_end1 )
		{
			VectorType normal(outputMesh.point(v_it1));
			VectorType laplacian(outputMesh.laplacian(v_it1));
			// Get Derivative
			VectorType normalLaplacian = (laplacian | normal) * normal;	//<lap, normal> normal
			VectorType derivative = laplacian - normalLaplacian;			//der = lap - normalLap
			// Adjust the Vertices
			PointType  adjVertex;
			adjVertex = normal + delt_t * derivative;	//delt -= der
			// Update the position of Vertices
			outputMesh.set_point( v_it1, adjVertex.normalize() );

		++v_it1;
		}

		// Calculate the energy
		double energy(0.0);
		e_it1 = outputMesh.edges_begin();
		while( e_it1 != e_end1 )
		{
			energy += outputMesh.calc_edge_sqr_length( e_it1 );
		++e_it1;
		}

		deltEnergy = energy0 - energy;
		energy0 = energy;

		// adjust speed of variation
		if(deltEnergy < 0)
		{
			std::cout << "Energy " << energy0 << " deltEnergy: " << deltEnergy 
				<< " delt: "<< delt_t <<  std::endl;
			delt_t = delt_t * 0.5;
		}

	}while( abs(deltEnergy) > tuetteThreshold ); // Condition of conversation

	std::cout << "Final Tutte Energy: " << energy0 << " delt: " << deltEnergy << std::endl;

/*************************************************************************************
* Spherical Conformal Mapping
* Procudure:
*	1. Compute Gauss map; initial Tuette energy E0
*	2. For each vertex, compute absolute derivative Dt
*	3. Update t(v)
*	4. Compute Mobius transformation
*	4. Compute the harmonic energy E
*	5. If |E - E0| < threshold, reture t. Else, assign E to E0 repeat steps 2-5
**************************************************************************************/

	//
	// Initialization
	//
	// Require Harmonic energy
	inputMesh.request_StringConstants();
	inputMesh.update_StringConstants(0);
	// Initial energy E0 = Tutte energy

	//
	// Area element
	OpenMesh::VPropHandleT< double > AREA_ELEMENT;
	outputMesh.add_property(AREA_ELEMENT);

		v_it0 = inputMesh.vertices_begin();
		v_it1 = outputMesh.vertices_begin();
		while( v_it1 != v_end1 )
		{
			outputMesh.property(AREA_ELEMENT, v_it1) = inputMesh.vertexAreaElement(v_it0);
			++v_it0;
		++v_it1;
		}
	
	//
	// AbsoluteDerivative
	std::cout << "Spherical Conformal Mapping..." << std::endl;

	do
	{
		// AbsoluteDerivative
		v_it1 = outputMesh.vertices_begin();
		while( v_it1 != v_end1 )
		{
			VectorType normal(outputMesh.point(v_it1));
			VectorType laplacian(outputMesh.laplacian(v_it1));
			// Get Derivative
			VectorType normalLaplacian = (laplacian | normal) * normal;	//<lap, normal> normal
			VectorType derivative = laplacian - normalLaplacian; //der = lap - normalLap
			// Adjust the Vertices
			PointType  adjVertex;
			adjVertex = normal + delt_t * derivative;	//delt -= der
			// Update the position of Vertices
			outputMesh.set_point( v_it1, adjVertex.normalize() );

		++v_it1;
		}

		// Initial
		VectorType massCenter(0,0,0);					// Mass center
		// Calculate the mass center
		double mass(0);
		v_it1 = outputMesh.vertices_begin();
		while( v_it1 != v_end1 )
		{
			double areaEle = outputMesh.property(AREA_ELEMENT, v_it1);
			VectorType element;
			element[0] = outputMesh.point(v_it1)[0];
			element[1] = outputMesh.point(v_it1)[1];
			element[2] = outputMesh.point(v_it1)[2];
			massCenter += element * areaEle;
			mass += areaEle;
		++v_it1;
		}
		massCenter /= mass;
	
		// Update the position of Vertices
		v_it1 = outputMesh.vertices_begin();
		while( v_it1 != v_end1 )
		{
			PointType adjVertex;
			adjVertex = outputMesh.point(v_it1) - (PointType)massCenter;
			outputMesh.set_point( v_it1, adjVertex.normalize());

		++v_it1;
		}

		// Calculate the energy
		double energy(0.0);
		e_it0 = inputMesh.edges_begin();
		e_it1 = outputMesh.edges_begin();
		while( e_it1 != e_end1 )
		{
			double tmp = inputMesh.getStringConstants(e_it0);
			energy += tmp * outputMesh.calc_edge_sqr_length( e_it1 );
		++e_it0;++e_it1;
		}

		deltEnergy = energy0 - energy;
		energy0 = energy;
		// adjust speed of variation
		if(deltEnergy < 0)
		{
			std::cout << "Energy " << energy0  << " delt: "<< delt_t <<  std::endl;
			delt_t *= 0.98;
		}
	}while( fabs(deltEnergy) > conformalThreshold ); 

	inputMesh.release_StringConstants();
	//outputMesh.remove_property(AREA_ELEMENT);
	std::cout << "Final Conformal Energy: " << energy0 << " delt:" << deltEnergy << std::endl;


	//
	// Output
	if ( ! OpenMesh::IO::write_mesh( outputMesh, argv[2]) )
	{
		std::cerr << "Error: cannot write mesh to Case" << argv[2] << std::endl;
		return EXIT_FAILURE;
	}

/***********************************************************************
* End
************************************************************************/
    //system("PAUSE");
	return EXIT_SUCCESS;
}