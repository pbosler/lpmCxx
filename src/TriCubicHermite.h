//
//  TriCubicHermite.h
//  LPM
//
//  Created by Peter Bosler on 11/5/14.
//  Copyright (c) 2014 Peter Bosler. All rights reserved.
//

/** @file TriCubicHermite.h 
	@brief TriCubicHermite class header file
	@author Peter Bosler, Sandia National Laboratories, Multiphysics Applications
*/
#ifndef __LPM__TriCubicHermite__
#define __LPM__TriCubicHermite__

#include "PolyMesh2d.h"
#include "Field.h"
#include "Particles.h"
#include "xyzVector.h"
#include <vector>

/** @class TriCubicHermite
	@brief Cubic Hermite bivariate polynomial interpolation on triangular faces with a center particle.
	Data values are required at each particle, and gradient vectors (or estimates) are required at each vertex to define the polynomial.
	@image html startingTriLPM.png "Data locations for defining a cubic Hermite interpolating polynomial within a triangular face."
	
	
	Implements a cubic Hermite interpolation scheme on a planar mesh of triangular faces.  
	Each triangular face will have its own bivariate cubic polynomial that interpolates the Field data (vector or scalar)
	within that face.  	
	Interpolated data are continuous and once continuously differentiable at each particle but not necessarily along 
	each edge.
	
	Each instance of this class is associated with a PolyMesh2d and a Field.  
	Many instances may point to the same mesh, but each instance should have its own Field.  
	
*/
class TriCubicHermite 
{
	public :
		/** @brief Constructor.  
		Allocates memory for polynomial coefficients and local mapping transformations for each face in a PolyMesh2d
		planar mesh of triangular faces.
			@param aMesh Planar triangle mesh associated with this interpolation object
			@param aField scalar or vector field to be interpolated
		*/
		TriCubicHermite( const PolyMesh2d* aMesh, const Field* aField );
		
		/** @brief Computes an interpolated scalar at a location using bivariate cubic Hermite interpolation on triangular faces.
			Requires all coefficients to have been previously defined (using TriCubicHermite::findAllCoefficientsScalar).

			@param loc : location of desired interpolated scalar data
			@return interpolated scalar value
		*/
		double interpolateScalar( const xyzVector loc ) const;
		
		/** @brief Computes a gradient at a location using the gradient of a bivariate cubic Hermite interpolating polynomial on triangular faces.
			Requires all coefficients to have been previously defined (using TriCubicHermite::findAllCoefficientsScalar).

			@param loc : location of desired gradient data
			@return gradient approximation
		*/
		xyzVector scalarGradient( const xyzVector loc ) const;
		
		/** @brief Computes a Laplacian at a location from the Laplacian of a bivariate cubic Hermite interpolating polynomial on triangular faces.
			Requires all coefficients to have been previously defined (using TriCubicHermite::findAllCoefficientsScalar).

			@param loc : location of desired gradient data
			@return scalar Laplacian approximation
		*/
		double scalarLaplacian( const xyzVector loc ) const;
		
		double scalarCurl2D( const xyzVector loc );
		
		/** @brief Interpolates a scalar field from a planar mesh to a new set of particles.
			Requires all coefficients to have been previously defined (using TriCubicHermite::findAllCoefficientsScalar).
			Returned Field's name and units must be specified separately.

			@param aParticles : new particle set 
			@return Field object containing interpolated scalar data on the same set of particles 
		*/
		Field interpolateScalar( const Particles& aParticles ) const;

		xyzVector interpolateVector( const xyzVector loc ) const;

		double vectorDivergence( const xyzVector loc );
		
		/** @brief Approximates the gradient of a scalar field at a new set of particles.
		Requires all coefficients to have been previously defined (using TriCubicHermite::findAllCoefficientsScalar).
		Returned Field's name and units must be specified separately.

		@param aParticles : new particle set 
		*/
		Field scalarGradient( const Particles& aParticles ) const;
		
		/** @brief Approximates the Laplacian of a scalar field at a new set of particles.
		Requires all coefficients to have been previously defined (using TriCubicHermite::findAllCoefficientsScalar).
		Returned Field's name and units must be specified separately.

			@param aParticles : new particle set 
		*/
		Field scalarLaplacian( const Particles& aParticles ) const;
		
		Field interpolateVector( const Particles& aParticles ) const;
		Field vectorDivergence( const Particles& aParticles ) const;
		Field scalarCurl2D( const Particles& aParticles ) const;
		
		void resetFaces();
		
		void setGradientFromField( const Field& gradField );
		
		/** @brief Defines the interpolating coefficients at each triangular face of a planar mesh object.

			@param gradField Field object containing the gradient vector at each vertex particle of the scalar field associated with this
			 TriCubicHermie interpolation object.  
			 May either be estimated (using TriCubicHermite::estimateScalarGradientAtVertices) or defined separately.
		*/
		void findAllCoefficientsScalar(const Field& gradField);
		
		/** @brief Defines the interpolating coefficients of two polynomials (one for the vector's x-component, one for the y)
			
			@param xGrad Field object containing the gradient vector of the source data's x-component
			@param yGrad Field object containing the gradient vector of the source data's y-component
		*/
		void findAllCoefficientsVector( const Field& xGrad, const Field& yGrad );
		
		Field estimateScalarGradientAtVertices( const int pow1 = 0, const int pow2 = 0);
		
	protected :
		std::vector< std::vector<double> > faceU;
		std::vector< std::vector<double> > faceV;
		std::vector< std::vector<double> > faceCoeffs;
		std::vector< bool > faceReady;
		bool faceMapsReady;
		bool allFacesReady;
		
		void generateUVMapAtFace( const int faceIndex);		
		void generateAllUVMaps();
		
		xyzVector estimateScalarGradientAtVertex( const int vertIndex, const int pow1, const int pow2 );
				
		const PolyMesh2d* planeMesh;
		const Field* data;
		
		Logger* log;
		
};

#endif