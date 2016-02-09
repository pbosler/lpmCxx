#include "TriCubicHermite.h"
#include "PolyMesh2d.h"
#include "Faces.h"
#include "Field.h"
#include "Edges.h"
#include "Particles.h"
#include <vector>
#include <assert.h>
#include "xyzVector.h"
#include <iostream>
#include <sstream>
#include "OutputMessage.h"
#include "Logger.h"
#include <cmath>

TriCubicHermite::TriCubicHermite( const PolyMesh2d* aMesh, const Field* aField )
{
	log = Logger::Instance();
	if ( aMesh->getFaceKind() != Faces::triangularFaces || aMesh->nDim() != 2 )
	{
		OutputMessage statusMsg("TriCubicHermite can only be used on triangulations of planar particles",
								OutputMessage::errorPriority, "TriCubicHermite constructor.");
		log->logMessage(statusMsg);
		return;
	}

	planeMesh = aMesh;
	data = aField;

	const int nTri = aMesh->nFaces();
	faceU = std::vector< std::vector<double> > ( nTri );
	faceV = std::vector< std::vector<double> > ( nTri );
	faceCoeffs = std::vector< std::vector<double> > ( nTri );
	for ( int i = 0; i < nTri; ++i)
	{
		faceU[i] = std::vector<double> ( 2, 0.0 );
		faceV[i] = std::vector<double> ( 2, 0.0 );
		faceCoeffs[i] = std::vector<double> (10 * aField->nDim(), 0.0 );
	}
	faceReady = std::vector< bool > ( nTri, false );
	faceMapsReady = false;
	allFacesReady = false;
};


double TriCubicHermite::interpolateScalar( const xyzVector loc ) const
{
	assert(data->nDim() == 1);
	const int faceIndex = planeMesh->locatePointInMeshFace( loc );
	if ( !faceReady[faceIndex] )
	{
		std::stringstream ss;
		ss << "interp coefficients not defined at face " << faceIndex;
		OutputMessage statusMsg(ss.str(), OutputMessage::errorPriority, "TriCubicHermite::interpolateScalar");
		log->logMessage(statusMsg);
		return 0.0;
	}
	if ( planeMesh->pointIsOutsideMesh( loc ) )
	{
		return 0.0;
	}
	else
	{	
		const std::vector<int> faceVerts = planeMesh->ccwVerticesAroundFace(faceIndex);
		const xyzVector v0 = planeMesh->particlePosition(faceVerts[0]);
	
		const double denom = faceU[faceIndex][0] * faceV[faceIndex][1] - faceV[faceIndex][0] * faceU[faceIndex][1];
	
		const double uLoc = ( faceV[faceIndex][1] * (loc.x - v0.x) - faceV[faceIndex][0] * (loc.y - v0.y) ) / denom;
		const double vLoc = (-faceU[faceIndex][1] * (loc.x - v0.x) + faceU[faceIndex][0] * (loc.y - v0.y) ) / denom;
						
		return faceCoeffs[faceIndex][0] + faceCoeffs[faceIndex][1] * vLoc + faceCoeffs[faceIndex][2] * vLoc * vLoc +
						faceCoeffs[faceIndex][3] * vLoc * vLoc * vLoc + faceCoeffs[faceIndex][4] * uLoc + 
						faceCoeffs[faceIndex][5] * uLoc * vLoc + faceCoeffs[faceIndex][6] * uLoc * vLoc * vLoc + 
						faceCoeffs[faceIndex][7] * uLoc * uLoc + faceCoeffs[faceIndex][8] * uLoc * uLoc * vLoc +
						faceCoeffs[faceIndex][9] * uLoc * uLoc * uLoc;
	}
};

xyzVector TriCubicHermite::interpolateVector( const xyzVector loc ) const
{
	assert(data->nDim() == 2);
	xyzVector result;
	const int faceIndex = planeMesh->locatePointInMeshFace( loc );
	if ( !faceReady[faceIndex] )
	{
		std::stringstream ss;
		ss << "interp coefficients not defined at face " << faceIndex;
		OutputMessage statusMsg(ss.str(), OutputMessage::errorPriority, "TriCubicHermite::interpolateVector");
		log->logMessage(statusMsg);
		return result;
	}
	if ( !planeMesh->pointIsOutsideMesh( loc ) )
	{
		const std::vector<int> faceVerts = planeMesh->ccwVerticesAroundFace(faceIndex);
		const xyzVector v0 = planeMesh->particlePosition(faceVerts[0]);

		const double denom = faceU[faceIndex][0] * faceV[faceIndex][1] - faceV[faceIndex][0] * faceU[faceIndex][1];
	
		const double uLoc = ( faceV[faceIndex][1] * (loc.x - v0.x) - faceV[faceIndex][0] * (loc.y - v0.y) ) / denom;
		const double vLoc = (-faceU[faceIndex][1] * (loc.x - v0.x) + faceU[faceIndex][0] * (loc.y - v0.y) ) / denom;

		const double xComp = faceCoeffs[faceIndex][0] + faceCoeffs[faceIndex][1] * vLoc + faceCoeffs[faceIndex][2] * vLoc * vLoc +
						faceCoeffs[faceIndex][3] * vLoc * vLoc * vLoc + faceCoeffs[faceIndex][4] * uLoc + 
						faceCoeffs[faceIndex][5] * uLoc * vLoc + faceCoeffs[faceIndex][6] * uLoc * vLoc * vLoc + 
						faceCoeffs[faceIndex][7] * uLoc * uLoc + faceCoeffs[faceIndex][8] * uLoc * uLoc * vLoc +
						faceCoeffs[faceIndex][9] * uLoc * uLoc * uLoc;

		const double yComp = faceCoeffs[faceIndex][10] + faceCoeffs[faceIndex][11] * vLoc + faceCoeffs[faceIndex][12] * vLoc * vLoc +
						faceCoeffs[faceIndex][13] * vLoc * vLoc * vLoc + faceCoeffs[faceIndex][14] * uLoc + 
						faceCoeffs[faceIndex][15] * uLoc * vLoc + faceCoeffs[faceIndex][16] * uLoc * vLoc * vLoc + 
						faceCoeffs[faceIndex][17] * uLoc * uLoc + faceCoeffs[faceIndex][18] * uLoc * uLoc * vLoc +
						faceCoeffs[faceIndex][19] * uLoc * uLoc * uLoc;					
		result = xyzVector(xComp, yComp);
	}
	return result;
};

xyzVector TriCubicHermite::scalarGradient( const xyzVector loc ) const
{
	assert(data->nDim() == 1);
	xyzVector result;
	const int faceIndex = planeMesh->locatePointInMeshFace( loc );
	if ( !faceReady[faceIndex] )
	{
		std::stringstream ss;
		ss << "interp coefficients not defined at face " << faceIndex;
		OutputMessage statusMsg(ss.str(), OutputMessage::errorPriority, "TriCubicHermite::scalarGradient");
		log->logMessage(statusMsg);
		return result;
	}
	if ( !planeMesh->pointIsOutsideMesh(loc) )
	{
		const std::vector<int> faceVerts = planeMesh->ccwVerticesAroundFace(faceIndex);
		const xyzVector v0 = planeMesh->particlePosition(faceVerts[0]);
	
		const double denom = faceU[faceIndex][0] * faceV[faceIndex][1] - faceV[faceIndex][0] * faceU[faceIndex][1];
		const double a =  faceV[faceIndex][1] / denom;
		const double b = -faceU[faceIndex][1] / denom;
		const double c = -faceV[faceIndex][0] / denom;
		const double d =  faceU[faceIndex][0] / denom;
	
		const double uLoc = ( faceV[faceIndex][1] * (loc.x - v0.x) - faceV[faceIndex][0] * (loc.y - v0.y) ) / denom;
		const double vLoc = (-faceU[faceIndex][1] * (loc.x - v0.x) + faceU[faceIndex][0] * (loc.y - v0.y) ) / denom;

		double p_u = faceCoeffs[faceIndex][4] + faceCoeffs[faceIndex][5] * vLoc + faceCoeffs[faceIndex][6] * vLoc * vLoc
					 + faceCoeffs[faceIndex][7] * 2.0 * uLoc + 2.0 * uLoc * faceCoeffs[faceIndex][8] * vLoc + 
					 3.0 * uLoc * uLoc * faceCoeffs[faceIndex][9];
		double p_v = faceCoeffs[faceIndex][1] + 2.0 * faceCoeffs[faceIndex][2] * vLoc + 
					 3.0 * faceCoeffs[faceIndex][3] * vLoc * vLoc + faceCoeffs[faceIndex][5] * uLoc + 
					 2.0 * faceCoeffs[faceIndex][6] * uLoc * vLoc + faceCoeffs[faceIndex][8] * uLoc * uLoc;

	
		result = xyzVector(a * p_u + b * p_v, c * p_u + d * p_v );
	}
	return result;
};


double TriCubicHermite::scalarLaplacian( const xyzVector loc ) const
{
	assert(data->nDim() == 1);
	double result;
	const int faceIndex = planeMesh->locatePointInMeshFace( loc );
	if ( !faceReady[faceIndex] )
	{
		std::stringstream ss;
		ss << "interp coefficients not defined at face " << faceIndex;
		OutputMessage statusMsg(ss.str(), OutputMessage::errorPriority, "TriCubicHermite::scalarLaplacian");
		log->logMessage(statusMsg);
		return 0.0;
	}
	if ( !planeMesh->pointIsOutsideMesh(loc) )
	{
		const std::vector<int> faceVerts = planeMesh->ccwVerticesAroundFace(faceIndex);
		const xyzVector v0 = planeMesh->particlePosition(faceVerts[0]);
	
		const double denom = faceU[faceIndex][0] * faceV[faceIndex][1] - faceV[faceIndex][0] * faceU[faceIndex][1];
		const double a =  faceV[faceIndex][1] / denom;
		const double b = -faceU[faceIndex][1] / denom;
		const double c = -faceV[faceIndex][0] / denom;
		const double d =  faceU[faceIndex][0] / denom;
	
		const double uLoc = ( faceV[faceIndex][1] * (loc.x - v0.x) - faceV[faceIndex][0] * (loc.y - v0.y) ) / denom;
		const double vLoc = (-faceU[faceIndex][1] * (loc.x - v0.x) + faceU[faceIndex][0] * (loc.y - v0.y) ) / denom;
	
		const double p_uu = 2.0 * faceCoeffs[faceIndex][7] + 2.0 * faceCoeffs[faceIndex][8] * vLoc + 
							6.0 * faceCoeffs[faceIndex][9] * uLoc;
		const double p_uv = faceCoeffs[faceIndex][5] + 2.0 * faceCoeffs[faceIndex][6] * vLoc + 
							2.0 * faceCoeffs[faceIndex][8] * uLoc;
		const double p_vv = 2.0 * faceCoeffs[faceIndex][2] + 6.0 * faceCoeffs[faceIndex][3] * vLoc +
							2.0 * faceCoeffs[faceIndex][6] * uLoc;
		result = a * a * p_uu + 2.0 * a * b * p_uv + b * b * p_vv + c * c * p_uu + 2.0 * c * d * p_uv + d * d * p_vv;
	}
	return result;	 
};

double TriCubicHermite::scalarCurl2D( const xyzVector loc )
{
	assert(data->nDim() == 2);
	return 0.0;
};

Field TriCubicHermite::interpolateScalar( const Particles& aParticles ) const
{
	assert(data->nDim() == 1);
	Field result("n/a","n/a", 1, aParticles.N());
	if ( !allFacesReady ) 
	{
		OutputMessage statusMsg("interp coefficients not defined ", OutputMessage::errorPriority, "(Field) TriCubicHermite::interpolateScalar");
		log->logMessage(statusMsg);
		return result;
	}
	
	for ( int i = 0; i < aParticles.N(); ++i)
	{
		double newVal = interpolateScalar( aParticles.physCoord(i) );
		result.insertScalarToField( i, newVal );	
	}
	return result;
};



double TriCubicHermite::vectorDivergence( const xyzVector loc )
{
	assert(data->nDim() == 2);
	return 0.0;
};

Field TriCubicHermite::scalarGradient( const Particles& aParticles ) const
{
	assert(data->nDim() == 1);
	Field result("n/a","n/a", 2, aParticles.N());
	if ( !allFacesReady ) 
	{
		OutputMessage statusMsg("interp coefficients not defined ", OutputMessage::errorPriority, "(Field) TriCubicHermite::scalarGradient");
		log->logMessage(statusMsg);
		return result;
	}
	
	for ( int i = 0; i < aParticles.N(); ++i )
	{
		xyzVector newGrad = scalarGradient( aParticles.physCoord(i) );
		result.insertVectorToField( i, newGrad );
	}
	return result;
};

Field TriCubicHermite::scalarLaplacian( const Particles& aParticles ) const
{
	assert(data->nDim() == 1);
	Field result("n/a","n/a", 1, aParticles.N());
	if ( !allFacesReady ) 
	{
		OutputMessage statusMsg("interp coefficients not defined ", OutputMessage::errorPriority, "(Field) TriCubicHermite::scalarLaplacian");
		log->logMessage(statusMsg);
		return result;
	}
	
	for ( int i = 0; i < aParticles.N(); ++i )
	{
		double newVal = scalarLaplacian( aParticles.physCoord(i) );
		result.insertScalarToField( i, newVal );
	}
	return result;
};

Field TriCubicHermite::interpolateVector( const Particles& aParticles ) const
{
	assert(data->nDim() == 2);
	Field result("n/a","n/a", 2, aParticles.N());
	if ( !allFacesReady )
	{
		OutputMessage statusMsg("interp coefficients not defined", OutputMessage::errorPriority, "(Field) TriCubicHermite::interpolatVector");
		log->logMessage(statusMsg);
		return result;
	}
	
	for ( int i = 0; i < aParticles.N(); ++i)
	{
		xyzVector newVec = interpolateVector( aParticles.physCoord(i) );
		result.insertVectorToField( i, newVec );
	}
	
	return result;
};

Field TriCubicHermite::vectorDivergence( const Particles& aParticles ) const
{
	assert(data->nDim() == 2);
	Field result("n/a","n/a", 1, aParticles.N());
	return result;
};

Field TriCubicHermite::scalarCurl2D( const Particles& aParticles ) const
{
	assert(data->nDim() == 2);
	Field result("n/a","n/a", 1, aParticles.N());
	return result;
};

void TriCubicHermite::findAllCoefficientsScalar( const Field& gradField )
{
	assert(data->nDim() == 1 );
	if ( !faceMapsReady )
		generateAllUVMaps();
	for ( int i = 0; i < planeMesh->nFaces(); ++i )
	{
		if ( !planeMesh->faceIsDivided(i) )
		{
			std::vector<int> triVerts = planeMesh->ccwVerticesAroundFace(i);
			xyzVector xy0 = planeMesh->particlePosition(triVerts[0]);
			xyzVector xyC = planeMesh->facePosition( i );
					
			faceCoeffs[i][0] = data->scalar[triVerts[0]]; 																// p00
			faceCoeffs[i][1] = faceV[i][0] * gradField.xComp[triVerts[0]] + faceV[i][1] * gradField.yComp[triVerts[0]]; // p01
			faceCoeffs[i][4] = faceU[i][0] * gradField.xComp[triVerts[0]] + faceU[i][1] * gradField.yComp[triVerts[0]]; // p10
			faceCoeffs[i][7] = 3.0 * (data->scalar[triVerts[1]] - faceCoeffs[i][4] - faceCoeffs[i][0]) - 
							   faceU[i][0] * gradField.xComp[triVerts[1]] - faceU[i][1] * gradField.yComp[triVerts[1]] 
							   + faceCoeffs[i][4];																		// p20
			faceCoeffs[i][9] = -2.0 * (data->scalar[triVerts[1]] - faceCoeffs[i][4] - faceCoeffs[i][0]) + 
			  				   faceU[i][0] * gradField.xComp[triVerts[1]] + faceU[i][1] * gradField.yComp[triVerts[1]]
			  				   - faceCoeffs[i][4];																		// p30
			faceCoeffs[i][2] = 3.0 * ( data->scalar[triVerts[2]] - faceCoeffs[i][1] - faceCoeffs[i][0]) -
							   faceV[i][0] * gradField.xComp[triVerts[2]] - faceV[i][1] * gradField.yComp[triVerts[2]]
							   + faceCoeffs[i][1];																		// p02
			faceCoeffs[i][3] = -2.0 * (data->scalar[triVerts[2]] - faceCoeffs[i][1] - faceCoeffs[i][0]) +
							   faceV[i][0] * gradField.xComp[triVerts[2]] + faceV[i][1] * gradField.yComp[triVerts[2]]
							   - faceCoeffs[i][1];																		// p03
		
			int centerParticle = planeMesh->faceParticleIndex(i);
			double denom  = faceU[i][0] * faceV[i][1] - faceV[i][0] * faceU[i][1];
			double uc =  faceV[i][1] * (xyC.x - xy0.x) - faceV[i][0] * (xyC.y - xy0.y) / denom;
			double vc = -faceU[i][1] * (xyC.x - xy0.x) + faceU[i][0] * (xyC.y - xy0.y) / denom;
			double b1 = data->scalar[centerParticle] - faceCoeffs[i][0] - faceCoeffs[i][1] * vc - faceCoeffs[i][2] * vc * vc
						- faceCoeffs[i][3] * vc * vc * vc - faceCoeffs[i][4] * uc - faceCoeffs[i][7] * uc * uc -
						faceCoeffs[i][9] * uc * uc * uc;
			double b2 = faceV[i][0] * gradField.xComp[triVerts[1]] + faceV[i][1] * gradField.yComp[triVerts[1]] - faceCoeffs[i][1];
			double b3 = faceU[i][0] * gradField.xComp[triVerts[2]] + faceU[i][1] * gradField.yComp[triVerts[2]] - faceCoeffs[i][4];			
			double hc = uc * uc * vc + uc * vc * vc - uc * vc;
			
			faceCoeffs[i][5] = (-b1 + uc * uc * vc * b2 + uc * vc * vc * b3)/hc;										// p11
			faceCoeffs[i][6] = ( b1 - uc * uc * vc * b2 - (uc * vc + uc * uc * vc) * b3)/hc;							// p12
			faceCoeffs[i][8] = ( b1 + (-uc * vc + uc * vc * vc)*b2 - uc * vc * vc * b3)/hc;								// p21
		}
	}
	allFacesReady = true;
};

void TriCubicHermite::findAllCoefficientsVector( const Field& xGrad, const Field& yGrad )
{
	assert(data->nDim() == 2 );
	if ( !faceMapsReady )
		generateAllUVMaps();
	for ( int i = 0; i < planeMesh->nFaces(); ++i )
	{
		if ( !planeMesh->faceIsDivided(i) )
		{
			std::vector<int> triVerts = planeMesh->ccwVerticesAroundFace(i);
			xyzVector xy0 = planeMesh->particlePosition(triVerts[0]);
			xyzVector xyC = planeMesh->facePosition( i );
			
			faceCoeffs[i][0] = data->xComp[triVerts[0]];														// p00x
			faceCoeffs[i][1] = faceV[i][0] * xGrad.xComp[triVerts[0]] + faceV[i][1] * xGrad.yComp[triVerts[0]]; // p01x
			faceCoeffs[i][4] = faceU[i][0] * xGrad.xComp[triVerts[0]] + faceU[i][1] * xGrad.yComp[triVerts[0]]; // p10x
			faceCoeffs[i][7] = 3.0 * (data->xComp[triVerts[1]] - faceCoeffs[i][4] - faceCoeffs[i][0]) - 
							   faceU[i][0] * xGrad.xComp[triVerts[1]] - faceU[i][1] * xGrad.yComp[triVerts[1]] 
							   + faceCoeffs[i][4];																// p20x
			faceCoeffs[i][9] = -2.0 * (data->xComp[triVerts[1]] - faceCoeffs[i][4] - faceCoeffs[i][0]) + 
			  				   faceU[i][0] * xGrad.xComp[triVerts[1]] + faceU[i][1] * xGrad.yComp[triVerts[1]]
			  				   - faceCoeffs[i][4];																// p30x
			faceCoeffs[i][2] = 3.0 * ( data->xComp[triVerts[2]] - faceCoeffs[i][1] - faceCoeffs[i][0]) -
							   faceV[i][0] * xGrad.xComp[triVerts[2]] - faceV[i][1] * xGrad.yComp[triVerts[2]]
							   + faceCoeffs[i][1];																// p02x
			faceCoeffs[i][3] = -2.0 * (data->xComp[triVerts[2]] - faceCoeffs[i][1] - faceCoeffs[i][0]) +
							   faceV[i][0] * xGrad.xComp[triVerts[2]] + faceV[i][1] * xGrad.yComp[triVerts[2]]
							   - faceCoeffs[i][1];																// p03x
			
			int centerParticle = planeMesh->faceParticleIndex(i);
			double denom  = faceU[i][0] * faceV[i][1] - faceV[i][0] * faceU[i][1];
			double uc =  faceV[i][1] * (xyC.x - xy0.x) - faceV[i][0] * (xyC.y - xy0.y) / denom;
			double vc = -faceU[i][1] * (xyC.x - xy0.x) + faceU[i][0] * (xyC.y - xy0.y) / denom;
			double b1 = data->xComp[centerParticle] - faceCoeffs[i][0] - faceCoeffs[i][1] * vc - faceCoeffs[i][2] * vc * vc
						- faceCoeffs[i][3] * vc * vc * vc - faceCoeffs[i][4] * uc - faceCoeffs[i][7] * uc * uc -
						faceCoeffs[i][9] * uc * uc * uc;
			double b2 = faceV[i][0] * xGrad.xComp[triVerts[1]] + faceV[i][1] * xGrad.yComp[triVerts[1]] - faceCoeffs[i][1];
			double b3 = faceU[i][0] * xGrad.xComp[triVerts[2]] + faceU[i][1] * xGrad.yComp[triVerts[2]] - faceCoeffs[i][4];			
			double hc = uc * uc * vc + uc * vc * vc - uc * vc;
			
			faceCoeffs[i][5] = (-b1 + uc * uc * vc * b2 + uc * vc * vc * b3)/hc;								// p11x
			faceCoeffs[i][6] = ( b1 - uc * uc * vc * b2 - (uc * vc + uc * uc * vc) * b3)/hc;					// p12x
			faceCoeffs[i][8] = ( b1 + (-uc * vc + uc * vc * vc)*b2 - uc * vc * vc * b3)/hc;						// p21x
			
			faceCoeffs[i][10] = data->yComp[triVerts[0]];														// p00y
			faceCoeffs[i][11] = faceV[i][0] * yGrad.xComp[triVerts[0]] + faceV[i][1] * yGrad.yComp[triVerts[0]]; // p01y
			faceCoeffs[i][14] = faceU[i][0] * yGrad.xComp[triVerts[0]] + faceU[i][1] * yGrad.yComp[triVerts[0]]; // p10y
			faceCoeffs[i][17] = 3.0 * (data->yComp[triVerts[1]] - faceCoeffs[i][14] - faceCoeffs[i][0]) - 
							   faceU[i][0] * yGrad.xComp[triVerts[1]] - faceU[i][1] * yGrad.yComp[triVerts[1]] 
							   + faceCoeffs[i][14];																// p20y
			faceCoeffs[i][19] = -2.0 * (data->yComp[triVerts[1]] - faceCoeffs[i][14] - faceCoeffs[i][10]) + 
			  				   faceU[i][0] * yGrad.xComp[triVerts[1]] + faceU[i][1] * yGrad.yComp[triVerts[1]]
			  				   - faceCoeffs[i][14];																// p30y
			faceCoeffs[i][12] = 3.0 * ( data->yComp[triVerts[2]] - faceCoeffs[i][11] - faceCoeffs[i][10]) -
							   faceV[i][0] * yGrad.xComp[triVerts[2]] - faceV[i][1] * yGrad.yComp[triVerts[2]]
							   + faceCoeffs[i][11];																// p02y
			faceCoeffs[i][13] = -2.0 * (data->yComp[triVerts[2]] - faceCoeffs[i][11] - faceCoeffs[i][10]) +
							   faceV[i][0] * yGrad.xComp[triVerts[2]] + faceV[i][1] * yGrad.yComp[triVerts[2]]
							   - faceCoeffs[i][11];																// p03y
			
			b1 = data->yComp[centerParticle] - faceCoeffs[i][10] - faceCoeffs[i][11] * vc - faceCoeffs[i][12] * vc * vc
						- faceCoeffs[i][13] * vc * vc * vc - faceCoeffs[i][14] * uc - faceCoeffs[i][17] * uc * uc -
						faceCoeffs[i][19] * uc * uc * uc;
			b2 = faceV[i][0] * yGrad.xComp[triVerts[1]] + faceV[i][1] * yGrad.yComp[triVerts[1]] - faceCoeffs[i][11];
			b3 = faceU[i][0] * yGrad.xComp[triVerts[2]] + faceU[i][1] * yGrad.yComp[triVerts[2]] - faceCoeffs[i][14];			
			hc = uc * uc * vc + uc * vc * vc - uc * vc;
			
			faceCoeffs[i][15] = (-b1 + uc * uc * vc * b2 + uc * vc * vc * b3)/hc;								// p11y
			faceCoeffs[i][16] = ( b1 - uc * uc * vc * b2 - (uc * vc + uc * uc * vc) * b3)/hc;					// p12y
			faceCoeffs[i][18] = ( b1 + (-uc * vc + uc * vc * vc)*b2 - uc * vc * vc * b3)/hc;					// p21y
		}
	}
}

/** @brief Generates the coordinate transformation from physical xy space to the interpolation coordinate space local to each triangular face.
*/
void TriCubicHermite::generateAllUVMaps()
{
	for ( int i = 0; i < planeMesh->nFaces(); ++i )
		generateUVMapAtFace(i);
	faceMapsReady = true;
};

/** @brief Computes the affine transformation from physical xy space to local uv coordinates at a triangular face.
@param faceIndex : face whose local transformation is to be defined.
*/
void TriCubicHermite::generateUVMapAtFace( const int faceIndex ) 
{
	if ( !planeMesh->faceIsDivided(faceIndex) )
	{
		std::vector<int> verts = planeMesh->ccwVerticesAroundFace(faceIndex);
		xyzVector v0 = planeMesh->particlePosition(verts[0]);
		xyzVector v1 = planeMesh->particlePosition(verts[1]);
		xyzVector v2 = planeMesh->particlePosition(verts[2]);
		
		faceU[faceIndex][0] = v1.x - v0.x;
		faceU[faceIndex][1] = v1.y - v0.y;
		
		faceV[faceIndex][0] = v2.x - v0.x;
		faceV[faceIndex][1] = v2.y - v0.y;
		
		faceReady[faceIndex] = true;
	}
};

xyzVector TriCubicHermite::estimateScalarGradientAtVertex( const int vertIndex, const int pow1, const int pow2 )
{
	const std::vector<int> faceIndices = planeMesh->ccwFacesAroundVertex( vertIndex );
	
	xyzVector result( 0.0, 0.0 );
	
	bool boundaryVertex = false;
	for ( int i = 0; i < faceIndices.size(); ++i)
	{
		if ( faceIndices[i] < 0 )
			boundaryVertex = true;
	}
	if ( !boundaryVertex )
	{
		xyzVector p;
		xyzVector pi, pj, pij;
		double w, w1, w2;
		for ( int i = 0; i < faceIndices.size(); ++i)
		{
			std::vector<int> vertIndices = planeMesh->ccwVerticesAroundFace( faceIndices[i] );
		
			assert( vertIndices.size() == 3 );
		
			int locInd;
			for ( int j = 0; j < 3; ++j )
			{
				if ( vertIndex == vertIndices[j] )
					locInd = j;
			}
			xyzVector v0 = planeMesh->particlePosition( vertIndex );
			v0.z = data->scalar[vertIndex];
			xyzVector v1 = planeMesh->particlePosition( vertIndices[ (locInd + 1)%3 ] );
			v1.z = data->scalar[ vertIndices[ (locInd+1)%3 ] ];
			xyzVector v2 = planeMesh->particlePosition( vertIndices[ (locInd + 2)%3 ] );
			v2.z = data->scalar[ vertIndices[ (locInd+2)%3 ] ];
		
			pi = v1 - v0;
			pj = v2 - v0;
			pij = pi.crossProduct( pj );
			double denom = pi.magnitude() * pj.magnitude();
			w1 = 1.0 / denom;
			w2 = pij.magnitude() / denom;
		
			w = pow(w1, pow1) * pow(w2, pow2);
			pij.scale(w);
			p += pij;
		}
		result = xyzVector( - p.x / p.z, - p.y / p.z );
	}	
	return result;
};

Field TriCubicHermite::estimateScalarGradientAtVertices(const int pow1, const int pow2)
{
	Field result( "estGradient", "units", 2, planeMesh->nParticles() );
	
	xyzVector gradVec;
	xyzVector zeroVec = xyzVector(0.0, 0.0);
	for ( int i = 0; i < planeMesh->nParticles(); ++i )
	{
		if ( planeMesh->particleArea(i) == 0.0 )
		{
			gradVec = estimateScalarGradientAtVertex( i, pow1, pow2 );
			result.insertVectorToField( i, gradVec );
		}
		else
		{
			result.insertVectorToField( i, zeroVec );
		}
	}
	return result;
};