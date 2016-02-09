#include "PolyMesh2d.h"
#include "xyzVector.h"
#include "Field.h"
#include "OutputMessage.h"
#include "Logger.h"
#include "TriCubicHermite.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>

int main ( int argc, char* argv[] )
{
	const bool verbose = false;
	const int procRank = 0;
	const int numProcs = 1;
	const OutputMessage::priority loggingLevel = OutputMessage::debugPriority;
	
	Logger* exeLog = Logger::Instance(loggingLevel, procRank, numProcs);
	
	const double maxR = 1.0;
	
	double nInv[7];
	double interpMaxErr[7];
	double interpAvgErr[7];
	double gradMaxErr[7];
	double gradAvgErr[7];
	double lapMaxErr[7];
	double lapAvgErr[7];
	
	
	const int scalarDim = 1;
	const int vectorDim = 2;
	
	std::stringstream ss;
	std::string s;
	
	for ( int iNest = 1; iNest < 8; ++iNest )
	{
		//
		// construct a mesh to use for testing
		//
		PolyMesh2d pMesh( iNest, PolyMesh2d::triHexSeed, maxR, procRank, numProcs);
		std::cout << "iNest = " << iNest << ", nParticles = " << pMesh.nParticles() << ", nFaces = " << pMesh.nFaces() << std::endl;
		
		if ( verbose )
		{
			LongMessage meshStatus("planar mesh built", OutputMessage::tracePriority, "main", pMesh.getInfo() );
			exeLog->logMessage(meshStatus);
		}
		
		nInv[iNest-1] = 1.0/pMesh.nFaces();
		
		//
		// construct scalar field to use for testing
		//
		Field gaussScalar("Gaussian", "n/a", scalarDim, pMesh.nParticles() );
		Field exactGradient( "grad Gaussian", "n/a", vectorDim, pMesh.nParticles() );
		Field exactLaplacian( "laplacian Gaussian", "n/a", scalarDim, pMesh.nParticles() );
		Particles* pp = pMesh.getParticles();
		const double b = 8.0;
		for (int i = 0; i < pMesh.nParticles(); ++i)
		{
			double x = pp->x[i];
			double y = pp->y[i];
			double expPart = exp(- b*b * ( x*x + y*y ) );
			gaussScalar.insertScalarToField( i, expPart );
			
			exactGradient.insertVectorToField( i, -2.0 * b * b * x * expPart, - 2.0 * b * b * y * expPart );
			
			exactLaplacian.insertScalarToField(i, 4.0 * b * b * expPart * (-1.0 + b*b * ( x*x + y*y )) );
		}
		
		std::cout << "exact Fields defined... " << std::endl;
		
		ss.str(s);
		ss << "gaussScalar_nest_" << iNest << ".m";
		gaussScalar.outputForMatlab( ss.str(), *pp );
		
		std::cout << "source data output complete... " << std::endl;
		
		//
		// initialize interpolation
		//
		TriCubicHermite gaussScalarInterp( &pMesh, &gaussScalar);
		//
		// compute gradient estimates for Hermite interpolation, set interpolation coefficients
		//
		const double w1 = 2.0;
		const double w2 = 1.0;
		Field gaussGrad = gaussScalarInterp.estimateScalarGradientAtVertices(w1, w2);
		gaussScalarInterp.findAllCoefficientsScalar( gaussGrad );
		
		std::cout << "interpolation coefficients set... " << std::endl;
		//
		// interpolate
		//
		Field interpolatedScalar("GaussInterp", "n/a", scalarDim, pMesh.nParticles() );
		Field interpolatedGradient("GaussInterpGrad", "n/a", vectorDim, pMesh.nParticles() );
		Field interpolatedLaplacian("GaussInterpLaplacian", "n/a", scalarDim, pMesh.nParticles() );
		
		for ( int i = 0; i < pMesh.nParticles(); ++i)
		{
			xyzVector loc( pp->x[i], pp->y[i] );
			
			interpolatedScalar.insertScalarToField( i, gaussScalarInterp.interpolateScalar( loc ) );			
			
			interpolatedGradient.insertVectorToField(i, gaussScalarInterp.scalarGradient( loc ) );
			
			interpolatedLaplacian.insertScalarToField( i, gaussScalarInterp.scalarLaplacian( loc ) );
		}
		std::cout << "interpolation complete... calculating error ..." << std::endl;
		
		Field interpError("interpolation error", "n/a", scalarDim, pMesh.nParticles() );
		Field gradError( "gradient estimation error", "n/a", vectorDim, pMesh.nParticles() );
		Field lapError("laplacian error", "n/a", scalarDim, pMesh.nParticles() );
		
		interpMaxErr[iNest-1] = 0.0;
		interpAvgErr[iNest-1] = 0.0;
		gradMaxErr[iNest-1] = 0.0;
		gradAvgErr[iNest-1] = 0.0;
		lapMaxErr[iNest-1] = 0.0;
		lapAvgErr[iNest-1] = 0.0;
		double interpDenom = 0.0;
		double gradDenom = 0.0;
		double lapDenom = 0.0;
		double sMax = 0.0;
		double gMax = 0.0;
		double lMax = 0.0;
		for ( int i = 0; i < pMesh.nParticles(); ++i)
		{
			interpError.scalar[i] = interpolatedScalar.scalar[i] - gaussScalar.scalar[i]; 
			double errI = std::abs(interpError.scalar[i]);
			interpAvgErr[iNest-1] += errI * errI * pp->area[i];
			interpDenom += gaussScalar.scalar[i] * gaussScalar.scalar[i]  * pp->area[i];
			if ( gaussScalar.scalar[i] > sMax )
				sMax = gaussScalar.scalar[i];
			if (  errI > interpMaxErr[iNest-1] )
				interpMaxErr[iNest-1] = errI;
			
			gradError.xComp[i] = interpolatedGradient.xComp[i] - exactGradient.xComp[i];
			gradError.yComp[i] = interpolatedGradient.yComp[i] - exactGradient.yComp[i];
			double errG = std::sqrt( gradError.xComp[i]*gradError.xComp[i] + gradError.yComp[i]*gradError.yComp[i] );
			gradAvgErr[iNest-1] += errG * errG * pp->area[i];
			gradDenom += ( exactGradient.xComp[i]*exactGradient.xComp[i] + exactGradient.yComp[i]*exactGradient.yComp[i]) * pp->area[i];
			if ( std::sqrt( exactGradient.xComp[i]*exactGradient.xComp[i] + exactGradient.yComp[i]*exactGradient.yComp[i]) > gMax )
				gMax = std::sqrt( exactGradient.xComp[i]*exactGradient.xComp[i] + exactGradient.yComp[i]*exactGradient.yComp[i]);
			if ( errG > gradMaxErr[iNest-1] )
				gradMaxErr[iNest-1] = errG;
			
			lapError.scalar[i] = interpolatedLaplacian.scalar[i] - exactLaplacian.scalar[i];
			double errL = std::abs( lapError.scalar[i] );
			lapAvgErr[iNest-1] += errL * errL * pp->area[i];
			lapDenom += exactLaplacian.scalar[i] * exactLaplacian.scalar[i] * pp->area[i];
			if ( exactLaplacian.scalar[i] > lMax )
				lMax = exactLaplacian.scalar[i];
			if ( errL > lapMaxErr[iNest-1] )
				lapMaxErr[iNest-1] = errL;
		}
		interpMaxErr[iNest-1] /= sMax;
		gradMaxErr[iNest-1] /= gMax;
		lapMaxErr[iNest-1] /= lMax;
		interpAvgErr[iNest-1] /= interpDenom;
		gradAvgErr[iNest-1] /= gradDenom;
		lapAvgErr[iNest-1] /= lapDenom;
		
		interpAvgErr[iNest-1] = std::sqrt( interpAvgErr[iNest-1] );
		gradAvgErr[iNest-1] = std::sqrt( gradAvgErr[iNest-1] );
		lapAvgErr[iNest-1] = std::sqrt( lapAvgErr[iNest-1] );
		
		std::cout << "error complete ..." << std::endl;
		
		
		
		
	}
	
	std::cout << "--- 1/N --- " << " --- interp max err ---- " << " --- interp avg err ---- " 
								<< " --- grad max err --- " << " --- grad avg err --- "
								<< " --- lap max err --- " << " --- lap avg err --- " << std::endl;
	for ( int i = 0; i < 7; ++i )
		std::cout << nInv[i] << " , "
				  << interpMaxErr[i] << " , "
				  << interpAvgErr[i] << " , "
				  << gradMaxErr[i] << " , "
				  << gradAvgErr[i] << " , "
				  << lapMaxErr[i] << " , "
				  << lapAvgErr[i] << std::endl;


return 0;
};