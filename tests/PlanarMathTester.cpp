/** @file PlanarMathTester.cpp
	@author Peter Bosler, Sandia National Laboratories
	@brief Unit test for PolyMesh2d, Field, and TriCubicHermite objects and implementations.
*/
#include "PolyMesh2d.h"
#include "xyzVector.h"
#include "Field.h"
#include "OutputMessage.h"
#include "Logger.h"
#include "TriCubicHermite.h"
#include <iostream>
#include <cmath>
#include <fstream>


int main ( int argc, char* argv[] )
{
	const int procRank = 0;
	const int numProcs = 1;
	const OutputMessage::priority loggingLevel = OutputMessage::debugPriority;
	
	Logger* exeLog = Logger::Instance(loggingLevel, procRank, numProcs);
	
	const int initNest = 5;
	const double maxR = 1.0;
	
	//
	// construct a mesh to use for testing
	//
	PolyMesh2d pMesh( initNest, PolyMesh2d::triHexSeed, maxR, procRank, numProcs);
	LongMessage meshStatus("planar mesh built", OutputMessage::tracePriority, "main", pMesh.getInfo() );
	//exeLog->logMessage(meshStatus);
	
	//
	// construct scalar fields to use for testing
	//
	const int scalarDim = 1;
	const double scalarVal = 2.0;
	const int vectorDim = 2;
	Field constScalar("constant", "n/a", scalarDim, pMesh.nParticles() );
	constScalar.initializeToScalarConstant( scalarVal );
	std::cout << "nParticles = " << pMesh.nParticles() << ", field.N = " << constScalar.N() << ". \n";
	constScalar.outputForMatlab("constScalarField.m", *pMesh.getParticles() );
	
	Field gaussScalar("Gaussian", "n/a", scalarDim, pMesh.nParticles() );
	Particles* pp = pMesh.getParticles();
	const double b = 8.0;
	for ( int i = 0; i < pMesh.nParticles(); ++i )
	{
		double x = pp->x[i];
		double y = pp->y[i];
		gaussScalar.insertScalarToField( i, exp( - b * b * ( x * x + y * y ) ) );
	}
	gaussScalar.outputForMatlab("gaussScalarField.m", *pMesh.getParticles() );
	
	//
	// initialize interpolation
	//
	TriCubicHermite constScalarInterp( &pMesh, &constScalar );
	//Field zeroVector("zero vector", "n/a", vectorDim, pMesh.nParticles() );
	Field zeroVector = constScalarInterp.estimateScalarGradientAtVertices( 0, 0);
	constScalarInterp.findAllCoefficientsScalar( zeroVector );

	TriCubicHermite gaussScalarInterp( &pMesh, &gaussScalar );
	Field gaussGrad = gaussScalarInterp.estimateScalarGradientAtVertices( 2 , 1 );
	gaussScalarInterp.findAllCoefficientsScalar( gaussGrad );
	
	//
	// interpolate from particles / planar mesh to uniform mesh
	//
	const int nUnif = 80;
	const double dx = 2.0 * maxR / (nUnif-1);
	double unifExact[nUnif][nUnif];
	double unifInterp[nUnif][nUnif];
	double unifLap[nUnif][nUnif];
	double unifGradX[nUnif][nUnif];
	double unifGradY[nUnif][nUnif];
	
	double unifExact2[nUnif][nUnif];
	double unifInterp2[nUnif][nUnif];
	double unifLap2[nUnif][nUnif];
	double unifGradX2[nUnif][nUnif];
	double unifGradY2[nUnif][nUnif];
	
	for ( int i = 0; i < nUnif; ++i )
	{
		double yi = -maxR + i * dx;
		for ( int j = 0; j < nUnif; ++j )
		{
			double xj = -maxR + j * dx;
			
			unifExact[i][j] = scalarVal;
			unifExact2[i][j] = exp( - b * b * ( xj * xj + yi * yi ));
			
			xyzVector loc(xj,yi);
			unifInterp[i][j] = constScalarInterp.interpolateScalar( loc );			
			unifLap[i][j] = constScalarInterp.scalarLaplacian( loc );
			
			unifInterp2[i][j] = gaussScalarInterp.interpolateScalar( loc );			
			unifLap2[i][j] = gaussScalarInterp.scalarLaplacian( loc );
			
			xyzVector grad = constScalarInterp.scalarGradient( loc );
			
			xyzVector grad2 = gaussScalarInterp.scalarGradient( loc );
			
			unifGradX[i][j] = grad.x;
			unifGradY[i][j] = grad.y;
			
			unifGradX2[i][j] = grad2.x;
			unifGradY2[i][j] = grad2.y;
		}
	}

	std::string constTestFilename = "constTests.m";
	std::ofstream constTestFile( constTestFilename.c_str() );
	constTestFile << "unifX = [";
	for ( int i = 0; i < nUnif - 1; ++i )
		constTestFile << -maxR + i * dx << ", ";
	constTestFile << -maxR + (nUnif-1)*dx << "];" << std::endl;
	constTestFile << "unifY = unifX;" << std::endl << std::endl;
	
	constTestFile << "unifExact = [";
	for ( int i = 0; i < nUnif-1; ++i )
	{
		for ( int j = 0; j < nUnif-1; ++j )
		{
			constTestFile << unifExact[i][j] << ", ";
		}
		constTestFile << unifExact[i][nUnif - 1] << "; ";
	}
	for ( int j = 0; j < nUnif-1; ++j )
		constTestFile << unifExact[nUnif-1][j] << ", ";
	constTestFile << unifExact[nUnif-1][nUnif-1] << "];" << std::endl << std::endl;
	
	constTestFile << "unifExact2 = [";
	for ( int i = 0; i < nUnif-1; ++i )
	{
		for ( int j = 0; j < nUnif-1; ++j )
		{
			constTestFile << unifExact2[i][j] << ", ";
		}
		constTestFile << unifExact2[i][nUnif - 1] << "; ";
	}
	for ( int j = 0; j < nUnif-1; ++j )
		constTestFile << unifExact2[nUnif-1][j] << ", ";
	constTestFile << unifExact2[nUnif-1][nUnif-1] << "];" << std::endl << std::endl;
	
	constTestFile << "unifInterp = [";
	for ( int i = 0; i < nUnif-1; ++i )
	{
		for ( int j = 0; j < nUnif-1; ++j )
		{
			constTestFile << unifInterp[i][j] << ", ";
		}
		constTestFile << unifInterp[i][nUnif - 1] << "; ";
	}
	for ( int j = 0; j < nUnif-1; ++j )
		constTestFile << unifInterp[nUnif-1][j] << ", ";
	constTestFile << unifInterp[nUnif-1][nUnif-1] << "];" << std::endl << std::endl;
	
	constTestFile << "unifInterp2 = [";
	for ( int i = 0; i < nUnif-1; ++i )
	{
		for ( int j = 0; j < nUnif-1; ++j )
		{
			constTestFile << unifInterp2[i][j] << ", ";
		}
		constTestFile << unifInterp2[i][nUnif - 1] << "; ";
	}
	for ( int j = 0; j < nUnif-1; ++j )
		constTestFile << unifInterp2[nUnif-1][j] << ", ";
	constTestFile << unifInterp2[nUnif-1][nUnif-1] << "];" << std::endl << std::endl;
	
	constTestFile << "unifLap = [";
	for ( int i = 0; i < nUnif-1; ++i )
	{
		for ( int j = 0; j < nUnif-1; ++j )
		{
			constTestFile << unifLap[i][j] << ", ";
		}
		constTestFile << unifLap[i][nUnif - 1] << "; ";
	}
	for ( int j = 0; j < nUnif-1; ++j )
		constTestFile << unifLap[nUnif-1][j] << ", ";
	constTestFile << unifLap[nUnif-1][nUnif-1] << "];" << std::endl << std::endl;
	
	constTestFile << "unifLap2 = [";
	for ( int i = 0; i < nUnif-1; ++i )
	{
		for ( int j = 0; j < nUnif-1; ++j )
		{
			constTestFile << unifLap2[i][j] << ", ";
		}
		constTestFile << unifLap2[i][nUnif - 1] << "; ";
	}
	for ( int j = 0; j < nUnif-1; ++j )
		constTestFile << unifLap2[nUnif-1][j] << ", ";
	constTestFile << unifLap2[nUnif-1][nUnif-1] << "];" << std::endl << std::endl;
	
	constTestFile << "unifGradX = [";
	for ( int i = 0; i < nUnif-1; ++i )
	{
		for ( int j = 0; j < nUnif-1; ++j )
		{
			constTestFile << unifGradX[i][j] << ", ";
		}
		constTestFile << unifGradX[i][nUnif - 1] << "; ";
	}
	for ( int j = 0; j < nUnif-1; ++j )
		constTestFile << unifGradX[nUnif-1][j] << ", ";
	constTestFile << unifGradX[nUnif-1][nUnif-1] << "];" << std::endl << std::endl;
	
	constTestFile << "unifGradY = [";
	for ( int i = 0; i < nUnif-1; ++i )
	{
		for ( int j = 0; j < nUnif-1; ++j )
		{
			constTestFile << unifGradY[i][j] << ", ";
		}
		constTestFile << unifGradY[i][nUnif - 1] << "; ";
	}
	for ( int j = 0; j < nUnif-1; ++j )
		constTestFile << unifGradY[nUnif-1][j] << ", ";
	constTestFile << unifGradY[nUnif-1][nUnif-1] << "];" << std::endl << std::endl;
	
	constTestFile << "unifGradX2 = [";
	for ( int i = 0; i < nUnif-1; ++i )
	{
		for ( int j = 0; j < nUnif-1; ++j )
		{
			constTestFile << unifGradX2[i][j] << ", ";
		}
		constTestFile << unifGradX2[i][nUnif - 1] << "; ";
	}
	for ( int j = 0; j < nUnif-1; ++j )
		constTestFile << unifGradX2[nUnif-1][j] << ", ";
	constTestFile << unifGradX2[nUnif-1][nUnif-1] << "];" << std::endl << std::endl;
	
	constTestFile << "unifGradY2 = [";
	for ( int i = 0; i < nUnif-1; ++i )
	{
		for ( int j = 0; j < nUnif-1; ++j )
		{
			constTestFile << unifGradY2[i][j] << ", ";
		}
		constTestFile << unifGradY2[i][nUnif - 1] << "; ";
	}
	for ( int j = 0; j < nUnif-1; ++j )
		constTestFile << unifGradY2[nUnif-1][j] << ", ";
	constTestFile << unifGradY2[nUnif-1][nUnif-1] << "];" << std::endl << std::endl;
return 0;
};