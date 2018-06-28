
#include "jello.h"
#include "physics.h"

const double restLen = 1.0 / 7.0;
const double restLenDiagonal2D = restLen * sqrt(2.0);
const double restLenDiagonal3D = sqrt(restLen * restLen + restLenDiagonal2D * restLenDiagonal2D);
const double restLenBend = 2.0 / 7.0;

// Hook Force
void hookForce(const point& a, const point& b, double kHook, point& f)
{
	// initialize temp force f
	f.x = 0;
	f.y = 0;
	f.z = 0;
	point vectorDiff;
	point normalizedvectorDiff;
	double vectorDiffLength;
	double elasticForceScalar;
	
	// let vectorDiff be the vector pointing from b to a
	pDIFFERENCE(a, b, vectorDiff);

	vectorDiffLength = sqrt(vectorDiff.x * vectorDiff.x + vectorDiff.y * vectorDiff.y + vectorDiff.z * vectorDiff.z);
	
	// normalized vector vectorDiff
	normalizedvectorDiff.x = vectorDiff.x / vectorDiffLength;
	normalizedvectorDiff.y = vectorDiff.y / vectorDiffLength;
	normalizedvectorDiff.z = vectorDiff.z / vectorDiffLength;
	
	// the elastic force exerted on a is
	elasticForceScalar = (-kHook) * (vectorDiffLength - restLen);
	
	pMULTIPLY(normalizedvectorDiff, elasticForceScalar, f);
}

// Hook force Diagonal
void hookForceDiagonal2D(const point& a, const point& b, double kHook, point& f)
{
	// initialize temp force f
	f.x = 0;
	f.y = 0;
	f.z = 0;

	point vectorDiff;
	point normalizedvectorDiff;
	double vectorDiffLength;
	double elasticForceScalar;

	// let vectorDiff be the vector pointing from b to a
	pDIFFERENCE(a, b, vectorDiff);
	vectorDiffLength = sqrt(vectorDiff.x * vectorDiff.x + vectorDiff.y * vectorDiff.y + vectorDiff.z * vectorDiff.z);

	// normalized vector vectorDiff
	normalizedvectorDiff.x = vectorDiff.x / vectorDiffLength;
	normalizedvectorDiff.y = vectorDiff.y / vectorDiffLength;
	normalizedvectorDiff.z = vectorDiff.z / vectorDiffLength;

	// the elastic force exerted on a is
	elasticForceScalar = (-kHook) * (vectorDiffLength - restLenDiagonal2D);

	pMULTIPLY(normalizedvectorDiff, elasticForceScalar, f);
}

// Hook Force Diagonal 3D
void hookForceDiagonal3D(const point& a, const point& b, double kHook, point& f)
{
	// initialize temp force f
	f.x = 0;
	f.y = 0;
	f.z = 0;
	point vectorDiff;
	point normalizedvectorDiff;
	double vectorDiffLength;
	double elasticForceScalar;

	// let vectorDiff be the vector pointing from b to a
	pDIFFERENCE(a, b, vectorDiff);
	vectorDiffLength = sqrt(vectorDiff.x * vectorDiff.x + vectorDiff.y * vectorDiff.y + vectorDiff.z * vectorDiff.z);
	
	// normalized vector vectorDiff
	normalizedvectorDiff.x = vectorDiff.x / vectorDiffLength;
	normalizedvectorDiff.y = vectorDiff.y / vectorDiffLength;
	normalizedvectorDiff.z = vectorDiff.z / vectorDiffLength;
	
	// the elastic force exerted on a is
	elasticForceScalar = (-kHook) * (vectorDiffLength - restLenDiagonal3D);
	pMULTIPLY(normalizedvectorDiff, elasticForceScalar, f);
}

// Hook Force for Bending
void hookForceBend(const point& a, const point& b, double kHook, point& f)
{
	// initialize temp force f
	f.x = 0;
	f.y = 0;
	f.z = 0;
	point vectorDiff;
	point normalizedvectorDiff;
	double vectorDiffLength;
	double elasticForceScalar;
	// let vectorDiff be the vector pointing from b to a
	pDIFFERENCE(a, b, vectorDiff);
	vectorDiffLength = sqrt(vectorDiff.x * vectorDiff.x + vectorDiff.y * vectorDiff.y + vectorDiff.z * vectorDiff.z);
	// normalized vector vectorDiff
	normalizedvectorDiff.x = vectorDiff.x / vectorDiffLength;
	normalizedvectorDiff.y = vectorDiff.y / vectorDiffLength;
	normalizedvectorDiff.z = vectorDiff.z / vectorDiffLength;
	// the elastic force exerted on a is
	elasticForceScalar = (-kHook) * (vectorDiffLength - restLenBend);
	pMULTIPLY(normalizedvectorDiff, elasticForceScalar, f);
}

// Hook Force for Collission Spring (Penalty method)
void hookForceCollision(const point& a, const point& b, double kHook, point& n, point& f)
{
	// initialize temp force f
	f.x = 0;
	f.y = 0;
	f.z = 0;
	point vectorDiff;
	double elasticForceScalar;
	// let vectorDiff be the vector pointing from b to a
	pDIFFERENCE(a, b, vectorDiff);
	pMULTIPLY(vectorDiff, -kHook, vectorDiff);
	DOTPRODUCTp(vectorDiff, n, elasticForceScalar);
	pMULTIPLY(n, elasticForceScalar, f);
}

// Damping Force
void dampingForce(const point& a, const point& b, const point& va, const point& vb, double kDamp, point& f)
{
	// initialize temp force f
	f.x = 0;
	f.y = 0;
	f.z = 0;
	point vectorDiff;
	point vectorDiffVelocity;
	point normalizedvectorDiff;
	double vectorDiffLength;
	double elasticForceScalar;
	// let vectorDiff be the vector pointing from b to a
	pDIFFERENCE(a, b, vectorDiff);
	vectorDiffLength = sqrt(vectorDiff.x * vectorDiff.x + vectorDiff.y * vectorDiff.y + vectorDiff.z * vectorDiff.z);
	// here va and vb are velocities of points a and b
	pDIFFERENCE(va, vb, vectorDiffVelocity);
	elasticForceScalar = (-kDamp) * ((vectorDiffVelocity.x * vectorDiff.x + vectorDiffVelocity.y * vectorDiff.y + vectorDiffVelocity.z * vectorDiff.z) / vectorDiffLength);
	// normalized vector vectorDiff
	normalizedvectorDiff.x = vectorDiff.x / vectorDiffLength;
	normalizedvectorDiff.y = vectorDiff.y / vectorDiffLength;
	normalizedvectorDiff.z = vectorDiff.z / vectorDiffLength;
	pMULTIPLY(normalizedvectorDiff, elasticForceScalar, f);
}

// Damping force for collision
void dampingForceCollision(const point& a, const point& b, const point& va, const point& vb, double kDamp, point& n, point& f)
{
	// initialize temp force f
	f.x = 0;
	f.y = 0;
	f.z = 0;
	point vectorDiff;
	point vectorDiffVelocity;
	point normalizedvectorDiff;
	double vectorDiffLength;
	double elasticForceScalar;
	double cosF;
	// let vectorDiff be the vector pointing from b to a
	pDIFFERENCE(a, b, vectorDiff);
	vectorDiffLength = sqrt(vectorDiff.x * vectorDiff.x + vectorDiff.y * vectorDiff.y + vectorDiff.z * vectorDiff.z);
	// here va and vb are velocities of points a and b
	pDIFFERENCE(va, vb, vectorDiffVelocity);
	elasticForceScalar = (-kDamp) * ((vectorDiffVelocity.x * vectorDiff.x + vectorDiffVelocity.y * vectorDiff.y + vectorDiffVelocity.z * vectorDiff.z) / vectorDiffLength);
	// normalized vector vectorDiff
	normalizedvectorDiff.x = vectorDiff.x / vectorDiffLength;
	normalizedvectorDiff.y = vectorDiff.y / vectorDiffLength;
	normalizedvectorDiff.z = vectorDiff.z / vectorDiffLength;
	pMULTIPLY(normalizedvectorDiff, elasticForceScalar, f);
	DOTPRODUCTp(f, n, cosF);
	pMULTIPLY(n, cosF, f);
}

// Structural Force between point p at (i, j, k) and its neighbour point b
void structualForce(world *jello, int i, int j, int k, point& a)
{
	point f;
	if (i > 0) // point p at (i, j, k) has its left neighbour (i-1, j, k), calculate the hook force and damping force from left
	{
		hookForce(jello->p[i][j][k], jello->p[i - 1][j][k], jello->kElastic, f);
		// accumulate the hook force at point p
		pSUM(f, a, a);
		dampingForce(jello->p[i][j][k], jello->p[i - 1][j][k], jello->v[i][j][k], jello->v[i - 1][j][k], jello->dElastic, f);
		// accumulate the damping force at point p
		pSUM(f, a, a);
	}
	if (i < 7) // point p at (i, j, k) has its right neighbour (i+1, j, k), calculate the hook force and damping force from right
	{
		hookForce(jello->p[i][j][k], jello->p[i + 1][j][k], jello->kElastic, f);
		// accumulate the hook force at point p
		pSUM(f, a, a);
		dampingForce(jello->p[i][j][k], jello->p[i + 1][j][k], jello->v[i][j][k], jello->v[i + 1][j][k], jello->dElastic, f);
		// accumulate the damping force at point p
		pSUM(f, a, a);
	}
	if (j > 0) // point p at (i, j, k) has its upper neighbour (i, j-1, k), calculate the hook force and damping force from upper
	{
		hookForce(jello->p[i][j][k], jello->p[i][j - 1][k], jello->kElastic, f);
		// accumulate the hook force at point p 
		pSUM(f, a, a);
		dampingForce(jello->p[i][j][k], jello->p[i][j - 1][k], jello->v[i][j][k], jello->v[i][j - 1][k], jello->dElastic, f);
		// accumulate the damping force at point p
		pSUM(f, a, a);
	}
	if (j < 7) // point p at (i, j, k) has its bottom neighbour (i, j+1, k), calculate the hook force and damping force from bottom
	{
		hookForce(jello->p[i][j][k], jello->p[i][j + 1][k], jello->kElastic, f);
		// accumulate the hook force at point p
		pSUM(f, a, a);
		dampingForce(jello->p[i][j][k], jello->p[i][j + 1][k], jello->v[i][j][k], jello->v[i][j + 1][k], jello->dElastic, f);
		// accumulate the damping force at point p
		pSUM(f, a, a);
	}
	if (k > 0) // point p at (i, j, k) has its front neighbour (i, j, k-1), calculate the hook force and damping force from front
	{
		hookForce(jello->p[i][j][k], jello->p[i][j][k - 1], jello->kElastic, f);
		// accumulate the hook force at point p
		pSUM(f, a, a);
		dampingForce(jello->p[i][j][k], jello->p[i][j][k - 1], jello->v[i][j][k], jello->v[i][j][k - 1], jello->dElastic, f);
		// accumulate the damping force at point p
		pSUM(f, a, a);
	}
	if (k < 7) // point p at (i, j, k) has its behind neighbour (i, j, k+1), calculate the hook force and damping force from behind
	{
		hookForce(jello->p[i][j][k], jello->p[i][j][k + 1], jello->kElastic, f);
		// accumulate the hook force at point p
		pSUM(f, a, a);
		dampingForce(jello->p[i][j][k], jello->p[i][j][k + 1], jello->v[i][j][k], jello->v[i - 1][j][k], jello->dElastic, f);
		// accumulate the damping force at point p
		pSUM(f, a, a);
	}
}

// Shear Force between point p at (i, j, k) and its diagonal neighbour point b
void shearForce(world *jello, int i, int j, int k, point& a)
{
	point f;
	if (i > 0)
	{
		// on 2D surface
		if (k > 0) // point p at (i, j, k) has its neighbour (i-1, j, k-1), calculate the hook force and damping force
		{
			hookForceDiagonal2D(jello->p[i][j][k], jello->p[i - 1][j][k - 1], jello->kElastic, f);
			// accumulate the hook force at point p
			pSUM(f, a, a);
			dampingForce(jello->p[i][j][k], jello->p[i - 1][j][k - 1], jello->v[i][j][k], jello->v[i - 1][j][k - 1], jello->dElastic, f);
			// accumulate the damping force at point p
			pSUM(f, a, a);
		}
		if (k < 7) // point p at (i, j, k) has its neighbour (i-1, j, k+1), calculate the hook force and damping force
		{
			hookForceDiagonal2D(jello->p[i][j][k], jello->p[i - 1][j][k + 1], jello->kElastic, f);
			// accumulate the hook force at point p
			pSUM(f, a, a);
			dampingForce(jello->p[i][j][k], jello->p[i - 1][j][k + 1], jello->v[i][j][k], jello->v[i - 1][j][k + 1], jello->dElastic, f);
			// accumulate the damping force at point p
			pSUM(f, a, a);
		}
		if (j > 0) // point p at (i, j, k) has its neighbour (i-1, j-1, k), calculate the hook force and damping force
		{
			hookForceDiagonal2D(jello->p[i][j][k], jello->p[i - 1][j - 1][k], jello->kElastic, f);
			// accumulate the hook force at point p
			pSUM(f, a, a);
			dampingForce(jello->p[i][j][k], jello->p[i - 1][j - 1][k], jello->v[i][j][k], jello->v[i - 1][j - 1][k], jello->dElastic, f);
			// accumulate the damping force at point p
			pSUM(f, a, a);
		}
		if (j < 7) // point p at (i, j, k) has its neighbour (i-1, j+1, k), calculate the hook force and damping force
		{
			hookForceDiagonal2D(jello->p[i][j][k], jello->p[i - 1][j + 1][k], jello->kElastic, f);
			// accumulate the hook force at point p
			pSUM(f, a, a);
			dampingForce(jello->p[i][j][k], jello->p[i - 1][j + 1][k], jello->v[i][j][k], jello->v[i - 1][j + 1][k], jello->dElastic, f);
			// accumulate the damping force at point p
			pSUM(f, a, a);
		}
		// diagonals in 3D cube
		if (j > 0 && k > 0) // point p at (i, j, k) has its neighbour (i-1, j-1, k-1), calculate the hook force and damping force
		{
			hookForceDiagonal3D(jello->p[i][j][k], jello->p[i - 1][j - 1][k - 1], jello->kElastic, f);
			// accumulate the hook force at point p
			pSUM(f, a, a);
			dampingForce(jello->p[i][j][k], jello->p[i - 1][j - 1][k - 1], jello->v[i][j][k], jello->v[i - 1][j - 1][k - 1], jello->dElastic, f);
			// accumulate the damping force at point p
			pSUM(f, a, a);
		}
		if (j > 0 && k < 7) // point p at (i, j, k) has its neighbour (i-1, j-1, k+1), calculate the hook force and damping force
		{
			hookForceDiagonal3D(jello->p[i][j][k], jello->p[i - 1][j - 1][k + 1], jello->kElastic, f);
			// accumulate the hook force at point p
			pSUM(f, a, a);
			dampingForce(jello->p[i][j][k], jello->p[i - 1][j - 1][k + 1], jello->v[i][j][k], jello->v[i - 1][j - 1][k + 1], jello->dElastic, f);
			// accumulate the damping force at point p
			pSUM(f, a, a);
		}
		if (j < 7 && k > 0) // point p at (i, j, k) has its neighbour (i-1, j+1, k-1), calculate the hook force and damping force
		{
			hookForceDiagonal3D(jello->p[i][j][k], jello->p[i - 1][j + 1][k - 1], jello->kElastic, f);
			// accumulate the hook force at point p
			pSUM(f, a, a);
			dampingForce(jello->p[i][j][k], jello->p[i - 1][j + 1][k - 1], jello->v[i][j][k], jello->v[i - 1][j + 1][k - 1], jello->dElastic, f);
			// accumulate the damping force at point p
			pSUM(f, a, a);
		}
		if (j < 7 && k < 7) // point p at (i, j, k) has its neighbour (i-1, j+1, k+1), calculate the hook force and damping force
		{
			hookForceDiagonal3D(jello->p[i][j][k], jello->p[i - 1][j + 1][k + 1], jello->kElastic, f);
			// accumulate the hook force at point p
			pSUM(f, a, a);
			dampingForce(jello->p[i][j][k], jello->p[i - 1][j + 1][k + 1], jello->v[i][j][k], jello->v[i - 1][j + 1][k + 1], jello->dElastic, f);
			// accumulate the damping force at point p
			pSUM(f, a, a);
		}
	}
	if (i < 7)
	{
		// on 2D surface
		if (k > 0) // point p at (i, j, k) has its neighbour (i+1, j, k-1), calculate the hook force and damping force
		{
			hookForceDiagonal2D(jello->p[i][j][k], jello->p[i + 1][j][k - 1], jello->kElastic, f);
			// accumulate the hook force at point p
			pSUM(f, a, a);
			dampingForce(jello->p[i][j][k], jello->p[i + 1][j][k - 1], jello->v[i][j][k], jello->v[i + 1][j][k - 1], jello->dElastic, f);
			// accumulate the damping force at point p
			pSUM(f, a, a);
		}
		if (k < 7) // point p at (i, j, k) has its neighbour (i+1, j, k+1), calculate the hook force and damping force
		{
			hookForceDiagonal2D(jello->p[i][j][k], jello->p[i + 1][j][k + 1], jello->kElastic, f);
			// accumulate the hook force at point p
			pSUM(f, a, a);
			dampingForce(jello->p[i][j][k], jello->p[i + 1][j][k + 1], jello->v[i][j][k], jello->v[i + 1][j][k + 1], jello->dElastic, f);
			// accumulate the damping force at point p
			pSUM(f, a, a);
		}
		if (j > 0) // point p at (i, j, k) has its neighbour (i+1, j-1, k), calculate the hook force and damping force
		{
			hookForceDiagonal2D(jello->p[i][j][k], jello->p[i + 1][j - 1][k], jello->kElastic, f);
			// accumulate the hook force at point p
			pSUM(f, a, a);
			dampingForce(jello->p[i][j][k], jello->p[i + 1][j - 1][k], jello->v[i][j][k], jello->v[i + 1][j - 1][k], jello->dElastic, f);
			// accumulate the damping force at point p
			pSUM(f, a, a);
		}
		if (j < 7) // point p at (i, j, k) has its neighbour (i+1, j+1, k), calculate the hook force and damping force
		{
			hookForceDiagonal2D(jello->p[i][j][k], jello->p[i + 1][j + 1][k], jello->kElastic, f);
			// accumulate the hook force at point p
			pSUM(f, a, a);
			dampingForce(jello->p[i][j][k], jello->p[i + 1][j + 1][k], jello->v[i][j][k], jello->v[i + 1][j + 1][k], jello->dElastic, f);
			// accumulate the damping force at point p
			pSUM(f, a, a);
		}
		// diagonals in 3D cube
		if (j > 0 && k > 0) // point p at (i, j, k) has its neighbour (i+1, j-1, k-1), calculate the hook force and damping force
		{
			hookForceDiagonal3D(jello->p[i][j][k], jello->p[i + 1][j - 1][k - 1], jello->kElastic, f);
			// accumulate the hook force at point p
			pSUM(f, a, a);
			dampingForce(jello->p[i][j][k], jello->p[i + 1][j - 1][k - 1], jello->v[i][j][k], jello->v[i + 1][j - 1][k - 1], jello->dElastic, f);
			// accumulate the damping force at point p
			pSUM(f, a, a);
		}
		if (j > 0 && k < 7) // point p at (i, j, k) has its neighbour (i+1, j-1, k+1), calculate the hook force and damping force
		{
			hookForceDiagonal3D(jello->p[i][j][k], jello->p[i + 1][j - 1][k + 1], jello->kElastic, f);
			// accumulate the hook force at point p
			pSUM(f, a, a);
			dampingForce(jello->p[i][j][k], jello->p[i + 1][j - 1][k + 1], jello->v[i][j][k], jello->v[i + 1][j - 1][k + 1], jello->dElastic, f);
			// accumulate the damping force at point p
			pSUM(f, a, a);
		}
		if (j < 7 && k > 0) // point p at (i, j, k) has its neighbour (i+1, j+1, k-1), calculate the hook force and damping force
		{
			hookForceDiagonal3D(jello->p[i][j][k], jello->p[i + 1][j + 1][k - 1], jello->kElastic, f);
			// accumulate the hook force at point p
			pSUM(f, a, a);
			dampingForce(jello->p[i][j][k], jello->p[i + 1][j + 1][k - 1], jello->v[i][j][k], jello->v[i + 1][j + 1][k - 1], jello->dElastic, f);
			// accumulate the damping force at point p
			pSUM(f, a, a);
		}
		if (j < 7 && k < 7) // point p at (i, j, k) has its neighbour (i+1, j+1, k+1), calculate the hook force and damping force
		{
			hookForceDiagonal3D(jello->p[i][j][k], jello->p[i + 1][j + 1][k + 1], jello->kElastic, f);
			// accumulate the hook force at point p
			pSUM(f, a, a);
			dampingForce(jello->p[i][j][k], jello->p[i + 1][j + 1][k + 1], jello->v[i][j][k], jello->v[i + 1][j + 1][k + 1], jello->dElastic, f);
			// accumulate the damping force at point p
			pSUM(f, a, a);
		}
	}
	if (j > 0)
	{
		// on 2D surface
		if (k > 0) // point p at (i, j, k) has its neighbour (i, j-1, k-1), calculate the hook force and damping force
		{
			hookForceDiagonal2D(jello->p[i][j][k], jello->p[i][j - 1][k - 1], jello->kElastic, f);
			// accumulate the hook force at point p
			pSUM(f, a, a);
			dampingForce(jello->p[i][j][k], jello->p[i][j - 1][k - 1], jello->v[i][j][k], jello->v[i][j - 1][k - 1], jello->dElastic, f);
			// accumulate the damping force at point p
			pSUM(f, a, a);
		}
		if (k < 7) // point p at (i, j, k) has its neighbour (i, j-1, k+1), calculate the hook force and damping force
		{
			hookForceDiagonal2D(jello->p[i][j][k], jello->p[i][j - 1][k + 1], jello->kElastic, f);
			// accumulate the hook force at point p
			pSUM(f, a, a);
			dampingForce(jello->p[i][j][k], jello->p[i][j - 1][k + 1], jello->v[i][j][k], jello->v[i][j - 1][k + 1], jello->dElastic, f);
			// accumulate the damping force at point p
			pSUM(f, a, a);
		}
	}
	if (j < 7)
	{
		// on 2D surface
		if (k > 0) // point p at (i, j, k) has its neighbour (i, j+1, k-1), calculate the hook force and damping force
		{
			hookForceDiagonal2D(jello->p[i][j][k], jello->p[i][j + 1][k - 1], jello->kElastic, f);
			// accumulate the hook force at point p
			pSUM(f, a, a);
			dampingForce(jello->p[i][j][k], jello->p[i][j + 1][k - 1], jello->v[i][j][k], jello->v[i][j + 1][k - 1], jello->dElastic, f);
			// accumulate the damping force at point p
			pSUM(f, a, a);
		}
		if (k < 7) // point p at (i, j, k) has its neighbour (i, j+1, k+1), calculate the hook force and damping force
		{
			hookForceDiagonal2D(jello->p[i][j][k], jello->p[i][j + 1][k + 1], jello->kElastic, f);
			// accumulate the hook force at point p
			pSUM(f, a, a);
			dampingForce(jello->p[i][j][k], jello->p[i][j + 1][k + 1], jello->v[i][j][k], jello->v[i][j + 1][k + 1], jello->dElastic, f);
			// accumulate the damping force at point p
			pSUM(f, a, a);
		}
	}
}

// Bending Force between point p at (i, j, k) and its following neighbour point b
void bendForce(world *jello, int i, int j, int k, point& a)
{
	point f;
	if (i > 1) // point p at (i, j, k) has its following neighbour (i-2, j, k), calculate the hook force and damping force
	{
		hookForceBend(jello->p[i][j][k], jello->p[i - 2][j][k], jello->kElastic, f);
		// accumulate the hook force at point p
		pSUM(f, a, a);
		dampingForce(jello->p[i][j][k], jello->p[i - 2][j][k], jello->v[i][j][k], jello->v[i - 2][j][k], jello->dElastic, f);
		// accumulate the damping force at point p
		pSUM(f, a, a);
	}
	if (i < 6) // point p at (i, j, k) has its following neighbour (i+2, j, k), calculate the hook force and damping force
	{
		hookForceBend(jello->p[i][j][k], jello->p[i + 2][j][k], jello->kElastic, f);
		// accumulate the hook force at point p
		pSUM(f, a, a);
		dampingForce(jello->p[i][j][k], jello->p[i + 2][j][k], jello->v[i][j][k], jello->v[i + 2][j][k], jello->dElastic, f);
		// accumulate the damping force at point p
		pSUM(f, a, a);
	}
	if (j > 1) // point p at (i, j, k) has its following neighbour (i, j-2, k), calculate the hook force and damping force
	{
		hookForceBend(jello->p[i][j][k], jello->p[i][j - 2][k], jello->kElastic, f);
		// accumulate the hook force at point p
		pSUM(f, a, a);
		dampingForce(jello->p[i][j][k], jello->p[i][j - 2][k], jello->v[i][j][k], jello->v[i][j - 2][k], jello->dElastic, f);
		// accumulate the damping force at point p
		pSUM(f, a, a);
	}
	if (j < 6) // point p at (i, j, k) has its following neighbour (i, j+2, k), calculate the hook force and damping force
	{
		hookForceBend(jello->p[i][j][k], jello->p[i][j + 2][k], jello->kElastic, f);
		// accumulate the hook force at point p
		pSUM(f, a, a);
		dampingForce(jello->p[i][j][k], jello->p[i][j + 2][k], jello->v[i][j][k], jello->v[i][j + 2][k], jello->dElastic, f);
		// accumulate the damping force at point p
		pSUM(f, a, a);
	}
	if (k > 1) // point p at (i, j, k) has its following neighbour (i, j, k-2), calculate the hook force and damping force
	{
		hookForceBend(jello->p[i][j][k], jello->p[i][j][k - 2], jello->kElastic, f);
		// accumulate the hook force at point p
		pSUM(f, a, a);
		dampingForce(jello->p[i][j][k], jello->p[i][j][k - 2], jello->v[i][j][k], jello->v[i][j][k - 2], jello->dElastic, f);
		// accumulate the damping force at point p
		pSUM(f, a, a);
	}
	if (k < 6) // point p at (i, j, k) has its following neighbour (i, j, k+2), calculate the hook force and damping force
	{
		hookForceBend(jello->p[i][j][k], jello->p[i][j][k + 2], jello->kElastic, f);
		// accumulate the hook force at point p
		pSUM(f, a, a);
		dampingForce(jello->p[i][j][k], jello->p[i][j][k + 2], jello->v[i][j][k], jello->v[i][j][k + 2], jello->dElastic, f);
		// accumulate the damping force at point p
		pSUM(f, a, a);
	}
}


// External Force field, add it to point p at (i, j, k) 
void externalForce(world *jello, int x, int y, int z, point& a)
{
	// the external force index in resolution array
	int i, j, k;

	// forces at 8 corners in a specific grid
	point f000, f001;
	point f010, f011;
	point f100, f101;
	point f110, f111;

	// the external force position in grid
	double px, py, pz;

	// the force field grid
	double grid;

	// the external force field value
	point externalForce;
	externalForce.x = 0;
	externalForce.y = 0;
	externalForce.z = 0;

	// the array resolution represents the force field, the external force point position is 
	// ((-2 + i * 4 / (jello->resolution - 1)), (-2 + j * 4 / (jello->resolution - 1)), (-2 + k * 4 / (jello->resolution - 1))) by i, j, k in the bounding box
	// since (-2 + i * 4 / (jello->resolution - 1) <= pt.x <= (-2 + (i + 1) * 4 / (jello->resolution - 1), we have
	i = int((jello->p[x][y][z].x + 2) * (jello->resolution - 1) / 4);
	j = int((jello->p[x][y][z].y + 2) * (jello->resolution - 1) / 4);
	k = int((jello->p[x][y][z].z + 2) * (jello->resolution - 1) / 4);

	// check if the index is at the wall of the bounding box
	if (i == (jello->resolution - 1))
	{
		i--;
	}
	if (j == (jello->resolution - 1))
	{
		j--;
	}
	if (k == (jello->resolution - 1))
	{
		k--;
	}
	// check if the point is inside the bounding box, read the force field value
	if (((i >= 0) && (i <= jello->resolution - 1)) && ((j >= 0) && (j <= jello->resolution - 1)) && ((j >= 0) && (j <= jello->resolution - 1)))
	{
		f000 = jello->forceField[(i * jello->resolution * jello->resolution + j * jello->resolution + k)];
		f001 = jello->forceField[(i * jello->resolution * jello->resolution + j * jello->resolution + (k + 1))];

		f010 = jello->forceField[(i * jello->resolution * jello->resolution + (j + 1) * jello->resolution + k)];
		f011 = jello->forceField[(i * jello->resolution * jello->resolution + (j + 1) * jello->resolution + (k + 1))];

		f100 = jello->forceField[((i + 1) * jello->resolution * jello->resolution + j * jello->resolution + k)];
		f101 = jello->forceField[((i + 1) * jello->resolution * jello->resolution + j * jello->resolution + (k + 1))];

		f110 = jello->forceField[((i + 1) * jello->resolution * jello->resolution + (j + 1) * jello->resolution + k)];
		f111 = jello->forceField[((i + 1) * jello->resolution * jello->resolution + (j + 1) * jello->resolution + (k + 1))];

		// 3D interpolation
		grid = 1.0 * 4 / (jello->resolution - 1);
		px = (jello->p[x][y][z].x - (-2 + 1.0 * 4 * i / (jello->resolution - 1))) / grid;
		py = (jello->p[x][y][z].y - (-2 + 1.0 * 4 * j / (jello->resolution - 1))) / grid;
		pz = (jello->p[x][y][z].z - (-2 + 1.0 * 4 * k / (jello->resolution - 1))) / grid;

		pMULTIPLY(f000, (1 - px) * (1 - py) * (1 - pz), f000);
		pMULTIPLY(f001, (1 - px) * (1 - py) * pz, f001);
		pMULTIPLY(f010, (1 - px) * py * (1 - pz), f010);
		pMULTIPLY(f011, (1 - px) * py * pz, f011);
		pMULTIPLY(f100, px * (1 - py) * (1 - pz), f100);
		pMULTIPLY(f101, px * (1 - py) * pz, f101);
		pMULTIPLY(f110, px * py * (1 - pz), f110);
		pMULTIPLY(f111, px * py * pz, f111);

		pSUM(externalForce, f000, externalForce);
		pSUM(externalForce, f001, externalForce);
		pSUM(externalForce, f010, externalForce);
		pSUM(externalForce, f011, externalForce);
		pSUM(externalForce, f100, externalForce);
		pSUM(externalForce, f101, externalForce);
		pSUM(externalForce, f110, externalForce);
		pSUM(externalForce, f111, externalForce);
		a.x = externalForce.x;
		a.y = externalForce.y;
		a.z = externalForce.z;
	}
}

// Detect collision force by the bounding box walls
void collisionForce(world *jello, int i, int j, int k, point& a)
{
	// collision detection, check if point (i, j, k) is outside the bounding box
	if ((jello->p[i][j][k].x <= -2) || (jello->p[i][j][k].x >= 2) ||
		(jello->p[i][j][k].y <= -2) || (jello->p[i][j][k].y >= 2) ||
		(jello->p[i][j][k].z <= -2) || (jello->p[i][j][k].z >= 2))
	{
		point f;
		// initialize the normal vector of the wall
		point normal;
		normal.x = 0;
		normal.y = 0;
		normal.z = 0;
		// initialize the collision point on the wall
		point wall;
		wall = jello->p[i][j][k];
		// since the wall is static, the velocity is (0, 0, 0)
		point vWall;
		vWall.x = 0;
		vWall.y = 0;
		vWall.z = 0;
		// using the penalty method, calculate the force direction and the collision point on the wall
		if (jello->p[i][j][k].x <= -2)
		{
			normal.x = 1;
			wall.x = -2;
		}
		if (jello->p[i][j][k].x >= 2)
		{
			normal.x = -1;
			wall.x = 2;
		}
		if (jello->p[i][j][k].y <= -2)
		{
			normal.y = 1;
			wall.y = -2;
		}
		if (jello->p[i][j][k].y >= 2)
		{
			normal.y = -1;
			wall.y = 2;
		}
		if (jello->p[i][j][k].z <= -2)
		{
			normal.z = 1;
			wall.z = -2;
		}
		if (jello->p[i][j][k].z >= 2)
		{
			normal.z = -1;
			wall.z = 2;
		}
		hookForceCollision(jello->p[i][j][k], wall, jello->kCollision, normal, f);
		pSUM(f, a, a);
		dampingForceCollision(jello->p[i][j][k], wall, jello->v[i][j][k], vWall, jello->dCollision, normal, f);
		pSUM(f, a, a);
	}
}


/* Computes acceleration to every control point of the jello cube,
which is in state given by 'jello'.
Returns result in array 'a'. */
void computeAcceleration(struct world * jello, struct point a[8][8][8])
{	
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			for (int k = 0; k < 8; k++)
			{
				a[i][j][k].x = 0;
				a[i][j][k].y = 0;
				a[i][j][k].z = 0;
				
				// calculate the internal force of the jello
				structualForce(jello, i, j, k, a[i][j][k]);
				shearForce(jello, i, j, k, a[i][j][k]);
				bendForce(jello, i, j, k, a[i][j][k]);

			
				// calculate the external force in the bounding box
				point externF;
				externalForce(jello, i, j, k, externF);
				pSUM(a[i][j][k], externF, a[i][j][k]);
				
				// calculate the collision force in the bounding box
				
				collisionForce(jello, i, j, k, a[i][j][k]);

				// calculate the acceleration using F = m * a, a = F / m
				
				pMULTIPLY(a[i][j][k], (1 / jello->mass), a[i][j][k]);
			}
		}
	}
}

/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world * jello)
{
	int i, j, k;
	point a[8][8][8];

	computeAcceleration(jello, a);

	for (i = 0; i <= 7; i++)
		for (j = 0; j <= 7; j++)
			for (k = 0; k <= 7; k++)
			{
				jello->p[i][j][k].x += jello->dt * jello->v[i][j][k].x;
				jello->p[i][j][k].y += jello->dt * jello->v[i][j][k].y;
				jello->p[i][j][k].z += jello->dt * jello->v[i][j][k].z;
				jello->v[i][j][k].x += jello->dt * a[i][j][k].x;
				jello->v[i][j][k].y += jello->dt * a[i][j][k].y;
				jello->v[i][j][k].z += jello->dt * a[i][j][k].z;

			}
}

/* performs one step of RK4 Integration */
/* as a result, updates the jello structure */
void RK4(struct world * jello)
{
	point F1p[8][8][8], F1v[8][8][8],
		F2p[8][8][8], F2v[8][8][8],
		F3p[8][8][8], F3v[8][8][8],
		F4p[8][8][8], F4v[8][8][8];

	point a[8][8][8];


	struct world buffer;

	int i, j, k;

	buffer = *jello; // make a copy of jello

	computeAcceleration(jello, a);

	for (i = 0; i <= 7; i++)
		for (j = 0; j <= 7; j++)
			for (k = 0; k <= 7; k++)
			{
				pMULTIPLY(jello->v[i][j][k], jello->dt, F1p[i][j][k]);
				pMULTIPLY(a[i][j][k], jello->dt, F1v[i][j][k]);
				pMULTIPLY(F1p[i][j][k], 0.5, buffer.p[i][j][k]);
				pMULTIPLY(F1v[i][j][k], 0.5, buffer.v[i][j][k]);
				pSUM(jello->p[i][j][k], buffer.p[i][j][k], buffer.p[i][j][k]);
				pSUM(jello->v[i][j][k], buffer.v[i][j][k], buffer.v[i][j][k]);
			}

	computeAcceleration(&buffer, a);

	for (i = 0; i <= 7; i++)
		for (j = 0; j <= 7; j++)
			for (k = 0; k <= 7; k++)
			{
				// F2p = dt * buffer.v;
				pMULTIPLY(buffer.v[i][j][k], jello->dt, F2p[i][j][k]);
				// F2v = dt * a(buffer.p,buffer.v);     
				pMULTIPLY(a[i][j][k], jello->dt, F2v[i][j][k]);
				pMULTIPLY(F2p[i][j][k], 0.5, buffer.p[i][j][k]);
				pMULTIPLY(F2v[i][j][k], 0.5, buffer.v[i][j][k]);
				pSUM(jello->p[i][j][k], buffer.p[i][j][k], buffer.p[i][j][k]);
				pSUM(jello->v[i][j][k], buffer.v[i][j][k], buffer.v[i][j][k]);
			}

	computeAcceleration(&buffer, a);

	for (i = 0; i <= 7; i++)
		for (j = 0; j <= 7; j++)
			for (k = 0; k <= 7; k++)
			{
				// F3p = dt * buffer.v;
				pMULTIPLY(buffer.v[i][j][k], jello->dt, F3p[i][j][k]);
				// F3v = dt * a(buffer.p,buffer.v);     
				pMULTIPLY(a[i][j][k], jello->dt, F3v[i][j][k]);
				pMULTIPLY(F3p[i][j][k], 0.5, buffer.p[i][j][k]);
				pMULTIPLY(F3v[i][j][k], 0.5, buffer.v[i][j][k]);
				pSUM(jello->p[i][j][k], buffer.p[i][j][k], buffer.p[i][j][k]);
				pSUM(jello->v[i][j][k], buffer.v[i][j][k], buffer.v[i][j][k]);
			}

	computeAcceleration(&buffer, a);


	for (i = 0; i <= 7; i++)
		for (j = 0; j <= 7; j++)
			for (k = 0; k <= 7; k++)
			{
				// F3p = dt * buffer.v;
				pMULTIPLY(buffer.v[i][j][k], jello->dt, F4p[i][j][k]);
				// F3v = dt * a(buffer.p,buffer.v);     
				pMULTIPLY(a[i][j][k], jello->dt, F4v[i][j][k]);

				pMULTIPLY(F2p[i][j][k], 2, buffer.p[i][j][k]);
				pMULTIPLY(F3p[i][j][k], 2, buffer.v[i][j][k]);
				pSUM(buffer.p[i][j][k], buffer.v[i][j][k], buffer.p[i][j][k]);
				pSUM(buffer.p[i][j][k], F1p[i][j][k], buffer.p[i][j][k]);
				pSUM(buffer.p[i][j][k], F4p[i][j][k], buffer.p[i][j][k]);
				pMULTIPLY(buffer.p[i][j][k], 1.0 / 6, buffer.p[i][j][k]);
				pSUM(buffer.p[i][j][k], jello->p[i][j][k], jello->p[i][j][k]);

				pMULTIPLY(F2v[i][j][k], 2, buffer.p[i][j][k]);
				pMULTIPLY(F3v[i][j][k], 2, buffer.v[i][j][k]);
				pSUM(buffer.p[i][j][k], buffer.v[i][j][k], buffer.p[i][j][k]);
				pSUM(buffer.p[i][j][k], F1v[i][j][k], buffer.p[i][j][k]);
				pSUM(buffer.p[i][j][k], F4v[i][j][k], buffer.p[i][j][k]);
				pMULTIPLY(buffer.p[i][j][k], 1.0 / 6, buffer.p[i][j][k]);
				pSUM(buffer.p[i][j][k], jello->v[i][j][k], jello->v[i][j][k]);
			}

	return;
}
