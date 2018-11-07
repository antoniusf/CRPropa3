#ifndef CRPROPA_GRID_H
#define CRPROPA_GRID_H

#include "crpropa/Referenced.h"
#include "crpropa/Vector3.h"
#include <vector>
#include <Eigen/Core>

namespace crpropa {

/** Lower and upper neighbor in a periodically continued unit grid */
inline void periodicClamp(double x, int n, int &lo, int &hi) {
	lo = ((int(floor(x)) % n) + n) % n;
	hi = (lo + 1) % n;
}

/** Lower and upper neighbor in a reflectively repeated unit grid */
inline void reflectiveClamp(double x, int n, int &lo, int &hi) {
	while ((x < 0) or (x > n))
		x = 2 * n * (x > n) - x;
	lo = floor(x);
	hi = lo + (lo < n-1);
}

 inline int periodicBoundary(int index, int n) {
   return ((index % n) + n) % n;
 }

 inline int reflectiveBoundary(int index, int n) {
   return 0;
 }

/** Symmetrical round */
inline double round(double r) {
    return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
}

/**
 * \addtogroup Core
 * @{
 */
/**
 @class Grid
 @brief Template class for fields on a periodic grid with trilinear interpolation

 The grid spacing is constant and equal along all three axes.
 Values are calculated by trilinear interpolation of the surrounding 8 grid points.
 The grid is periodically (default) or reflectively extended.
 The grid sample positions are at 1/2 * size/N, 3/2 * size/N ... (2N-1)/2 * size/N.
 */
template<typename T>
class Grid: public Referenced {
	std::vector<T> grid;
	size_t Nx, Ny, Nz; /**< Number of grid points */
	Vector3d origin; /**< Origin of the volume that is represented by the grid. */
	Vector3d gridOrigin; /**< Grid origin */
	double spacing; /**< Distance between grid points, determines the extension of the grid */
	bool reflective; /**< If set to true, the grid is repeated reflectively instead of periodically */
        Eigen::Matrix<float, 64, 64> makePolynomeCoefficients; /**< This matrix is used to compute the actual interpolation polynome for each cell from the field values and their derivatives. */
	int component;

public:
	/** Constructor for cubic grid
	 @param	origin	Position of the lower left front corner of the volume
	 @param	N		Number of grid points in one direction
	 @param spacing	Spacing between grid points
	 */
	Grid(Vector3d origin, size_t N, double spacing) {
		setOrigin(origin);
		setGridSize(N, N, N);
		setSpacing(spacing);
		setReflective(false);
		component = 0;
		loadMakeInterpolationMatrix();
	}

	/** Constructor for non-cubic grid
	 @param	origin	Position of the lower left front corner of the volume
	 @param	Nx		Number of grid points in x-direction
	 @param	Ny		Number of grid points in y-direction
	 @param	Nz		Number of grid points in z-direction
	 @param spacing	Spacing between grid points
	 */
	Grid(Vector3d origin, size_t Nx, size_t Ny, size_t Nz, double spacing) {
		setOrigin(origin);
		setGridSize(Nx, Ny, Nz);
		setSpacing(spacing);
		setReflective(false);
		component = 0;
		loadMakeInterpolationMatrix();
	}

	void loadMakeInterpolationMatrix() {
	    const int temp[64][64] = 
		{{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{-3, 3, 0, 0, 0, 0, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{2, -2, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, -3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{9, -9, -9, 9, 0, 0, 0, 0, 6, 3, -6, -3, 0, 0, 0, 0, 6, -6, 3, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{-6, 6, 6, -6, 0, 0, 0, 0, -3, -3, 3, 3, 0, 0, 0, 0, -4, 4, -2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, -2, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{2, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{-6, 6, 6, -6, 0, 0, 0, 0, -4, -2, 4, 2, 0, 0, 0, 0, -3, 3, -3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, -1, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{4, -4, -4, 4, 0, 0, 0, 0, 2, 2, -2, -2, 0, 0, 0, 0, 2, -2, 2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, 0, 0, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, -1, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9, -9, -9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 3, -6, -3, 0, 0, 0, 0, 6, -6, 3, -3, 0, 0, 0, 0, 4, 2, 2, 1, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -6, 6, 6, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, -3, 3, 3, 0, 0, 0, 0, -4, 4, -2, 2, 0, 0, 0, 0, -2, -2, -1, -1, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -6, 6, 6, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, -2, 4, 2, 0, 0, 0, 0, -3, 3, -3, 3, 0, 0, 0, 0, -2, -1, -2, -1, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, -4, -4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, -2, -2, 0, 0, 0, 0, 2, -2, 2, -2, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0},
		{-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, -3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{9, -9, 0, 0, -9, 9, 0, 0, 6, 3, 0, 0, -6, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, -6, 0, 0, 3, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 2, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{-6, 6, 0, 0, 6, -6, 0, 0, -3, -3, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, 4, 0, 0, -2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, -2, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, -1, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9, -9, 0, 0, -9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 3, 0, 0, -6, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, -6, 0, 0, 3, -3, 0, 0, 4, 2, 0, 0, 2, 1, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -6, 6, 0, 0, 6, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, -3, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, 4, 0, 0, -2, 2, 0, 0, -2, -2, 0, 0, -1, -1, 0, 0},
		{9, 0, -9, 0, -9, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 3, 0, -6, 0, -3, 0, 6, 0, -6, 0, 3, 0, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 2, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 9, 0, -9, 0, -9, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 3, 0, -6, 0, -3, 0, 6, 0, -6, 0, 3, 0, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 2, 0, 2, 0, 1, 0},
		{-27, 27, 27, -27, 27, -27, -27, 27, -18, -9, 18, 9, 18, 9, -18, -9, -18, 18, -9, 9, 18, -18, 9, -9, -18, 18, 18, -18, -9, 9, 9, -9, -12, -6, -6, -3, 12, 6, 6, 3, -12, -6, 12, 6, -6, -3, 6, 3, -12, 12, -6, 6, -6, 6, -3, 3, -8, -4, -4, -2, -4, -2, -2, -1},
		{18, -18, -18, 18, -18, 18, 18, -18, 9, 9, -9, -9, -9, -9, 9, 9, 12, -12, 6, -6, -12, 12, -6, 6, 12, -12, -12, 12, 6, -6, -6, 6, 6, 6, 3, 3, -6, -6, -3, -3, 6, 6, -6, -6, 3, 3, -3, -3, 8, -8, 4, -4, 4, -4, 2, -2, 4, 4, 2, 2, 2, 2, 1, 1},
		{-6, 0, 6, 0, 6, 0, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 0, -3, 0, 3, 0, 3, 0, -4, 0, 4, 0, -2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, -2, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, -6, 0, 6, 0, 6, 0, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 0, -3, 0, 3, 0, 3, 0, -4, 0, 4, 0, -2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, -2, 0, -1, 0, -1, 0},
		{18, -18, -18, 18, -18, 18, 18, -18, 12, 6, -12, -6, -12, -6, 12, 6, 9, -9, 9, -9, -9, 9, -9, 9, 12, -12, -12, 12, 6, -6, -6, 6, 6, 3, 6, 3, -6, -3, -6, -3, 8, 4, -8, -4, 4, 2, -4, -2, 6, -6, 6, -6, 3, -3, 3, -3, 4, 2, 4, 2, 2, 1, 2, 1},
		{-12, 12, 12, -12, 12, -12, -12, 12, -6, -6, 6, 6, 6, 6, -6, -6, -6, 6, -6, 6, 6, -6, 6, -6, -8, 8, 8, -8, -4, 4, 4, -4, -3, -3, -3, -3, 3, 3, 3, 3, -4, -4, 4, 4, -2, -2, 2, 2, -4, 4, -4, 4, -2, 2, -2, 2, -2, -2, -2, -2, -1, -1, -1, -1},
		{2, 0, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{-6, 6, 0, 0, 6, -6, 0, 0, -4, -2, 0, 0, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, -3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, -1, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{4, -4, 0, 0, -4, 4, 0, 0, 2, 2, 0, 0, -2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -6, 6, 0, 0, 6, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, -2, 0, 0, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, -3, 3, 0, 0, -2, -1, 0, 0, -2, -1, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, -4, 0, 0, -4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, -2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 2, -2, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0},
		{-6, 0, 6, 0, 6, 0, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, 0, -2, 0, 4, 0, 2, 0, -3, 0, 3, 0, -3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, -1, 0, -2, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, -6, 0, 6, 0, 6, 0, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, 0, -2, 0, 4, 0, 2, 0, -3, 0, 3, 0, -3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, -1, 0, -2, 0, -1, 0},
		{18, -18, -18, 18, -18, 18, 18, -18, 12, 6, -12, -6, -12, -6, 12, 6, 12, -12, 6, -6, -12, 12, -6, 6, 9, -9, -9, 9, 9, -9, -9, 9, 8, 4, 4, 2, -8, -4, -4, -2, 6, 3, -6, -3, 6, 3, -6, -3, 6, -6, 3, -3, 6, -6, 3, -3, 4, 2, 2, 1, 4, 2, 2, 1},
		{-12, 12, 12, -12, 12, -12, -12, 12, -6, -6, 6, 6, 6, 6, -6, -6, -8, 8, -4, 4, 8, -8, 4, -4, -6, 6, 6, -6, -6, 6, 6, -6, -4, -4, -2, -2, 4, 4, 2, 2, -3, -3, 3, 3, -3, -3, 3, 3, -4, 4, -2, 2, -4, 4, -2, 2, -2, -2, -1, -1, -2, -2, -1, -1},
		{4, 0, -4, 0, -4, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, -2, 0, -2, 0, 2, 0, -2, 0, 2, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 4, 0, -4, 0, -4, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, -2, 0, -2, 0, 2, 0, -2, 0, 2, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0},
		{-12, 12, 12, -12, 12, -12, -12, 12, -8, -4, 8, 4, 8, 4, -8, -4, -6, 6, -6, 6, 6, -6, 6, -6, -6, 6, 6, -6, -6, 6, 6, -6, -4, -2, -4, -2, 4, 2, 4, 2, -4, -2, 4, 2, -4, -2, 4, 2, -3, 3, -3, 3, -3, 3, -3, 3, -2, -1, -2, -1, -2, -1, -2, -1},
		{8, -8, -8, 8, -8, 8, 8, -8, 4, 4, -4, -4, -4, -4, 4, 4, 4, -4, 4, -4, -4, 4, -4, 4, 4, -4, -4, 4, 4, -4, -4, 4, 2, 2, 2, 2, -2, -2, -2, -2, 2, 2, -2, -2, 2, 2, -2, -2, 2, -2, 2, -2, 2, -2, 2, -2, 1, 1, 1, 1, 1, 1, 1, 1}};

	    for (int i = 0; i < 64; i++)
		for (int j = 0; j < 64; j++)
		    makePolynomeCoefficients(i, j) = temp[i][j];
	}

	void setOrigin(Vector3d origin) {
		this->origin = origin;
		this->gridOrigin = origin + Vector3d(spacing/2);
	}

	/** Resize grid, also enlarges the volume as the spacing stays constant */
	void setGridSize(size_t Nx, size_t Ny, size_t Nz) {
		this->Nx = Nx;
		this->Ny = Ny;
		this->Nz = Nz;
		grid.resize(Nx * Ny * Nz);
		setOrigin(origin);
	}

	void setSpacing(double spacing) {
		this->spacing = spacing;
		setOrigin(origin);
	}

	void setReflective(bool b) {
		reflective = b;
	}

	Vector3d getOrigin() const {
		return origin;
	}
	size_t getNx() const {
		return Nx;
	}

	size_t getNy() const {
		return Ny;
	}

	size_t getNz() const {
		return Nz;
	}

	double getSpacing() const {
		return spacing;
	}

	bool isReflective() const {
		return reflective;
	}

	/** Inspector & Mutator */
	T &get(size_t ix, size_t iy, size_t iz) {
		return grid[ix * Ny * Nz + iy * Nz + iz];
	}

	/** Inspector */
	const T &get(size_t ix, size_t iy, size_t iz) const {
		return grid[ix * Ny * Nz + iy * Nz + iz];
	}

	const T &periodicGet(size_t ix, size_t iy, size_t iz) const {
	        ix = periodicBoundary(ix, Nx);
	        iy = periodicBoundary(iy, Ny);
	        iz = periodicBoundary(iz, Nz);
		return grid[ix * Ny * Nz + iy * Nz + iz];
	}

	// This is a somewhat ugly intermediary solution just to get the interpolation working. pytricubic can only operate on single floats/doubles (this is related to the fact that it uses Eigen's matrix multiplication), so I'm defining getComponent as a workaround. It's ugly bc it only works for float and Vector3f (... oh the mess when all of these template constructs break down...), but this should not be a problem currently since no other instantiations of Grid are used.
	const float getComponent(size_t ix, size_t iy, size_t iz) const {
	  throw std::runtime_error("Grid.h: Was instanciated with a type other than float or Vector3f. We can't handle this case currently, I'm sorry.");
	}

	T getValue(size_t ix, size_t iy, size_t iz) {
		return grid[ix * Ny * Nz + iy * Nz + iz];
	}

	void setValue(size_t ix, size_t iy, size_t iz, T value) {
		grid[ix * Ny * Nz + iy * Nz + iz] = value;
	}

	/** Return a reference to the grid values */
	std::vector<T> &getGrid() {
		return grid;
	}

	/** Position of the grid point of a given index */
	Vector3d positionFromIndex(int index) const {
		int ix = index / (Ny * Nz);
		int iy = (index / Nz) % Ny;
		int iz = index % Nz;
		return Vector3d(ix, iy, iz) * spacing + gridOrigin;
	}

	/** Value of a grid point that is closest to a given position */
	T closestValue(const Vector3d &position) const {
		Vector3d r = (position - gridOrigin) / spacing;
		int ix = round(r.x);
		int iy = round(r.y);
		int iz = round(r.z);
		if (reflective) {
			while ((ix < 0) or (ix > Nx))
				ix = 2 * Nx * (ix > Nx) - ix;
			while ((iy < 0) or (iy > Ny))
				iy = 2 * Ny * (iy > Ny) - iy;
			while ((iz < 0) or (iz > Nz))
				iz = 2 * Nz * (iz > Nz) - iz;
		} else {
			ix = ((ix % Nx) + Nx) % Nx;
			iy = ((iy % Ny) + Ny) % Ny;
			iz = ((iz % Nz) + Nz) % Nz;
		}
		return get(ix, iy, iz);
	}

	/** Interpolate the grid at a given position */
	T interpolate(const Vector3d &position);

	float interpolateComponent(const Vector3d &position) const {
	    // position on a unit grid
	    Vector3d r = (position - gridOrigin) / spacing;

	    // indices of lower and upper neighbors
	    int ix, iX, iy, iY, iz, iZ;
	    if (reflective) {
		    reflectiveClamp(r.x, Nx, ix, iX);
		    reflectiveClamp(r.y, Ny, iy, iY);
		    reflectiveClamp(r.z, Nz, iz, iZ);
	    } else {
		    periodicClamp(r.x, Nx, ix, iX);
		    periodicClamp(r.y, Ny, iy, iY);
		    periodicClamp(r.z, Nz, iz, iZ);
	    }

	    // linear fraction to lower and upper neighbors
	    double fx = r.x - floor(r.x);
	    double fX = 1 - fx;
	    double fy = r.y - floor(r.y);
	    double fY = 1 - fy;
	    double fz = r.z - floor(r.z);
	    double fZ = 1 - fz;

    // The following code (beginning at the end of this comment block, and
    // ending just above the "END PYTRICUBIC" mark) was taken and modified
    // from Daniel Guterding's pytricubic library, under the following
    // license:
    //
    // MIT License
    //
    // Copyright (c) 2018 Daniel Guterding
    // 
    // Permission is hereby granted, free of charge, to any person obtaining a copy
    // of this software and associated documentation files (the "Software"), to deal
    // in the Software without restriction, including without limitation the rights
    // to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    // copies of the Software, and to permit persons to whom the Software is
    // furnished to do so, subject to the following conditions:
    // 
    // The above copyright notice and this permission notice shall be included in all
    // copies or substantial portions of the Software.
    // 
    // THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    // IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    // FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    // AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    // LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    // OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    // SOFTWARE.

	    // compute interpolation polynomial for the current grid cell
	    // TODO: caching
	    Eigen::Matrix<float, 64, 1> x;
	    x <<
		// values of f(x,y,z) at each corner.
		getComponent(ix, iy, iz),
		getComponent(ix + 1, iy, iz), getComponent(ix, iy + 1, iz),
		getComponent(ix + 1, iy + 1, iz), getComponent(ix, iy, iz + 1), getComponent(ix + 1, iy, iz + 1),
		getComponent(ix, iy + 1, iz + 1), getComponent(ix + 1, iy + 1, iz + 1),
		// values of df/dx at each corner.
		0.5f * (getComponent(ix + 1, iy, iz) - getComponent(ix - 1, iy, iz)),
		0.5f * (getComponent(ix + 2, iy, iz) - getComponent(ix, iy, iz)),
		0.5f * (getComponent(ix + 1, iy + 1, iz) - getComponent(ix - 1, iy + 1, iz)),
		0.5f * (getComponent(ix + 2, iy + 1, iz) - getComponent(ix, iy + 1, iz)),
		0.5f * (getComponent(ix + 1, iy, iz + 1) - getComponent(ix - 1, iy, iz + 1)),
		0.5f * (getComponent(ix + 2, iy, iz + 1) - getComponent(ix, iy, iz + 1)),
		0.5f * (getComponent(ix + 1, iy + 1, iz + 1) - getComponent(ix - 1, iy + 1, iz + 1)),
		0.5f * (getComponent(ix + 2, iy + 1, iz + 1) - getComponent(ix, iy + 1, iz + 1)),
		// values of df/dy at each corner.
		0.5f * (getComponent(ix, iy + 1, iz) - getComponent(ix, iy - 1, iz)),
		0.5f * (getComponent(ix + 1, iy + 1, iz) - getComponent(ix + 1, iy - 1, iz)),
		0.5f * (getComponent(ix, iy + 2, iz) - getComponent(ix, iy, iz)),
		0.5f * (getComponent(ix + 1, iy + 2, iz) - getComponent(ix + 1, iy, iz)),
		0.5f * (getComponent(ix, iy + 1, iz + 1) - getComponent(ix, iy - 1, iz + 1)),
		0.5f * (getComponent(ix + 1, iy + 1, iz + 1) - getComponent(ix + 1, iy - 1, iz + 1)),
		0.5f * (getComponent(ix, iy + 2, iz + 1) - getComponent(ix, iy, iz + 1)),
		0.5f * (getComponent(ix + 1, iy + 2, iz + 1) - getComponent(ix + 1, iy, iz + 1)),
		// values of df/dz at each corner.
		0.5f * (getComponent(ix, iy, iz + 1) - getComponent(ix, iy, iz - 1)),
		0.5f * (getComponent(ix + 1, iy, iz + 1) - getComponent(ix + 1, iy, iz - 1)),
		0.5f * (getComponent(ix, iy + 1, iz + 1) - getComponent(ix, iy + 1, iz - 1)),
		0.5f * (getComponent(ix + 1, iy + 1, iz + 1) - getComponent(ix + 1, iy + 1, iz - 1)),
		0.5f * (getComponent(ix, iy, iz + 2) - getComponent(ix, iy, iz)),
		0.5f * (getComponent(ix + 1, iy, iz + 2) - getComponent(ix + 1, iy, iz)),
		0.5f * (getComponent(ix, iy + 1, iz + 2) - getComponent(ix, iy + 1, iz)),
		0.5f * (getComponent(ix + 1, iy + 1, iz + 2) - getComponent(ix + 1, iy + 1, iz)),
		// values of d2f/dxdy at each corner.
		0.25f * (getComponent(ix + 1, iy + 1, iz) - getComponent(ix - 1, iy + 1, iz) - getComponent(ix + 1, iy - 1, iz) + getComponent(ix - 1, iy - 1, iz)),
		0.25f * (getComponent(ix + 2, iy + 1, iz) - getComponent(ix, iy + 1, iz) - getComponent(ix + 2, iy - 1, iz) + getComponent(ix, iy - 1, iz)),
		0.25f * (getComponent(ix + 1, iy + 2, iz) - getComponent(ix - 1, iy + 2, iz) - getComponent(ix + 1, iy, iz) + getComponent(ix - 1, iy, iz)),
		0.25f * (getComponent(ix + 2, iy + 2, iz) - getComponent(ix, iy + 2, iz) - getComponent(ix + 2, iy, iz) + getComponent(ix, iy, iz)),
		0.25f * (getComponent(ix + 1, iy + 1, iz + 1) - getComponent(ix - 1, iy + 1, iz + 1) - getComponent(ix + 1, iy - 1, iz + 1) + getComponent(ix - 1, iy - 1, iz + 1)),
		0.25f * (getComponent(ix + 2, iy + 1, iz + 1) - getComponent(ix, iy + 1, iz + 1) - getComponent(ix + 2, iy - 1, iz + 1) + getComponent(ix, iy - 1, iz + 1)),
		0.25f * (getComponent(ix + 1, iy + 2, iz + 1) - getComponent(ix - 1, iy + 2, iz + 1) - getComponent(ix + 1, iy, iz + 1) + getComponent(ix - 1, iy, iz + 1)),
		0.25f * (getComponent(ix + 2, iy + 2, iz + 1) - getComponent(ix, iy + 2, iz + 1) - getComponent(ix + 2, iy, iz + 1) + getComponent(ix, iy, iz + 1)),
		// values of d2f/dxdz at each corner.
		0.25f * (getComponent(ix + 1, iy, iz + 1) - getComponent(ix - 1, iy, iz + 1) - getComponent(ix + 1, iy, iz - 1) + getComponent(ix - 1, iy, iz - 1)),
		0.25f * (getComponent(ix + 2, iy, iz + 1) - getComponent(ix, iy, iz + 1) - getComponent(ix + 2, iy, iz - 1) + getComponent(ix, iy, iz - 1)),
		0.25f * (getComponent(ix + 1, iy + 1, iz + 1) - getComponent(ix - 1, iy + 1, iz + 1) - getComponent(ix + 1, iy + 1, iz - 1) + getComponent(ix - 1, iy + 1, iz - 1)),
		0.25f * (getComponent(ix + 2, iy + 1, iz + 1) - getComponent(ix, iy + 1, iz + 1) - getComponent(ix + 2, iy + 1, iz - 1) + getComponent(ix, iy + 1, iz - 1)),
		0.25f * (getComponent(ix + 1, iy, iz + 2) - getComponent(ix - 1, iy, iz + 2) - getComponent(ix + 1, iy, iz) + getComponent(ix - 1, iy, iz)),
		0.25f * (getComponent(ix + 2, iy, iz + 2) - getComponent(ix, iy, iz + 2) - getComponent(ix + 2, iy, iz) + getComponent(ix, iy, iz)),
		0.25f * (getComponent(ix + 1, iy + 1, iz + 2) - getComponent(ix - 1, iy + 1, iz + 2) - getComponent(ix + 1, iy + 1, iz) + getComponent(ix - 1, iy + 1, iz)),
		0.25f * (getComponent(ix + 2, iy + 1, iz + 2) - getComponent(ix, iy + 1, iz + 2) - getComponent(ix + 2, iy + 1, iz) + getComponent(ix, iy + 1, iz)),
		// values of d2f/dydz at each corner.
		0.25f * (getComponent(ix, iy + 1, iz + 1) - getComponent(ix, iy - 1, iz + 1) - getComponent(ix, iy + 1, iz - 1) + getComponent(ix, iy - 1, iz - 1)),
		0.25f * (getComponent(ix + 1, iy + 1, iz + 1) - getComponent(ix + 1, iy - 1, iz + 1) - getComponent(ix + 1, iy + 1, iz - 1) + getComponent(ix + 1, iy - 1, iz - 1)),
		0.25f * (getComponent(ix, iy + 2, iz + 1) - getComponent(ix, iy, iz + 1) - getComponent(ix, iy + 2, iz - 1) + getComponent(ix, iy, iz - 1)),
		0.25f * (getComponent(ix + 1, iy + 2, iz + 1) - getComponent(ix + 1, iy, iz + 1) - getComponent(ix + 1, iy + 2, iz - 1) + getComponent(ix + 1, iy, iz - 1)),
		0.25f * (getComponent(ix, iy + 1, iz + 2) - getComponent(ix, iy - 1, iz + 2) - getComponent(ix, iy + 1, iz) + getComponent(ix, iy - 1, iz)),
		0.25f * (getComponent(ix + 1, iy + 1, iz + 2) - getComponent(ix + 1, iy - 1, iz + 2) - getComponent(ix + 1, iy + 1, iz) + getComponent(ix + 1, iy - 1, iz)),
		0.25f * (getComponent(ix, iy + 2, iz + 2) - getComponent(ix, iy, iz + 2) - getComponent(ix, iy + 2, iz) + getComponent(ix, iy, iz)),
		0.25f * (getComponent(ix + 1, iy + 2, iz + 2) - getComponent(ix + 1, iy, iz + 2) - getComponent(ix + 1, iy + 2, iz) + getComponent(ix + 1, iy, iz)),
		// values of d3f/dxdydz at each corner.
		0.125f * (getComponent(ix + 1, iy + 1, iz + 1) - getComponent(ix - 1, iy + 1, iz + 1) - getComponent(ix + 1, iy - 1, iz + 1) + getComponent(ix - 1, iy - 1, iz + 1) - getComponent(ix + 1, iy + 1, iz - 1) + getComponent(ix - 1, iy + 1, iz - 1) + getComponent(ix + 1, iy - 1, iz - 1) - getComponent(ix - 1, iy - 1, iz - 1)),
		0.125f * (getComponent(ix + 2, iy + 1, iz + 1) - getComponent(ix, iy + 1, iz + 1) - getComponent(ix + 2, iy - 1, iz + 1) + getComponent(ix, iy - 1, iz + 1) - getComponent(ix + 2, iy + 1, iz - 1) + getComponent(ix, iy + 1, iz - 1) + getComponent(ix + 2, iy - 1, iz - 1) - getComponent(ix, iy - 1, iz - 1)),
		0.125f * (getComponent(ix + 1, iy + 2, iz + 1) - getComponent(ix - 1, iy + 2, iz + 1) - getComponent(ix + 1, iy, iz + 1) + getComponent(ix - 1, iy, iz + 1) - getComponent(ix + 1, iy + 2, iz - 1) + getComponent(ix - 1, iy + 2, iz - 1) + getComponent(ix + 1, iy, iz - 1) - getComponent(ix - 1, iy, iz - 1)),
		0.125f * (getComponent(ix + 2, iy + 2, iz + 1) - getComponent(ix, iy + 2, iz + 1) - getComponent(ix + 2, iy, iz + 1) + getComponent(ix, iy, iz + 1) - getComponent(ix + 2, iy + 2, iz - 1) + getComponent(ix, iy + 2, iz - 1) + getComponent(ix + 2, iy, iz - 1) - getComponent(ix, iy, iz - 1)),
		0.125f * (getComponent(ix + 1, iy + 1, iz + 2) - getComponent(ix - 1, iy + 1, iz + 2) - getComponent(ix + 1, iy - 1, iz + 2) + getComponent(ix - 1, iy - 1, iz + 2) - getComponent(ix + 1, iy + 1, iz) + getComponent(ix - 1, iy + 1, iz) + getComponent(ix + 1, iy - 1, iz) - getComponent(ix - 1, iy - 1, iz)),
		0.125f * (getComponent(ix + 2, iy + 1, iz + 2) - getComponent(ix, iy + 1, iz + 2) - getComponent(ix + 2, iy - 1, iz + 2) + getComponent(ix, iy - 1, iz + 2) - getComponent(ix + 2, iy + 1, iz) + getComponent(ix, iy + 1, iz) + getComponent(ix + 2, iy - 1, iz) - getComponent(ix, iy - 1, iz)),
		0.125f * (getComponent(ix + 1, iy + 2, iz + 2) - getComponent(ix - 1, iy + 2, iz + 2) - getComponent(ix + 1, iy, iz + 2) + getComponent(ix - 1, iy, iz + 2) - getComponent(ix + 1, iy + 2, iz) + getComponent(ix - 1, iy + 2, iz) + getComponent(ix + 1, iy, iz) - getComponent(ix - 1, iy, iz)),
		0.125f * (getComponent(ix + 2, iy + 2, iz + 2) - getComponent(ix, iy + 2, iz + 2) - getComponent(ix + 2, iy, iz + 2) + getComponent(ix, iy, iz + 2) - getComponent(ix + 2, iy + 2, iz) + getComponent(ix, iy + 2, iz) + getComponent(ix + 2, iy, iz) - getComponent(ix, iy, iz));

	    // Convert voxel values and partial derivatives to interpolation coefficients.
	    // pytricubic caches these coefficients in case the next call queries the same voxel.
	    // We can't do this currently since interpolate is const.
	    Eigen::Matrix<float, 64, 1> _coefs;
	    _coefs = makePolynomeCoefficients * x;


	    // evaluate the interpolation polynomial at (fx, fy, fz)
	    float result(0.);
	    int ijkn(0);
	    float fzpow(1);

	    for (int k = 0; k < 4; ++k) {
		    float fypow(1);
		    for (int j = 0; j < 4; ++j) {
			    result += fypow * fzpow * (_coefs[ijkn] + fx * (_coefs[ijkn + 1] + fx * (_coefs[ijkn + 2] + fx * _coefs[ijkn + 3])));
			    ijkn += 4;
			    fypow *= fy;
		    }
		    fzpow *= fz;
	    }
	    return result;

	    // END PYTRICUBIC
	}
};

template <>
inline const float Grid<float>::getComponent(size_t ix, size_t iy, size_t iz) const {
  float result = periodicGet(ix, iy, iz);
  return result;
}

template <>
inline const float Grid<Vector3f>::getComponent(size_t ix, size_t iy, size_t iz) const {
  Vector3f result = periodicGet(ix, iy, iz);
  if (component == 0) {
    return result.x;
  } else if (component == 1) {
    return result.y;
  } else if (component == 2) {
    return result.z;
  }
}

template <>
inline float Grid<float>::interpolate(const Vector3d &position){
  component = 0;
  float result = interpolateComponent(position);
  return result;
}

template <>
inline Vector3f Grid<Vector3f>::interpolate(const Vector3d &position) {
  Vector3f result = Vector3f(0.0f);
  component = 0;
  result.x = interpolateComponent(position);
  component = 1;
  result.y = interpolateComponent(position);
  component = 2;
  result.z = interpolateComponent(position);
  return result;
}

typedef Grid<Vector3f> VectorGrid;
typedef Grid<float> ScalarGrid;
/** @}*/

} // namespace crpropa

#endif // CRPROPA_GRID_H
