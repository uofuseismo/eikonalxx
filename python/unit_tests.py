#!/usr/bin/env python3
import pyEikonalXX

def test_geometry2d():
    """
    Test 2D geometry class.
    """
    geo2d = pyEikonalXX.Geometry2D()
    dx = 4.
    dz = 5.
    nx = 6
    nz = 7
    x0 = 1.
    z0 =-1.
    geo2d.dx = dx
    geo2d.dz = dz
    geo2d.nx = nx
    geo2d.nz = nz
    geo2d.x0 = x0
    geo2d.z0 = z0
    assert geo2d.nx == nx, 'nx failed'
    assert geo2d.nz == nz, 'nz failed'
    assert abs(geo2d.dx - dx) < 1.e-14, 'dx failed'
    assert abs(geo2d.dz - dz) < 1.e-14, 'dz failed'
    assert abs(geo2d.x0 - x0) < 1.e-14, 'x0 failed'
    assert abs(geo2d.z0 - z0) < 1.e-14, 'z0 failed'
    print("Passed geometry2d")

def test_geometry3d():
    """ 
    Test 3D geometry class.
    """
    geo3d = pyEikonalXX.Geometry3D()
    dx = 4.
    dy = 3.
    dz = 5.
    nx = 6 
    ny = 8
    nz = 7 
    x0 = 1.
    y0 = 0.5
    z0 =-1.
    geo3d.dx = dx
    geo3d.dy = dy
    geo3d.dz = dz
    geo3d.nx = nx
    geo3d.ny = ny
    geo3d.nz = nz
    geo3d.x0 = x0
    geo3d.y0 = y0
    geo3d.z0 = z0
    assert geo3d.nx == nx, 'nx failed'
    assert geo3d.ny == ny, 'ny failed'
    assert geo3d.nz == nz, 'nz failed'
    assert abs(geo3d.dx - dx) < 1.e-14, 'dx failed'
    assert abs(geo3d.dy - dy) < 1.e-14, 'dy failed'
    assert abs(geo3d.dz - dz) < 1.e-14, 'dz failed'
    assert abs(geo3d.x0 - x0) < 1.e-14, 'x0 failed'
    assert abs(geo3d.y0 - y0) < 1.e-14, 'y0 failed'
    assert abs(geo3d.z0 - z0) < 1.e-14, 'z0 failed'
    print("Passed geometry3d")

def test_solverOptions():
    """
    Test the solver options.
    """
    options = pyEikonalXX.SolverOptions()
    tol = 1.e-2
    n_gauss_sweeps = 3
    eps = 4
    algorithm = pyEikonalXX.SolverAlgorithm.level_set_method
    verbosity = pyEikonalXX.Verbosity.debug
 
    options.tolerance = tol
    options.number_of_sweeps = n_gauss_sweeps
    options.factored_eikonal_equation_solver_radius = eps
    options.algorithm = algorithm
    options.verbosity = verbosity

    assert abs(options.tolerance - tol) < 1.e-14, 'tolerance failed'
    assert options.number_of_sweeps == n_gauss_sweeps, 'n sweeps failed'
    assert options.factored_eikonal_equation_solver_radius == eps, 'factored eikonal solver radius failed' 
    assert options.algorithm == algorithm, 'algorithm failed'
    assert options.verbosity == verbosity, 'verbosity failed'
    print("Passed solverOptions")

def test_velocityModel2d():
    """
    Test the 2D velocity model.
    """
    geo2d = pyEikonalXX.Geometry2D()
    geo2d.dx = 1
    geo2d.dz = 2
    geo2d.nx = 11
    geo2d.nz = 13
    geo2d.x0 = 0
    model2d = pyEikonalXX.VelocityModel2D()
    model2d.initialize(geo2d)

    assert model2d.is_initialized(), 'initialized failed'
    geo_back = model2d.get_geometry()
    assert geo2d.nx == geo_back.nx, 'nx failed'
    assert geo2d.nz == geo_back.nz, 'nz failed'
    assert abs(geo2d.dx - geo_back.dx) < 1.e-14, 'dx failed'
    assert abs(geo2d.dz - geo_back.dz) < 1.e-14, 'dz failed'
    assert abs(geo2d.x0 - geo_back.x0) < 1.e-14, 'x0 failed'
    assert abs(geo2d.z0 - geo_back.z0) < 1.e-14, 'z0 failed'

    print("Passed velocityModel2d")

def test_source2d():
    """
    Test the 2D source.
    """
    geo2d = pyEikonalXX.Geometry2D()
    geo2d.dx = 100
    geo2d.dz = 101.
    geo2d.nx = 55
    geo2d.nz = 25
    geo2d.x0 = 1
    geo2d.z0 = 2
    source2d = pyEikonalXX.Source2D()
    source2d.geometry = geo2d
    x_src = geo2d.x0 + geo2d.nx/4*geo2d.dx
    z_src = geo2d.z0 + geo2d.nz/8*geo2d.dz
    source2d.x = x_src
    source2d.z = z_src

    geo_back = source2d.geometry
    assert geo2d.nx == geo_back.nx, 'nx failed'
    assert geo2d.nz == geo_back.nz, 'nz failed'
    assert abs(geo2d.dx - geo_back.dx) < 1.e-14, 'dx failed'
    assert abs(geo2d.dz - geo_back.dz) < 1.e-14, 'dz failed'
    assert abs(geo2d.x0 - geo_back.x0) < 1.e-14, 'x0 failed'
    assert abs(geo2d.z0 - geo_back.z0) < 1.e-14, 'z0 failed'
    assert abs(source2d.x - x_src) < 1.e-14, 'x failed'
    assert abs(source2d.z - z_src) < 1.e-14, 'z failed'
    source2d.set_z_to_free_surface()
    assert abs(source2d.z - geo2d.z0) < 1.e-14, 'z to free surface failed'
    print("Passed source2d")

if __name__ == "__main__":
    test_geometry2d()
    test_geometry3d()
    test_solverOptions()
    test_velocityModel2d()
    test_source2d()
