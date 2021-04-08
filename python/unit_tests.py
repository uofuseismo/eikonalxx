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

if __name__ == "__main__":
    test_geometry2d()
    test_geometry3d()
