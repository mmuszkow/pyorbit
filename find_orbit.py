#!/usr/bin/env python3

from pyorbit import Earth, Satellite
from math import sqrt, pi, sin, cos, asin, acos, degrees, radians, atan2
import datetime
import sys
import numpy as np
import scipy.optimize

def perp_error(params, txyz):
    a, b, c, d = params
    length = np.sqrt(a**2 + b**2 + c**2)
    return np.mean([np.abs(a*x + b*y + c*z + d) / length for t, x, y, z in txyz])

def ellipse_error(params, r, theta):
    a, e = params
    return r - (a * (1 - e**2)/(1 - e*np.cos(theta)))

def jac(p, r, theta):
    a, e = p
    da = (1 - e**2)/(1 - e*np.cos(theta))
    de = (-2*a*e*(1-e*np.cos(theta)) + a*(1-e**2)*np.cos(theta))/(1 - e*np.cos(theta))**2
    return -da,  -de
    #return np.array((-da, -de)).T

class N2yo:
    def __init__(self, filename):
        self.earth = Earth()
        self.satellites = {}
        with open(filename) as f:
            for line in f.readlines():
                vals = line.strip().split()
                norad_id = int(vals[0])
                dt = datetime.datetime.strptime(vals[1], '%Y-%m-%d_%H:%M:%S')
                t = dt.timestamp()
                lon = float(vals[2])
                lat = float(vals[3])
                alt = float(vals[4])
                if not norad_id in self.satellites: self.satellites[norad_id] = []
                x, y, z = self.earth.lon_lat2xyz(lon, lat, alt)
                self.satellites[norad_id].append( (t, x, y, z) )

    def find_rotation_angles(self, plane1, plane2):
        # rot = {{cos(Z),-sin(Z),0},{sin(Z),cos(Z),0},{0,0,1}}.{{cos(Y),0,sin(Y)},{0,1,0},{-sin(Y),0,cos(Y)}}.{{1,0,0},{0,cos(X),-sin(X)},{0,sin(X),cos(X)}}
        # rot * A = B
        # rot = B * A^-1
        rot_mat = plane1 * plane2.I
        if abs(rot_mat[2,0] != 1):
            ay1 = -asin(rot_mat[2,0])
            ay2 = pi - ay1
            ax1 = atan2(rot_mat[2,1]/cos(ay1), rot_mat[2,2]/cos(ay1))
            ax2 = atan2(rot_mat[2,1]/cos(ay2), rot_mat[2,2]/cos(ay2))
            az1 = atan2(rot_mat[1,0]/cos(ay1), rot_mat[0,0]/cos(ay1))
            az2 = atan2(rot_mat[1,0]/cos(ay2), rot_mat[0,0]/cos(ay2))
            return [[ax1, ay1, az1], [ax2, ay2, az2]]
        else:
            az = 0
            if rot_mat[2,0] == -1:
                ay = pi/2
                ax = az + atan2(rot_mat[0,1], rot_mat[0,2])
                return [[ax, ay, az]]
            else:
                ay = -pi/2
                ax = -az + atan2(-rot_mat[0,1], -rot_mat[0,2])
                return [[ax, ay, az]]

    def rotate(self, txyz, rot):
        return [(t, \
            x * cos(rot[1]) * cos(rot[2]) + \
            y * (cos(rot[2]) * sin(rot[0]) * sin(rot[1]) - cos(rot[0]) * sin(rot[2])) + \
            z * (cos(rot[0]) * cos(rot[2]) * sin(rot[1]) + sin(rot[0]) * sin(rot[2])), \
            x * cos(rot[1]) * sin(rot[2]) + \
            z * (cos(rot[0]) * sin(rot[1]) * sin(rot[2]) - cos(rot[2]) * sin(rot[0])) + \
            y * (cos(rot[0]) * cos(rot[2]) + sin(rot[0]) * sin(rot[1]) * sin(rot[2])), \
            z * cos(rot[0]) * cos(rot[1]) + y * sin(rot[0]) * cos(rot[1]) - x * sin(rot[1])) \
            for t, x, y, z in txyz]

    def find_orbit(self, norad_id):
        data = self.satellites[norad_id]

        # "unrotate"
        unrotated = []
        for t, x, y, z in data:
            earth_rot = (t/(24*60*60)) * 2*pi
            rx = x * cos(earth_rot) - y * sin(earth_rot)
            ry = y * cos(earth_rot) + x * sin(earth_rot)
            unrotated.append((t, rx, ry, z))

        # find plane best fitting all points (perpendicular)
        cons = ({'type': 'eq', 'fun': lambda p: p[0]**2 + p[1]**2 + p[2]**2 - 1}) # keep vector len=1
        plane = scipy.optimize.minimize(perp_error, [0, 0, 1, 0], args=unrotated, constraints=cons).x
        plane = np.matrix(plane[0:3]).T

        # find rotation angles
        xy_plane = np.matrix([[0], [0], [1]])
        rot = self.find_rotation_angles(xy_plane, plane) 

        # compute points coordinates on XY plane
        flat = self.rotate(unrotated, rot[0])

        # transformate x,y coordinates to ellipse coordinates
        r = [sqrt(x**2 + y**2) for t, x, y, z in flat]
        theta = [atan2(x, y) for t, x, y, z in flat]

        # fit ellipse
        ellipse = scipy.optimize.leastsq(ellipse_error, (7200, 0.01), Dfun=jac, args=(r, theta), col_deriv=True)    
        a = ellipse[0][0]
        e = ellipse[0][1]

        # compute period - Kepler again
        period = 2*pi*sqrt( ((a*1000)**3) / 3.9860044189e14)

        # rotation
        rot = (pi, acos(plane[2,0]), atan2(plane[1,0], plane[0,0]))

        mina = 0
        maxa = 2*pi
        lerr = 0
        while True:
            middle = (mina+maxa)/2
            langle = (mina+middle)/2
            rangle = (maxa+middle)/2
            lsat = Satellite(a, e, rot, period, langle)
            rsat = Satellite(a, e, rot, period, rangle)
            lerr = 0
            rerr = 0
            for t, x, y, z in data:
                lx, ly, lz = lsat.position(t)
                rx, ry, rz = rsat.position(t)
                lerr += sqrt((x-lx)**2 + (y-ly)**2 + (z-lz)**2)
                rerr += sqrt((x-rx)**2 + (y-ry)**2 + (z-rz)**2)
            lerr /= len(data)
            rerr /= len(data)
            if lerr < rerr:
                maxa = rangle
            elif rerr < lerr:
                mina = langle
            else:
                break
        
        print('err=', lerr, file=sys.stderr)

        #flat2 = []
        #sat0 = Satellite(a, e, (0, 0, 0), period, 0)
        #for t in range(0, int(period/2)+61, 60):
            #x, y, z = sat0.position(t)
            #flat2.append((t, x, y, z))
        #unrotated2 = self.rotate(flat2, rot)
        #self.earth.outlines.append([(x, y, z) for t, x, y, z in unrotated])
        #self.earth.outlines.append([(x, y, z) for t, x, y, z in flat])
        #self.earth.outlines.append([(x, y, z) for t, x, y, z in unrotated2])
        #self.earth.outlines.append([(x, y, z) for t, x, y, z in flat2])
        #self.earth.outlines.append([(x, y, z)for t, x, y, z in data])
        #sat = Satellite(a, e, rot, period, mina)
        #self.earth.outlines.append([sat.position(t) for t, x, y, z in data])
        #self.earth.draw3d()

        return (norad_id, a, e, rot, period, mina)

n2yo = N2yo('n2yo.txt')
for norad_id in n2yo.satellites:
    norad_id, a, e, (rx, ry, rz), period, angle0 = n2yo.find_orbit(norad_id)
    print(norad_id, a, e, rx, ry, rz, period, angle0)

