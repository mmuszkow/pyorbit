#!/usr/bin/env python3

from math import sin, cos, radians, pi, sqrt, atan, asin, acos, degrees, atan2
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D
import shapefile
import datetime

# flat XY-ellipse rotated in 3d space
class Orbit:
    def __init__(self, a, e, rot):
        self.a = a
        self.e = e
        self.rx = rot[0]
        self.ry = rot[1]
        self.rz = rot[2]

    # z,y,x-rotation matrices
    # {{cos(Z),-sin(Z),0},{sin(Z),cos(Z),0},{0,0,1}}.{{cos(Y),0,sin(Y)},{0,1,0},{-sin(Y),0,cos(Y)}}.{{1,0,0},{0,cos(X),-sin(X)},{0,sin(X),cos(X)}}.{{x},{y},{z}}
    def xyz(self, angle_rad):
        r = (self.a * (1 - self.e**2)) / (1 + self.e * cos(angle_rad)) # Kepler's orbit equation
        ex = r * cos(angle_rad)
        ey = r * sin(angle_rad)
        ez = 0
        x = ex * cos(self.ry) * cos(self.rz) + \
            ey * (cos(self.rz) * sin(self.rx) * sin(self.ry) - cos(self.rx) * sin(self.rz)) + \
            ez * (cos(self.rx) * cos(self.rz) * sin(self.ry) + sin(self.rx) * sin(self.rz))
        y = ex * cos(self.ry) * sin(self.rz) + \
            ez * (cos(self.rx) * sin(self.ry) * sin(self.rz) - cos(self.rz) * sin(self.rx)) + \
            ey * (cos(self.rx) * cos(self.rz) + sin(self.rx) * sin(self.ry) * sin(self.rz))
        z = ez * cos(self.rx) * cos(self.ry) + \
            ey * sin(self.rx) * cos(self.ry) - \
            ex * sin(self.ry)
        return (x, y, z)

class Satellite(Orbit):
    def __init__(self, a, e, rot, period, angle0, name=None):
        Orbit.__init__(self, a, e, rot)
        self.name = name
        self.period = period
        self.angle0 = angle0

    def position(self, t):
        angle = self.angle0 + (t/self.period)*2*pi
        x, y, z = self.xyz(angle)
        earth_rot = -(t/(24*60*60)) * 2*pi
        rx = x * cos(earth_rot) - y * sin(earth_rot)
        ry = y * cos(earth_rot) + x * sin(earth_rot)
        return (rx, ry, z)

    def __str__(self):
        return 'a=%f, e=%f, rot=(%f,%f,%f), period=%f, angle0=%f' % \
                (self.a, self.e, self.rx, self.ry, self.rz, self.period, self.angle0)

# check mpl_toolkits.basemap.Basemap for more advanced stuff
class Earth:
    def __init__(self):
        # ellipsoidal model
        self.R_EQUAT = 6371.137       # WSG-84
        self.R_POLAR = 6356.752314245 # WSG-84

        # continents outlines
        sf = shapefile.Reader('ne_110m_coastline.shp')
        self.outlines = []
        for s in sf.shapes():
            outline = [self.lon_lat2xyz(p[0], p[1], 0) for p in s.points]
            self.outlines.append(outline)

        self.satellites = []
        self.markers = []

    def lon_lat2xyz(self, lon, lat, alt):
        theta = radians(lon+180)
        phi = radians(lat+90)
        x = -(self.R_EQUAT+alt) * sin(phi) * cos(theta)
        y = -(self.R_EQUAT+alt) * sin(phi) * sin(theta)
        z = -(self.R_POLAR+alt) * cos(phi)
        return (x, y, z) 

    def xyz2lon_lat(self, x, y, z):
        r = sqrt(x**2 + y**2 + z**2)
        phi = -acos(z/r) + (pi/2)
        theta = atan2(y, x)
        return (degrees(theta), degrees(phi))

    def central_angle(self, sat, t):
        x, y, z = sat.position(t)
        lon, lat = self.xyz2lon_lat(x, y, z)
        r_h = sqrt(x**2 + y**2 + z**2) # sat distance from Earth center
        x, y, z = self.lon_lat2xyz(lon, lat, 0)
        r = sqrt(x**2 + y**2 + z**2)  # sphere radius in this specific point
        dlon = degrees(acos(r/r_h))
        
        # latitude is a bit more tricky cause we use the spheroid model
        # I use this system of equations to determine the 2 tangent points
        # x**2/A**2 + y**2/B**2 = 1, Y = a*X + b, y = a*x + b
        #A = self.R_EQUAT
        #B = self.R_POLAR
        #X = sqrt(r_h**2-z**2)
        #Y = z
        #x1 = A**2*(B**2*X + Y*sqrt(-A**2*B**2 + A**2*Y**2 + B**2*X**2))/(A**2*Y**2 + B**2*X**2)
        #y1 = B**2*(A**2*Y - X*sqrt(-A**2*B**2 + A**2*Y**2 + B**2*X**2))/(A**2*Y**2 + B**2*X**2)
        #x2 = A**2*(B**2*X - Y*sqrt(-A**2*B**2 + A**2*Y**2 + B**2*X**2))/(A**2*Y**2 + B**2*X**2)
        #y2 = B**2*(A**2*Y + X*sqrt(-A**2*B**2 + A**2*Y**2 + B**2*X**2))/(A**2*Y**2 + B**2*X**2)
        #print(self.xyz2lon_lat(x1, y1, z), (lon, lat))
        # the difference is around 1 deg, it's negligible so we use the simpler method

        return dlon
        
    # this is used to avoid avoid artifacts when 
    # crossing -180/180 parallel for 2d plot
    def _safe_cross(self, lon_lat):
        paths = [[]]
        for i in range(len(lon_lat)-1):
            paths[-1].append(lon_lat[i])
            if abs(lon_lat[i][0] - lon_lat[i+1][0]) >= 180:
                if lon_lat[i+1][0] > 0:
                    paths[-1].append((-180, lon_lat[i][1]))
                    paths.append([(180, lon_lat[i+1][1])])
                else:
                    paths[-1].append((180, lon_lat[i][1]))
                    paths.append([(-180, lon_lat[i+1][1])])
        paths[-1].append(lon_lat[-1])
        return paths

    def draw2d(self):
        fig, ax = plt.subplots(figsize=(9,7))
        ax.set_xlim(-180, 180)
        ax.set_ylim(-90, 90)
        ax.grid(True)

        # draw countries borders
        for outline in self.outlines:
            lon_lat = [self.xyz2lon_lat(p[0], p[1], p[2]) for p in outline]
            for path in self._safe_cross(lon_lat):
                ax.plot([p[0] for p in path], [p[1] for p in path], color='gray', linewidth=1)
        
        # draw satellites orbits predictions and access range
        patches = []
        for satellite in self.satellites:
            now = int(datetime.datetime.utcnow().timestamp())
            
            # draw orbit prediction
            trange = 10*int(satellite.period/2) # draw this period both ways
            xyz = [satellite.position(t) for t in range(now-trange, now+trange, 60)]
            lon_lat = [self.xyz2lon_lat(p[0], p[1], p[2]) for p in xyz]
            color = None # we want all path parts and current position in the same color
            for path in self._safe_cross(lon_lat):
                pplot = ax.plot([p[0] for p in path], [p[1] for p in path], color=color)
                if not color: color = pplot[0].get_color()
            
            x, y, z = satellite.position(now)
            
            # draw access range
            alpha = self.central_angle(satellite, now)
            lon, lat = self.xyz2lon_lat(x, y, z)
            path = [(lon + alpha * sin(radians(t)), lat + alpha * cos(radians(t))) for t in range(0, 360)]
            ax.plot([p[0] for p in path], [p[1] for p in path], '--', color=color, linewidth=1)

            # draw current position
            arrow0 = self.xyz2lon_lat(x, y, z)
            x, y, z = satellite.position(now + 1)
            arrow1 = self.xyz2lon_lat(x, y, z)
            ax.annotate('', xytext=arrow0, xy=arrow1, \
                        arrowprops=dict(ec=color, fc='white', arrowstyle='simple'))
            if satellite.name: patches.append(mpatches.Patch(color=color, label=satellite.name))
         
        ax.scatter([p[0] for p in self.markers], [p[1] for p in self.markers])
        plt.legend(loc='center', bbox_to_anchor=(0.5, 1.05), ncol=len(patches), handles=patches)
        plt.show()

    def draw3d(self):
        fig = plt.figure(figsize=(9,7))
        ax = fig.add_subplot(111, projection='3d')
        ax.set_aspect(1)
        
        for outline in self.outlines:
            xs = [p[0] for p in outline]
            ys = [p[1] for p in outline]
            zs = [p[2] for p in outline]
            ax.plot(xs, ys, zs, color='gray', linewidth=1)

        for satellite in self.satellites:
            now = int(datetime.datetime.utcnow().timestamp())
            trange = int(satellite.period/2) # draw this period both ways
            xyz = [satellite.position(t) for t in range(now-trange, now+trange, 60)]
            xs = [p[0] for p in xyz]
            ys = [p[1] for p in xyz]
            zs = [p[2] for p in xyz]
            ax.plot(xs, ys, zs)
            x, y, z = satellite.position(now)
            ax.scatter([x], [y], [z])

        xyz = [self.lon_lat2xyz(p[0], p[1], 0) for p in self.markers]
        xs = [p[0] for p in xyz]
        ys = [p[1] for p in xyz]
        zs = [p[2] for p in xyz]
        ax.scatter(xs, ys, zs)

        fig.tight_layout() # to fill empty space
        plt.show()

if __name__ == '__main__':
    earth = Earth()
    with open('orbits.txt') as f:
        for line in f.readlines():
            val = line.split()
            earth.satellites.append(Satellite( \
                    float(val[1]), float(val[2]), (float(val[3]), float(val[4]), float(val[5])), \
                    float(val[6]), float(val[7]), int(val[0])))
    earth.draw2d()
    earth.draw3d()

