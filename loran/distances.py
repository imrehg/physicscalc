from __future__ import division
from numpy import *

def deg2rad(deg):
    return deg/180*pi

def m2mile(m):
    return m / 1609.344

class Location:
    """ A geographical location """

    def __init__(self, loc):
        self.lat = loc[0]
        self.lon = loc[1]
        self.radlat = deg2rad(self.lat)
        self.radlon = deg2rad(self.lon)

    def tostring(self):
        return "%f,%f" %(self.lat, self.lon)

class Geoid:
    """ Set up a geoid for distance calculations """

    def __init__(self, radius = 6378137, rflattening=298.257222101):
        self.radius = radius
        self.rflattening = rflattening

    def _reducedlat(self, lat):
        """ Reduced latitude """
        return arctan((self.rflattening-1)/self.rflattening * tan(lat))

    def _centralangle(self, loc1, loc2):
        """ Central angle calculation according to the haversine formula """
        df = loc2.radlat - loc1.radlat
        dl = loc2.radlon - loc1.radlon
        ds = 2*arcsin(sqrt(sin(df/2)**2 
                           + cos(loc1.radlat)*cos(loc2.radlat)*sin(dl/2)**2))
        return ds

    def distance(self, loc1, loc2):
        """ 
        Distance between two locations on the geoid
        Lambert's formulae
        """
        ds = self._centralangle(loc1, loc2)
        P = (loc1.radlat + loc2.radlat)/2
        Q = (loc1.radlat - loc2.radlat)/2
        X = (ds - sin(ds)) * sin(P)**2 * cos(Q)**2 / cos(ds/2)**2
        Y = (ds + sin(ds)) * cos(P)**2 * sin(Q)**2 / sin(ds/2)**2
        return self.radius*(ds - (X+Y)/(2*self.rflattening))


mylocation = Location([25.01861, 121.538481])
# loran1 = Location([37.064379, 122.323876])
geo = Geoid()
# print geo.distance(mylocation, loran1)/1600
# l1 = Location([25.04945, 121.43558])
# l2 = Location([25.05232, 121.54295])
# print geo.distance(l2, l1)

data = loadtxt("data.csv", delimiter=',', skiprows=1,
               dtype={'names': ('chain', 'GRI', 'station', 'locname', 'loccoord', 'distance', 'power'),
                      'formats': ('a30', 'i4', 'a1', 'a40', 'a30', 'a1', 'i4')})
base = "http://maps.google.com/maps/api/staticmap?"
center = "center=%s" %mylocation.tostring()
center = "center=31.58,128.56"
zoom = "zoom=4"
size = "size=1000x500"
maptype = "maptype=roadmap"
sensor = "sensor=false"
params = [size, maptype, sensor]
stationlist = {}
for station in data:
    if stationlist.has_key(station['locname']):
        continue

    loran = Location([float(l) for l in station['loccoord'].split(' ')])
    print "%s (%d/%s) : %f" %(station['locname'], station['GRI'], station['station'],
                              m2mile(geo.distance(mylocation, loran)))
    stationlist[station['locname']] = m2mile(geo.distance(mylocation, loran))
    params.append("markers=color:blue|label:%s|%s" %(station['locname'][0], loran.tostring()))
params.append("markers=color:yellow|label:T|%s" %(mylocation.tostring()))

mymap = ''.join([base, '&'.join(params)])
import webbrowser
webbrowser.open(mymap)
