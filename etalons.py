import numpy as np
import json

WAVELENGTH = 532.0 * 1e-3 # um
BEAM_ANGLE = 4.1

c = 299792.   # um/ns

# www.refractiveindex.info

indexies_of_refraction = {
        "fused_silica":1.4585,
        "bk7":1.5195
        }

class Etalon:

    _delta = 0 # 0.0318

    def __init__(self, thickness:"mm", n):
        self.thickness = thickness * 1e3  # convert to um
        self._n = n
    
    @property
    def vpf(self) -> "km/s":
        """Return velocity per fringe shift due to etalon 
        units of um/ns (or km/s) 
        """
        return WAVELENGTH / ( 2. * self.tau * ( 1. + self._delta))

    @property
    def tau(self) -> "ns":
        """Time delay of etalon
        time delay of a double pass through the etalon
        """
        return 2. * self.thickness / c * (self._n - 1/self._n)

    @property
    def index_of_refraction(self):
        return self._n

    def find_offset(self) -> "mm":
        """Amount to move mirror in interferometer to account for shift in focal 
        plane due to etalon

        Returns
        -------
        (float) : units of mm
        """
        ray = Ray(0, BEAM_ANGLE)
        interface = Interface(1, self.index_of_refraction)
        propogate = Propogate(self.thickness)

        ray_1 = propogate()*interface()*ray()
        ray_2 = propogate()*ray()
        gap = ray_2[0,0] - ray_1[0,0]

        return gap/np.tan(BEAM_ANGLE*np.pi/180.) * 1e-3

def json_to_etalon(obj):
    return Etalon(
            thickness=obj["thickness"], 
            n=indexies_of_refraction[obj["material"]]
            )

class Interface:
    def __init__(self, n1, n2):
        """ Ray matrix for ray intersecting and interface with change of index
        args
        ----
        n1 (float) : index of refraction before
        n2 (float) : new index of refraction
        """
        self._n1 = n1
        self._n2 = n2
        
    def __call__(self):
        """
        """
        return self._matrix()

    def _matrix(self):
        return np.matrix([[1., 0], [0., self._n1/self._n2]])

    def operate_on_ray(self, ray):
        pass

class Propogate:
    def __init__(self, d):
        """Ray matrix for the propagation of a ray some distance, d
        ARGS
        ----
        d (float) : distance to translate ray
        """
        self._d = d

    def __call__(self):
        """
        """
        return self._matrix()

    def _matrix(self):
        return np.matrix([[1., self._d], [0., 1.]])

    def operate_on_ray(self, ray):
        pass

class Ray:
    def __init__(self, y, deg):
        """Ray to propagate through etalon system
        ARGS
        ----
        y (float): initial starting position of ray
        deg (float): initial starting angle (in degrees) for ray 
        """
        self._y = y
        self._theta = deg*np.pi/180. # convert to radians

    def _matrix(self):
        return np.matrix([[self._y],[self._theta]])

    def __call__(self):
        return self._matrix()


def dump_all_to_file(etalons):
    with open("etalons.txt", "w") as f:
        f.write("{thickness:_<20} {vpf:_<20} {offset}\n\n\n".format(
            thickness="thickness",
            vpf="VPF",
            offset="offset"
            ))
        for etalon in etalons:
            f.write("{thickness:_<20} {vpf:_<20} {offset}\n\n".format(
                thickness="{:.4} mm".format(etalon.thickness*1e-3),
                vpf="{:.4} km/s".format(etalon.vpf),
                offset="{:.4} mm".format(etalon.find_offset())
                ))

if __name__ == "__main__":

    with open("etalon_conf.json") as f:
        etalons = json.load(f, object_hook=json_to_etalon)

    for etalon in etalons:
        print("{thickness:_<20} {vpf:_<20} {tau:_<20} {offset}\n\n".format(
            thickness="{:.4} mm".format(etalon.thickness*1e-3),
            vpf="{:.4} km/s".format(etalon.vpf),
            tau="{:.4} ns".format(etalon.tau),
            offset="{:.4} mm".format(etalon.find_offset())
            ))

