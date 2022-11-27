import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def list_split(listA, n):
    """Solitting a list with division"""
    for x in range(0, len(listA), n):
        every_chunk = listA[x: n+x]

        if len(every_chunk) < n:
            every_chunk = every_chunk + \
                [None for y in range(n-len(every_chunk))]
        yield every_chunk



def get_brillouin_zone_3d(cell):
    """
    Generate the Brillouin Zone of a given cell. The BZ is the Wigner-Seitz cell
    of the reciprocal lattice, which can be constructed by Voronoi decomposition
    to the reciprocal lattice.  A Voronoi diagram is a subdivision of the space
    into the nearest neighborhoods of a given set of points. 
    """

    cell = np.asarray(cell, dtype=float)
    assert cell.shape == (3, 3)

    px, py, pz = np.tensordot(cell, np.mgrid[-1:2, -1:2, -1:2], axes=[0, 0])
    points = np.c_[px.ravel(), py.ravel(), pz.ravel()]

    from scipy.spatial import Voronoi
    vor = Voronoi(points)

    bz_facets = []
    bz_ridges = []
    bz_vertices = []

    for pid, rid in zip(vor.ridge_points, vor.ridge_vertices):
        # WHY 13 ????
        # The Voronoi ridges/facets are perpendicular to the lines drawn between the
        # input points. The 14th input point is [0, 0, 0].
        if(pid[0] == 13 or pid[1] == 13):
            bz_ridges.append(vor.vertices[np.r_[rid, [rid[0]]]])
            bz_facets.append(vor.vertices[rid])
            bz_vertices += rid

    bz_vertices = list(set(bz_vertices))

    return vor.vertices[bz_vertices], bz_ridges, bz_facets

# if __name__ == "__main__":
#     cell = np.array([[6.28318521724735,	-0.0,	0.0],
#                      [-0.0,	6.28318521724735, 0.0],
#                      [0.0, 0.0, 6.28318521724735,]])
#     icell = np.linalg.inv(cell).T                
#     b1, b2, b3 = np.linalg.norm(icell, axis=1)   

#     v, e, f = get_brillouin_zone_3d(icell)

#     fig = plt.figure(
#         figsize=(2, 2), dpi=300
#     )
#     ax = plt.subplot(111, projection='3d')

#     for xx in e:
#         ax.plot(xx[:, 0], xx[:, 1], xx[:, 2], color='k', lw=1.0)

#     ax.set_xlim(-b1, b1)
#     ax.set_ylim(-b2, b2)
#     ax.set_zlim(-b3, b3)

#     plt.show()
    