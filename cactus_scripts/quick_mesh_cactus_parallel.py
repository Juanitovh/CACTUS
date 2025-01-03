#!/usr/bin/env python

import argparse
import concurrent.futures
import pyvista
import numpy as np
import sys
from progressbar import ProgressBar



import tqdm
import openmesh as om

from libraries import RotatingStrand as RS
from libraries import WattsonFunction as WF
from libraries import ReadWriteStrands as RWS
from libraries import ascii_art 
from libraries import CactusMath as CM
import pyvista as pv

import multiprocessing
from multiprocessing import get_context



def get_circle(x, y, z, r, n):
    points = np.zeros((n, 3))
    theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
    for i, angle in enumerate(theta):
        points[i, 0] = np.cos(angle) * r + x
        points[i, 1] = np.sin(angle) * r + y
        points[i, 2] = z
    return points.T

# Define Functions
def create_circles_from_strand(strand, n_mesh):
    """
    Creates circular cross-sections along the strand.
    """
    vertices = []
    v1 = strand[0:3, 1] - strand[0:3, 0]
    vz = np.array([0, 0, 1])
    points0 = get_circle(0, 0, 0, strand[3, 0], n_mesh)
    mat1 = RS.rotation_matrix_from_vectors(vz, v1)

    circle1 = mat1.dot(points0).T + strand[0:3, 0]
    vertices.append(circle1)

    for i in range(1, len(strand[0]) - 1):
        v1 = strand[0:3, i] - strand[0:3, i - 1]
        v2 = strand[0:3, i + 1] - strand[0:3, i]

        v12 = (v1 / np.linalg.norm(v1) + v2 / np.linalg.norm(v2))

        points0 = get_circle(0, 0, 0, strand[3, i], n_mesh)
        mat1 = RS.rotation_matrix_from_vectors(vz, v12)

        circlei = mat1.dot(points0).T + strand[0:3, i]
        vertices.append(circlei)

    v1 = strand[0:3, -1] - strand[0:3, -2]

    points0 = get_circle(0, 0, 0, strand[3, -1], n_mesh)
    mat1 = RS.rotation_matrix_from_vectors(vz, v1)

    circle1 = mat1.dot(points0).T + strand[0:3, -1]
    vertices.append(circle1)

    return vertices


def volume_cyl(p1, p2):
    """
    Calculate the approximate volume of a cylindrical segment between two points.
    """
    r = (p1[-1] + p2[-1]) / 2.0
    dista = np.linalg.norm(p1[0:3] - p2[0:3])
    return np.pi * r * r * dista


def process_strand(elements):
    i, strand, g_ratio, n_circles , decimate = elements
    """
    Process a single strand and return the resulting mesh components and approximate volume.
    """
    mesh_sphere = None
    volumen_approx = 0

    mesh = om.TriMesh()

    if decimate >= 0 and decimate <= 1:
        n_percentage = len(strand) ** (1 -decimate)
        n_take = int(len(strand) // n_percentage)
        strand = strand[::n_take]

    if len(strand) < 4:
        return None, None, f"Error: Strand {i} is empty. Decimation too strong."

    strand[:, 3] *= g_ratio

    if len(strand) == 1:
        center = strand[0][0:3]
        rad = strand[0][-1]
        sphere = pyvista.Sphere(radius=rad, center=center, theta_resolution=7, phi_resolution=7)
        return sphere, volumen_approx, None

    vertices = create_circles_from_strand(strand.T, n_circles)
    vers = np.concatenate(vertices)
    list_vert = [mesh.add_vertex(v) for v in vers]

    n_control_points = len(strand)
    for i in range(n_control_points - 1):
        volumen_approx += volume_cyl(strand[i], strand[i + 1])
        for j in range(n_circles):
            v1 = list_vert[n_circles * i + j]
            v2 = list_vert[n_circles * i + (j + 1) % n_circles]
            v11 = list_vert[n_circles * (i + 1) + j]
            v22 = list_vert[n_circles * (i + 1) + (j + 1) % n_circles]
            mesh.add_face(v1, v2, v11)
            mesh.add_face(v22, v11, v2)

    v1 = list_vert[0]
    for i in range(1, n_circles - 1):
        v2 = list_vert[i]
        v3 = list_vert[i + 1]
        mesh.add_face(v1, v3, v2)

    v1 = list_vert[-n_circles]
    for i in range(-n_circles + 1, -1):
        v2 = list_vert[i]
        v3 = list_vert[i + 1]
        mesh.add_face(v3, v1, v2)



    

    faces = [[3] + [vh.idx() for vh in mesh.fv(f)] for f in mesh.faces()]
    

    #faces = np.hstack([ np.concatenate([[3] , v[::-1]]) for v in res.faces] )
    res2 = pv.PolyData( vers, faces)

    return res2 #, volumen_approx, None


def parallel_process_strands(strands, g_ratio, n_circles , decimate):
    """
    Process all strands in parallel and combine the results.
    """
    # Use multiprocessing with "spawn" context

    
    inputs = [[i, strand, g_ratio, n_circles , decimate]  for i, strand in enumerate(strands)]
    n_pool = min(len(inputs), 10)
    results = []
    with get_context("spawn").Pool(n_pool) as pool:
        for result in tqdm.tqdm(pool.imap_unordered(process_strand, inputs), total=len(inputs)):
            results.append(result)

    total_volume = 0
    """ for input in inputs: """
    """     result = process_strand(input) """
    """     results.append(result) """

    final_mesh = CM.paste_binary(results)

    return final_mesh , total_volume


# Main Script
if __name__ == "__main__":

    ascii_art.display_title()
    parser = argparse.ArgumentParser(description="Parallel Strand Processing")
    parser.add_argument("-file", help="Path to the input file containing strands")
    parser.add_argument("-output", default="output_mesh.ply", help="Output mesh file")
    parser.add_argument("-g_ratio", type=float, default=1.0, help="Growth ratio multiplier")
    parser.add_argument("-n_circles", type=int, default=12, help="Number of circles for strand meshing")
    parser.add_argument("-decimate", type=float, default=0.0, help="Decimation factor (0 to 1)")
    parser.add_argument("-strand_id", type=str, default=None, help="Strand ID to process (default: all)")
    args = parser.parse_args()
    print("Reading strands from file...")

    strands, box_size = RWS.read_generic_list(args.file)
    print(f"Read {len(strands)} strands from file.")

    if args.strand_id != None:
        try:
            ids = args.strand_id.split(",")
            ids = [int(s) for s in ids]
            strands = [strands[i] for i in ids]
        except:
            print("Error: Invalid strand ID.")
            sys.exit()

    mesh_sphere, total_volume = parallel_process_strands(strands, args.g_ratio, args.n_circles , args.decimate)

    print(f"Total Volume: {total_volume}")
    mesh_sphere.save(args.output, binary = False , recompute_normals = False)
