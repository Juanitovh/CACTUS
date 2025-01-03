#!/usr/bin/env python

import argparse
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pyvista
import openmesh as om
import progressbar
from libraries import RotatingStrand as RS
from libraries import WattsonFunction as WF
from libraries import ReadWriteStrands as RWS

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Process strand data to generate meshes.")
    parser.add_argument("-file", type=str, required=True, help="Input file containing strand data.")
    parser.add_argument("-file_out", type=str, default="test2.ply", help="Output mesh file.")
    parser.add_argument("-g_ratio", type=float, default=1.0, help="Growth ratio for scaling radii.")
    parser.add_argument("-strand_id", type=str, help="Comma-separated list of strand IDs to process.")
    parser.add_argument("-decimate", type=float, help="Percentage to decimate" , default = 0.0)
    return parser.parse_args()

def tetra_volume(p1, p2, p3):
    v321 = p3[0] * p2[1] * p1[2]
    v231 = p2[0] * p3[1] * p1[2]
    v312 = p3[0] * p1[1] * p2[2]
    v132 = p1[0] * p3[1] * p2[2]
    v213 = p2[0] * p1[1] * p3[2]
    v123 = p1[0] * p2[1] * p3[2]
    return (1 / 6.0) * (-v321 + v231 + v312 - v132 - v213 + v123)

def get_circle(x, y, z, r, n):
    points = np.zeros((n, 3))
    theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
    for i, angle in enumerate(theta):
        points[i, 0] = np.cos(angle) * r + x
        points[i, 1] = np.sin(angle) * r + y
        points[i, 2] = z
    return points.T

def create_circles_from_strand(strand, n_mesh):
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
    r = (p1[-1] + p2[-1]) / 2.0
    dist = np.linalg.norm(p1[0:3] - p2[0:3])
    return np.pi * r * r * dist

def main():
    args = parse_arguments()
    file = args.file
    file_out = args.file_out
    g_ratio = args.g_ratio
    strand_id = args.strand_id

    if strand_id:
        strand_id = [int(i) for i in strand_id.split(",")]
    else:
        strand_id = None

    # Read strands
    strands, box_size = RWS.read_generic_list(file)

    max_rads_all = [np.max(v[:, 3]) for v in strands]
    max_rad = np.max(max_rads_all)
    xs = np.concatenate([v[:, 0] for v in strands])
    ys = np.concatenate([v[:, 1] for v in strands])
    zs = np.concatenate([v[:, 2] for v in strands])

    bounding_box = [
        np.array([np.min(xs), np.min(ys), np.min(zs)]) - max_rad,
        np.array([np.max(xs), np.max(ys), np.max(zs)]) + max_rad,
    ]
    print("Bounding box:", np.array(bounding_box))

    # Parameters
    n_circles = 8
    n_strands = len(strands)
    mesh = om.TriMesh()
    volumen_approx = 0

    if strand_id is None:
        my_range = range(int(n_strands))
    else:
        my_range = strand_id


    for i in progressbar.progressbar(my_range):
        strand = strands[i]

        if 0 <= args.decimate <=1: 
            n_percentage =len(strand)**(1-args.decimate)
            n_take = int(len(strand)//n_percentage)
            strand = strand[::n_take]

        if len(strand) <  4:
            sys.exit("Error: Strand is empty. decimation too strong")

        strand[:, 3] *= g_ratio


        if len(strand) == 1:
            center = strand[0][0:3]
            rad = strand[0][-1]
            sphere = pyvista.Sphere(radius=rad, center=center, theta_resolution=7, phi_resolution=7)
            mesh_sphere = sphere if mesh_sphere is None else mesh_sphere + sphere
            continue

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
        
        v1 = list_vert[ 0 ]
        for i in range(1,n_circles-1):
            v2 = list_vert[ i ]
            v3 = list_vert[ i+1 ]
            face1 = mesh.add_face(v1,v3,v2)
            
        v1 = list_vert[ -n_circles ]
        
        for i in range(-n_circles+1,-1):
            v2 = list_vert[ i ]
            v3 = list_vert[ i+1 ]
            #face1 = mesh.add_face(v3,v2,v1)
            face1 = mesh.add_face(v3,v1,v2)



    if strand_id is None:
        om.write_mesh(file_out, mesh)
    else:
        om.write_mesh(file_out, mesh)

if __name__ == "__main__":
    print("decimation ")
    main()
