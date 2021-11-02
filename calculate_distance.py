#!/usr/bin/env python3

import os
import sys
import math
import json

from typing import Tuple
from io import StringIO
from pathlib import Path

import numpy as np

import skbio
from skbio import DistanceMatrix
from skbio.tree import nj
from skbio.tree import TreeNode
from ete3 import Tree, TreeStyle, TextFace

#//Jaccard Index
#//                shared AB
#//  -------------------------------------
#//  exclusive A + shared AB + exclusive B

def get_matrix(matrix_file: Path) -> np.ndarray:
    assert matrix_file.exists()
    assert matrix_file.is_file()
    print(f"get_matrix :: matrix_file {matrix_file}")
    npz = np.load(matrix_file)
    assert "matrix" in npz
    matrix = npz["matrix"]
    # print(f"matrix", matrix)
    print(f"get_matrix :: matrix type", type(matrix))
    print(f"get_matrix :: matrix.dtype", matrix.dtype)
    print(f"get_matrix :: matrix.shape", matrix.shape)
    return matrix

def calc_distance(matrix_file: Path, matrix: np.ndarray, fill_diagonal: bool=True) -> Tuple[Path,np.ndarray]:
    print(f"calc_distance :: matrix.shape", matrix.shape)
    # matrix = matrix[:3,:3,:]
    
    # shape = matrix.shape
    # res = np.ndarray(shape=(shape[0],shape[0]), dtype=np.float64)
    # for k in range(shape[0]):
    #     # print("calc_distance :: k", k)
    #     for l in range(k+1, shape[0]):
    #         # print("calc_distance ::  l", l)
    #         (k_count, l_count, s_count) = matrix[k,l,:]
    #         (k_count, l_count, s_count) = (k_count.item(), l_count.item(), s_count.item())
    #         # dists = Distances(k_count, l_count, s_count)
    #         # dist  = dists.D_jaccard()
    #         total = (k_count + l_count + s_count)
    #         dist  = 1 - (s_count / total)
            
    #         val   = matrix[k,l,:]
    #         # print(val)
    #         # print(val.shape)
    #         # print(val.dtype)
    #         dista = 1 - (val[2] / val.sum())

    #         res[k,l] = dist
    #         res[l,k] = dist
    #         print("calc_distance ::   k_count", k, k_count)
    #         print("calc_distance ::   l_count", l, l_count)
    #         print("calc_distance ::   s_count", '-', s_count)
    #         print("calc_distance ::   total  ", '-', total)
    #         print("calc_distance ::   dist   ", '-', dist)
    #         print("calc_distance ::   dista  ", '-', dista)
    #         # sys.exit(0)

    # print("calc_distance :: res.shape", res.shape)
    # print("calc_distance :: res.dtype", res.dtype)
    # # print("calc_distance :: res", res)


    # print("res   ", res   , res   .shape, res.dtype)

    shared:np.ndarray = matrix[:,:,2].astype(np.float64)
    total :np.ndarray = matrix.sum(axis=2).astype(np.float64)
    dist  :np.ndarray = 1.0 - (shared / total)

    if fill_diagonal:
        np.fill_diagonal(dist, 0.0)

    # print("calc_distance :: shared", shared, shared.shape, shared.dtype)
    # print("calc_distance :: total ", total , total .shape, total.dtype)
    # print("calc_distance :: dist  ", dist  , dist  .shape, dist.dtype)
    print("calc_distance :: dist  ", dist  .shape, dist.dtype)

    basefile    = Path(f"{matrix_file}.dist.jaccard")
    matrix_file = Path(f"{basefile}.npz")
    with matrix_file.open(mode="wb") as fhd:
        np.savez(fhd, distance=dist)

    return basefile, dist

def cluster_distance(matrix_file: Path, basefile: Path, distance: np.ndarray) -> np.ndarray:
    # http://scikit-bio.org/docs/0.2.1/generated/skbio.tree.nj.html
    # http://scikit-bio.org/docs/0.2.1/generated/skbio.tree.TreeNode.html?highlight=treenode
    # http://scikit-bio.org/docs/0.2.1/generated/generated/skbio.stats.distance.DistanceMatrix.html?highlight=distancematrix
    # http://etetoolkit.org/docs/latest/tutorial/tutorial_drawing.html#the-programmable-tree-drawing-engine
    #
    # data = [[0,  5,  9,  9,  8],
    #         [5,  0, 10, 10,  9],
    #         [9, 10,  0,  8,  7],
    #         [9, 10,  8,  0,  3],
    #         [8,  9,  7,  3,  0]]
    # ids = list('abcde')

    header_file = Path(f"{matrix_file}.json")
    with header_file.open(mode="rt") as fhd:
        header = json.load(fhd)

    project_name = header["project_name"]
    num_samples  = len(header["data"])
    # ids          = [d["header"]["input_file_name"] for d in header["data"]]
    ids          = [d["header"]["input_file_name"] for d in header["data"]]

    assert num_samples == distance.shape[0]
    dm   = DistanceMatrix(distance, ids)
    # print(dm, type(dm))
    # newick_tree, tree = nj(dm, result_constructor=lambda x: (x, TreeNode.read(StringIO(x), format='newick')))
    newick_tree:str = nj(dm, result_constructor=str)
    # tree: TreeNode = nj(dm)
    # dmc: np.ndarray = dm.condensed_form()
    dmr: np.ndarray = dm.redundant_form()
    # treea = tree.to_array()
    # print(dmc)
    # print(dmr)
    # print(treea, type(treea))

    # print("cluster_distance :: tree",tree, type(tree))
    # print("cluster_distance :: tree.ascii_art",tree.ascii_art())
    # print("cluster_distance :: tree.to_taxonomy",tree.to_taxonomy())
    
    newick_file = Path(f"{basefile}.newick")
    with newick_file.open("wt") as fhd:
        fhd.write(newick_tree)

    dm_file = Path(f"{basefile}.lsmat")
    with dm_file.open("wt") as fhd:
        dm.write(fhd, format='lsmat')

    dmr_file = Path(f"{basefile}.mat.npz")
    with dmr_file.open("wb") as fhd:
        np.savez(dmr_file, redundant_matrix=dmr)

    ete_tree = Tree(newick_tree)

    tree_file = Path(f"{basefile}.tree")
    with tree_file.open("wt") as fhd:
        fhd.write(str(ete_tree))


    png_file  = f"{basefile}.png"
    font_size = 12
    height    = font_size*4*(num_samples+5)
    width     = height // 2
    
    circular_style = TreeStyle()
    circular_style.mode = "c" # draw tree in circular mode
    circular_style.scale = 20

    leafy_style = TreeStyle()
    # leafy_style.show_branch_length = True
    # leafy_style.show_branch_support = True

    scale_style = TreeStyle()
    # scale_style.scale =  120 # 120 pixels per branch length unit

    tree_style = leafy_style
    tree_style.scale = 60
    tree_style.show_leaf_name = True
    tree_style.title.add_face(TextFace(project_name, fsize=20), column=0)

    ete_tree.render(png_file, h=height, w=width, dpi=72, units="px", tree_style=tree_style)

def load(matrix_file: Path):
    matrix             = get_matrix(matrix_file)
    basefile, distance = calc_distance(matrix_file, matrix, fill_diagonal=True)
    cluster            = cluster_distance(matrix_file, basefile, distance)

def main():
    matrix_file = Path(sys.argv[1])
    load(matrix_file)

if __name__ == "__main__":
    main()