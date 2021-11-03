#!/usr/bin/env python3

import sys
import json

from typing import Iterable, Tuple, Dict, List
from pathlib import Path

import numpy as np

from skbio import DistanceMatrix
from skbio.tree import nj

from ete3 import Tree, TreeStyle, TextFace

#//Jaccard Index
#//                shared AB
#//  -------------------------------------
#//  exclusive A + shared AB + exclusive B

def read_names_file(names_file: Path) -> Dict[str, str]:
    assert names_file.exists()
    with names_file.open(mode="rt") as fhd:
        rows: List[str]           = fhd.readlines()
        cols: Iterable[List[str]] = (r.split("\t") for r in rows)
        names = {c[0].strip(): c[1].strip() for c in cols if len(c) == 2}
    return names

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
    total :np.ndarray = matrix[:,:,0:2].sum(axis=2).astype(np.float64)
    dist  :np.ndarray = 1.0 - (shared / (total-shared))
    # Note
    #   Jaccard sharedAB / (exclusiveA + exclusiveB + sharedAB)
    # can alse be written as
    #           sharedAB / (totalA + totalB - sharedAB)
    # with totalA = exclusiveA + sharedAB
    #      totalB = exclusiveB + sharedAB
    # therefore
    #      totalA + totalB - sharedAB
    #   =  exclusiveA + sharedAB + exclusiveB + sharedAB - sharedAB
    #   =  exclusiveA + sharedAB + exclusiveB

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

def cluster_distance(
        matrix_file              : Path      ,
        basefile                 : Path      ,
        distance                 : np.ndarray,
        names_file               : Path=None ,
        load_header              : bool=True ,
        save_matrix_redundant_tsv: bool=True ,
        save_matrix_redundant_np : bool=True ,
        save_matrix_condensed_tsv: bool=True ,
        save_matrix_condensed_np : bool=True ,
        save_tree_newick         : bool=True ,
        save_tree_ascii          : bool=True ,
        save_tree_png            : bool=True ,
    ) -> np.ndarray:
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

    if load_header:
        header_file = Path(f"{matrix_file}.json")
        with header_file.open(mode="rt") as fhd:
            header = json.load(fhd)

        project_name   = header["project_name"]
        num_samples    = len(header["data"])
        ids: List[str] = [d["header"]["input_file_name"] for d in header["data"]]

        assert num_samples == distance.shape[0]
    
    else:
        project_name   = matrix_file
        ids: List[str] = [str(d+1) for d in range(distance.shape[0])]

    if names_file:
        names = read_names_file(names_file)
        ids   = [names.get(i, i) for i in ids]


    dm   = DistanceMatrix(distance, ids)
    # print(dm, type(dm))
    # newick_tree, tree = nj(dm, result_constructor=lambda x: (x, TreeNode.read(StringIO(x), format='newick')))


    dmr: np.ndarray = dm.redundant_form()
    if save_matrix_redundant_np or save_matrix_redundant_tsv:
        if save_matrix_redundant_np:
            dmr_file = Path(f"{basefile}.mat.redundant.np")
            with dmr_file.open("wb") as fhd:
                np.save(fhd, dmr, allow_pickle=False)

        if save_matrix_redundant_tsv:
            dmr_file = Path(f"{basefile}.mat.redundant.lsmat")
            with dmr_file.open("wt") as fhd:
                dm.write(fhd, format='lsmat')


    if save_matrix_condensed_np or save_matrix_condensed_tsv:
        dmc: np.ndarray = dm.condensed_form()
        if save_matrix_condensed_np:
            dmc_file = Path(f"{basefile}.mat.condensed.np")
            with dmc_file.open("wb") as fhd:
                np.save(fhd, dmc, allow_pickle=False)

        if save_matrix_condensed_tsv:
            dmc_file = Path(f"{basefile}.mat.condensed.txt")
            with dmc_file.open("wt") as fhd:
                np.savetxt(fhd, dmc)


    if save_tree_newick or save_tree_ascii or save_tree_png:
        newick_tree:str = nj(dm, result_constructor=str)
        # tree: TreeNode = nj(dm)

        # treea = tree.to_array()
        # print(dmc)
        # print(dmr)
        # print(treea, type(treea))

        # print("cluster_distance :: tree",tree, type(tree))
        # print("cluster_distance :: tree.ascii_art",tree.ascii_art())
        # print("cluster_distance :: tree.to_taxonomy",tree.to_taxonomy())
        
        if save_tree_newick:
            newick_file = Path(f"{basefile}.newick")
            with newick_file.open("wt") as fhd:
                fhd.write(newick_tree)

        if save_tree_ascii or save_tree_png:
            ete_tree = Tree(newick_tree)

            if save_tree_ascii:
                tree_file = Path(f"{basefile}.tree")
                with tree_file.open("wt") as fhd:
                    fhd.write(str(ete_tree))

            if save_tree_png:
                png_file  = f"{basefile}.png"
                font_size = 12
                height    = font_size*4*(num_samples+5)
                width     = height // 2
                
                circular_style            = TreeStyle()
                circular_style.mode       = "c" # draw tree in circular mode
                circular_style.scale      = 20
                # leafy_style.show_branch_length = True
                # leafy_style.show_branch_support = True

                leafy_style               = TreeStyle()

                tree_style                = leafy_style
                tree_style.scale          = 60
                tree_style.show_leaf_name = True
                tree_style.title.add_face(TextFace(project_name, fsize=20), column=0)

                ete_tree.render(png_file, h=height, w=width, dpi=72, units="px", tree_style=tree_style)

    return dmr

def load(matrix_file: Path, names_file: Path=None):
    if names_file is None:
        names_file_t = Path(f"{matrix_file}.names.tsv")
        if names_file_t.exists():
            names_file = names_file_t

    matrix             = get_matrix(matrix_file)
    basefile, distance = calc_distance(matrix_file, matrix, fill_diagonal=True)
    redundant_matrix   = cluster_distance(matrix_file, basefile, distance, names_file=names_file)

def main():
    matrix_file = Path(sys.argv[1])
    load(matrix_file)

if __name__ == "__main__":
    main()