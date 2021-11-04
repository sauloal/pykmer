#!/usr/bin/env python3

import sys
import json

from typing import Iterable, Tuple, Dict, List
from pathlib import Path

import numpy as np
import pandas as pd

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

def get_matrix(matrix_file: Path) -> Tuple[np.ndarray, List[str]]:
    assert matrix_file.exists()
    assert matrix_file.is_file()
    print(f"get_matrix :: matrix_file {matrix_file}")
    # matrix = np.loadtxt(matrix_file, dtype=np.float64, delimiter="\t")
    # matrix = np.loadtxt(matrix_file, delimiter="\t", skiprows=1)
    matrix = pd.read_csv(matrix_file, sep="\t", header=0, index_col=0)
    names: List[str] = list(matrix.columns)
    matrix = matrix.astype(np.float64).to_numpy()
    print(f"matrix", matrix)
    print(f"names", names)
    print(f"get_matrix :: matrix type", type(matrix))
    print(f"get_matrix :: matrix.dtype", matrix.dtype)
    print(f"get_matrix :: matrix.shape", matrix.shape)
    return matrix, names

def cluster_distance(
        matrix_file              : Path      ,
        basefile                 : Path      ,
        distance                 : np.ndarray,
        names                    : List[str] ,
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

    if names_file:
        names_ids = read_names_file(names_file)
        names     = [names.get(i, i) for i in names_ids]

    num_samples = len(names)
    project_name = matrix_file

    dm   = DistanceMatrix(distance, names)
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

    distance, names  = get_matrix(matrix_file)
    redundant_matrix = cluster_distance(matrix_file, matrix_file, distance, names=names, names_file=names_file)

def main():
    matrix_file = Path(sys.argv[1])
    load(matrix_file)

if __name__ == "__main__":
    main()