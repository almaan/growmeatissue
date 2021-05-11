from typing import Dict, List, Tuple

import model as m
import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering


def add_cluster(
    population: List[m.Cell],
    n_clusters: int,
) -> Tuple[np.ndarray, np.ndarray]:

    genome_profile = np.vstack(
        [x.genome.get_gene_count()[np.newaxis, :] for x in population]
    )
    cluster = AgglomerativeClustering(
        n_clusters=n_clusters, affinity="euclidean", linkage="ward"
    )

    idx = cluster.fit_predict(genome_profile)
    cluster_genome = np.zeros((n_clusters, population[0].genome.G))

    for clu in np.unique(idx):
        pos = idx == clu
        cluster_genome[clu, :] = genome_profile[pos, :].mean(axis=0).round()

    return idx, cluster_genome


def get_cluster_lineage(
    population: List[m.Cell],
    cluster_idx: np.ndarray,
    lineage: Dict[int, int],
) -> np.ndarray:

    n_cluster = np.unique(cluster_idx).shape[0]
    pos_to_id = {x: population[x].id for x in range(len(population))}
    id_to_pos = {v: k for k, v in pos_to_id.items()}
    cmat = np.zeros((n_cluster, n_cluster))

    for k1 in range(n_cluster):
        c_members = np.where(cluster_idx == k1)[0]
        for mb in c_members:
            i = pos_to_id[mb]
            k2 = k1
            while (k2 == k1) and (i > 0):
                i = lineage[i]
                if i in id_to_pos.keys():
                    k2 = cluster_idx[id_to_pos[i]]

            cmat[k1, k2] += 1

    return cmat


def assemble_tissue_data(
    population: List[m.Cell],
) -> Dict[str, np.ndarray]:

    expression = np.array([x.sample_expression() for x in population])
    crd = np.array([x.crd for x in population])
    genome = np.vstack([x.genome.get_gene_count()[np.newaxis, :] for x in population])
    benign = np.array([int(np.all(x.genome.get_gene_count() == 1)) for x in population])

    res = dict(
        expression=expression,
        coordinates=crd,
        genome_profile=genome,
        cell_status=benign,
    )

    return res


def in_silico_st(
    tissue_data: Dict[str, np.ndarray],
    domain_side_size: float,
    depth: int = 5000,
    spot_radius: float = 2.5,
    n_spots_x: int = 10,
    n_spots_y: int = 10,
) -> Dict[str, np.ndarray]:

    radius_squared = spot_radius ** 2

    xx = np.linspace(spot_radius, domain_side_size - spot_radius, n_spots_x)

    yy = np.linspace(spot_radius, domain_side_size - spot_radius, n_spots_y)

    xx, yy = np.meshgrid(xx, yy)
    xx = xx.flatten()
    yy = yy.flatten()

    n_genes = tissue_data["expression"].shape[1]

    exp_matrix = np.zeros((xx.shape[0], n_genes))
    gen_matrix = np.zeros((xx.shape[0], n_genes))
    ben_matrix = np.zeros(xx.shape[0])
    cell_count = np.zeros(xx.shape[0])
    cell_id_list = []

    make_cluster_matrix = "cluster_index" in tissue_data.keys()
    if make_cluster_matrix:
        n_clusters = np.unique(tissue_data["cluster_index"]).shape[0]
        cluster_matrix = np.zeros((xx.shape[0], n_clusters))

    for spot, (x, y) in enumerate(zip(xx, yy)):

        in_spot = []
        for k in range(tissue_data["coordinates"].shape[0]):
            delta = (tissue_data["coordinates"][k, 0] - x) ** 2 + (
                tissue_data["coordinates"][k, 1] - y
            ) ** 2
            if delta < radius_squared:
                in_spot.append(k)

        in_spot = np.array(in_spot)

        cell_id_list.append(in_spot.tolist())

        joint_expr = tissue_data["expression"][in_spot, :].sum(axis=0)
        gene_prob = joint_expr / joint_expr.sum()
        ns = np.min((depth, joint_expr.sum()))
        exp_matrix[spot, :] = np.random.multinomial(ns, gene_prob)
        cell_count[spot] = len(in_spot)

        av_gen = tissue_data["genome_profile"][in_spot].mean(axis=0)
        gen_matrix[spot, :] = av_gen

        state_ben = int(tissue_data["cell_status"][in_spot].mean() >= 0.9)
        ben_matrix[spot] = state_ben

        if make_cluster_matrix:
            clu, cnt = np.unique(tissue_data["cluster_index"], return_counts=True)
            cluster_matrix[spot, clu.astype(int)] = cnt / cnt.sum()

    crd = np.hstack((xx[:, np.newaxis], yy[:, np.newaxis]))

    res = dict(
        expression=exp_matrix,
        coordinates=crd,
        genome_profile=gen_matrix,
        spot_status=ben_matrix,
        cluster_index=(cluster_matrix if make_cluster_matrix else None),
        cell_count=cell_count,
        cell_by_spot=cell_id_list,
    )

    return res


def format_genome(
    genome: m.Genome,
) -> pd.DataFrame:

    genome_start = genome.pos
    len_last_gene = int(np.diff(genome_start).mean())
    genome_end = np.append(genome.pos[1::], genome.pos[-1] + len_last_gene)
    chrom = np.array(["chr1" for x in range(genome_start.shape[0])])
    tmp = np.hstack(
        (chrom[:, np.newaxis], genome_start[:, np.newaxis], genome_end[:, np.newaxis])
    )

    gene_data = pd.DataFrame(
        tmp,
        columns=["chr", "start", "end"],
    )

    return gene_data
