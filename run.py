import numpy as np
import pandas as pd
import os.path as osp

import models as m
import formating as f



def run(n_initial_cells,
        n_genes,
        action_rates,
        mutation_rate,
        n_clusters: int = 10,
        genome_size: int = 1e5,
        sigma: float = 10,
        alpha: float = 1,
        max_iter: int = 500,
        max_cells: Union[int,float] = 1e4,
        domain_side_size: float = 200,
        init_coordinates: Optional[np.ndarray] = None,
        init_event_size: Optional[int] = None,
        fraction_drop: float = 0.05,
        depth: int = 5000,
        spot_radius:float = 2.5,
        n_spots_x: int = 10,
        n_spots_y: int = 10,
        )->Dict[str,np.ndarray]:



    alpha = np.ones(n_genes) * alpha
    theta = np.random.dirichlet(alpha)
    genome = m.Genome(theta,
                      genome_size = int(genome_size),
                      sigma = sigma,
                      )
    population,lineage,cell_map = m.grow_tissue(n_initial_cells,
                                              genome,
                                              action_rates,
                                              mutation_rate,
                                              max_iter,
                                              max_cells,
                                              domain_side_size,
                                              init_coordinates,
                                              init_event_size,
                                              )

    #TODO: save lineage


    cluster_idx,cluster_genome = f.add_cluster(population,
                                               n_clusters,
                                             )



    cluster_relations = f.get_cluster_lineage(population,
                                            cluster_idx,
                                            lineage)


    #TODO: save clustering relations

    population,dropped = m.add_background(population,
                                          cell_map,
                                          fraction_drop,
                                          )

    cluster_idx = np.array([x for k,x in enumerate(cluster_idx) if k not in dropped])
    cluster_idx = np.append(cluster_idx,
                            (max(cluster_idx)+1)* np.ones(len(cluster_idx)-len(population)),
                            )

    tissue_data = f.assemble_tissue_data(population)

    tissue_data["cluster_index"] = cluster_idx

    st_data = f.in_silico_st(tissue_data,
                             domain_side_size,
                             depth = depth,
                             spot_radius=spot_radius,
                             n_spots_x,
                             n_spots_y,
                             )


    lineage = np.array(( (c,p) for c,p in lineage.items()))

    res = dict(tissue_data = tissue_data,
               lineage = lineage,
               st_data = st_data,
               )

    return res


def main():


    run(*args)
