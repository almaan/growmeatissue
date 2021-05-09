#!/usr/bin/env python3
import numpy as np
import pandas as pd
import os.path as osp
import os

from typing import *
import model as m
import formating as f

import argparse as arp
import toml

import shutil
import logging

import utils as ut



def run(n_initial_cells,
        n_genes,
        action_rates,
        mutation_rate,
        transitions,
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
                      genome_length = int(genome_size),
                      sigma = sigma,
                      )

    logger.info("Initialize : Tissue Growth")
    population,lineage,cell_map = m.grow_tissue(n_initial_cells,
                                                genome = genome,
                                                action_rates = action_rates,
                                                mutation_rate = mutation_rate,
                                                transitions=transitions,
                                                max_iter = max_iter,
                                                max_cells = max_cells,
                                                domain_side_size = domain_side_size,
                                                init_coordinates = init_coordinates,
                                                init_event_size = init_event_size,
                                              )
    n_aberrant_cells = len(population)

    logger.info("Completed : Tissue growth | Total Cells : {}".format(n_aberrant_cells))

    logger.info("Intialized : Clustering | Clusters : {}".format(n_clusters))

    cluster_idx,cluster_genome = f.add_cluster(population,
                                               n_clusters,
                                             )


    cluster_lineage = f.get_cluster_lineage(population,
                                            cluster_idx,
                                            lineage)


    #TODO: save clustering relations
    logger.info("Completed : Clustering")

    logger.info("Initialized : Add normal cells")
    population,dropped = m.add_background(population,
                                          genome,
                                          cell_map,
                                          fraction_drop,
                                          )

    n_aberrant_cells -= len(dropped)
    n_normal_cells = len(population) - n_aberrant_cells
    
    cluster_idx = np.array([x for k,x in enumerate(cluster_idx) if k not in dropped])
    cluster_idx = np.append(cluster_idx,
                            (max(cluster_idx)+1)* np.ones(len(population)- len(cluster_idx)),
                            )

    logger.info("Completed : Add normal cells | Total Cells : {}".format(len(population)))
    logger.info("Initialized : Tissue Assembly")
    tissue_data = f.assemble_tissue_data(population)

    tissue_data["cluster_index"] = cluster_idx
    logger.info("Completed : Tissue Assembly")
    logger.info("Initialized : In silico Spatial Transcriptomics")
    st_data = f.in_silico_st(tissue_data,
                             domain_side_size,
                             depth = depth,
                             spot_radius=spot_radius,
                             n_spots_x = n_spots_x,
                             n_spots_y = n_spots_y,
                             )

    logger.info("Completed : In silico Spatial Transcriptomics")
    lineage = np.array([(c,p) for k,(c,p) in\
                        enumerate(lineage.items()) if\
                       k not in dropped])

    lineage = np.vstack((lineage,-1*np.ones((n_normal_cells,2))))


    # make genome data

    genome_data = f.format_genome(genome)

    res = dict(tissue_data = tissue_data,
               lineage = lineage,
               cluster_lineage = cluster_lineage,
               st_data = st_data,
               genome_data =genome_data,
               )

    return res


def main():

    prs = arp.ArgumentParser()

    aa = prs.add_argument

    aa("-d","--design_file")
    aa("-o","--out_dir")
    aa("-vo","--visual_output",default = False,action ="store_true")
    aa("-rs","--random_seed",default = None,required=False)

    args = prs.parse_args()

    if args.random_seed is not None:
        try:
            random_seed = int(args.random_seed)
            np.random.seed(random_seed)
        except:
            log.warning("Failed to set random seed. Will continue")


    ds = toml.load(args.design_file)

    if "transitions" in ds["population"].keys():
        tmat = np.zeros((2,2))
        tmat[1,0] = ds["population"]["transitions"]["s0s1"]
        tmat[0,1] = ds["population"]["transitions"]["s1s0"]
        tmat[0,0] = 1 - tmat[1,0]
        tmat[1,1] = 1 - tmat[0,1]
    else:
        tmat = None

    if "rates_state_0" in ds["population"].keys() and\
       "rates_state_1" in ds["population"].keys():

        rates_0 = ds["population"]["rates_state_0"]
        rates_1 = ds["population"]["rates_state_1"]

        action_rates = np.zeros((4,2))
        action_rates[:,0] = np.array((rates_0["stay_rate"],
                                      rates_0["spawn_rate"],
                                      rates_0["move_rate"],
                                      rates_0["death_rate"],
                                    ))

        action_rates[:,1] = np.array((rates_1["stay_rate"],
                                    rates_1["spawn_rate"],
                                    rates_1["move_rate"],
                                    rates_1["death_rate"],
                                    ))

        mutation_rate = np.array((rates_0["mutation_rate"],
                                rates_1["mutation_rate"]
                                ))

    else:
        action_rates = None
        mutation_rate = None

    if "initial_coordinates" in ds["population"].keys():
        crds = ds["population"]["initial_coordinates"]
        n = len(crds)
        initial_coordinates = [(c["x"],c["y"]) for c in crds]
    else:
        initial_coordinates = None


    dargs = dict(n_initial_cells = ds["population"]["n_initial_cells"],
                n_genes = ds["genome"]["n_genes"],
                n_clusters= ds["population"].get("n_clusters",None),
                genome_size = ds["genome"].get("genome_size",None),
                sigma = ds["genome"].get("event_spread",None),
                alpha = ds["genome"].get("concentration",None),
                max_iter = ds["population"].get("max_iterations",None),
                max_cells = ds["population"].get("max_cells",None),
                domain_side_size = ds["population"].get("domain_side_size"),
                init_event_size = ds["genome"].get("initial_event_size",None),
                fraction_drop = ds["population"].get("fraction_drop",None),
                depth = ds["spatial"].get("depth",None),
                spot_radius = ds["spatial"].get("spot_radius",None),
                n_spots_x = ds["spatial"].get("n_spots_x",None),
                n_spots_y = ds["spatial"].get("n_spots_x",None),
                )

    drop_dargs = [k for k,v in dargs.items() if v is None]
    for k in drop_dargs:
        dargs.pop(k)

    dargs["action_rates"] = action_rates
    dargs["mutation_rate"] = mutation_rate
    dargs["init_coordinates"]= initial_coordinates
    dargs["transitions"] = tmat


    res = run(**dargs)

    args.out_dir = osp.join(args.out_dir,"synthetic")

    if osp.exists(args.out_dir):
        shutil.rmtree(args.out_dir)

    os.makedirs(args.out_dir,exist_ok=True)


    if args.visual_output:
        visual_out_dir = osp.join(args.out_dir,"visual")
        os.mkdir(visual_out_dir)

    sc_out_dir = osp.join(args.out_dir,"single_cell_data")
    os.mkdir(sc_out_dir)

    sc_data = res["tissue_data"]

    if args.visual_output:
        fig,ax = ut.visualize_tissue(sc_data["genome_profile"],
                            sc_data["coordinates"],
                            )
        fig.savefig(osp.join(visual_out_dir,"tissue.png"))

    n_cells,n_genes = sc_data["expression"].shape
    cell_names = pd.Index(["Cell_{}".format(x) for x in range(n_cells)])
    gene_names = pd.Index(["Gene_{}".format(x) for x in range(n_genes)])

    sc_expr = pd.DataFrame(sc_data["expression"],
                           columns = gene_names,
                           index = cell_names,
                           )

    sc_genome_profile = pd.DataFrame(sc_data["genome_profile"],
                                     columns = gene_names,
                                     index = cell_names,
                                     )

    sc_meta_data = pd.DataFrame(np.hstack((sc_data["coordinates"],
                                           sc_data["cluster_index"][:,np.newaxis])),
                                columns = ["x","y","cluster"],
                                index = cell_names,
                                )

    sc_expr.to_csv(osp.join(sc_out_dir,"sc-expression.tsv"),sep="\t")
    sc_genome_profile.to_csv(osp.join(sc_out_dir,"sc-genome_profile.tsv"),sep="\t")
    sc_meta_data.to_csv(osp.join(sc_out_dir,"sc-meta_data.tsv"),sep="\t")

    st_out_dir = osp.join(args.out_dir,"spatial_data")
    os.mkdir(st_out_dir)
    st_data = res["st_data"]
    spot_names = pd.Index(["Spot_{}".format(x) for x in range(st_data["expression"].shape[0])])

    if args.visual_output:
        fig,ax = ut.visualize_array(sc_data["genome_profile"],
                           sc_data["coordinates"],
                           array_x = st_data["coordinates"][:,0],
                           array_y = st_data["coordinates"][:,1],
                           )
        fig.savefig(osp.join(visual_out_dir,"tissue_w_array.png"))


    st_expr = pd.DataFrame(st_data["expression"],
                           columns = gene_names,
                           index = spot_names,
                           )

    st_genome_profile = pd.DataFrame(st_data["genome_profile"],
                                     columns = gene_names,
                                     index = spot_names,
                                     )



    st_meta_data = pd.DataFrame(st_data["coordinates"],
                                columns = ["x","y"],
                                index = spot_names,
                                )

    st_meta_data["cell_count"] = st_data["cell_count"]


    if st_data["cluster_index"] is not None:
        cluster_names = pd.Index(["Cluster_{}".format(x) for x \
                         in range(st_data["cluster_index"].shape[1])])

        cluster_data = pd.DataFrame(st_data["cluster_index"],
                                    columns = cluster_names,
                                    index = spot_names,
                                    )
        cluster_data.to_csv(osp.join(args.out_dir,"cluster_data.tsv"),sep="\t")


    st_annotation = pd.DataFrame(st_data["spot_status"],
                                 index = spot_names,
                                 columns = ["annotation"],
                                 )

    st_annotation.annotation = st_annotation.annotation.map({1:"benign",0:"aberrant"})

    for df,name in zip([st_expr.T,st_genome_profile.T,st_meta_data,st_annotation],
                  ["st-expression","st-genome_profile","st-meta","st-annotation"]
                  ):

        df.to_csv(osp.join(st_out_dir,name + ".tsv"),sep="\t")

    with open(osp.join(st_out_dir,"cell_by_spot.tsv"),"w+") as f:
        for k,s in enumerate( st_data["cell_by_spot"] ):
            if not isinstance(s,list):
                spot_list = [s]
            else:
                spot_list = s

            spot_list = [str(x) for x in spot_list]
            f.write(f"Spot_{k} : " + ",".join(spot_list) + "\n")


    #TODO: could drop eliminate cluster?
    cluster_names = pd.Index(["Cluster_{}".format(x) for x in range(res["cluster_lineage"].shape[0])])
    cluster_lineage = pd.DataFrame(res["cluster_lineage"],
                                   index=cluster_names,
                                   columns = cluster_names,
                                   )

    lineage = pd.DataFrame(res["lineage"],
                           columns = ["child","parent"],
                           index = cell_names,
                           )

    genome_data = res["genome_data"]
    genome_data.index = gene_names

    for df,name in zip([cluster_lineage,lineage,genome_data],
                       ["cluster_lineage","lineage","genome"],
                       ):

        df.to_csv(osp.join(args.out_dir,name + ".tsv"),sep="\t")


if __name__ == "__main__":
    FORMAT = '[%(name)s | %(asctime)s |  %(levelname)s] >> %(message)s'
    logging.basicConfig(format=FORMAT)
    logger = logging.getLogger('synthdata')
    logger.setLevel(logging.DEBUG)

    main()
