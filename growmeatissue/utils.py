import model as m
import matplotlib.pyplot as plt
from typing import *
import numpy as np

class Counter:
    def __init__(self,
                 start:int = 0):
        self.cnt = start -1
    def __call__(self):
        self.cnt += 1
        return self.cnt 


def clean_ax(ax: plt.Axes)->None:
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])

    for sp in ax.spines.values():
        sp.set_visible(False)

def visualize_tissue(sc_genome: np.ndarray,
                     sc_crds: np.ndarray,
                     ):


    abber = np.sum(sc_genome != 1, axis=1)

    fig,ax  = plt.subplots(1,1,figsize=(10,10))

    sc = ax.scatter(sc_crds[:,0],
                    sc_crds[:,1],
                    c = abber,
                    cmap = plt.cm.plasma,
                    s = 2,
                    )

    ax.set_aspect("equal")

    clean_ax(ax)

    return fig,ax


def visualize_array(sc_genome,
                    sc_crd,
                    array_x: np.ndarray,
                    array_y: np.ndarray,
                    ):

    fig,ax = visualize_tissue(sc_genome,
                              sc_crd,
                              )

    ax.scatter(array_x,
               array_y,
               s = 30,
               c = "red",
               alpha = 0.3,
               )

    return fig,ax
