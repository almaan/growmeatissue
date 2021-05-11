from __future__ import annotations

import copy
from typing import Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import utils as ut
from matplotlib import rcParams
from scipy.stats import truncnorm


class Genome:
    def __init__(
        self,
        theta,
        genome_length: int,
        expr_ub=1000,
        expr_lb=100,
        sigma: int = 2,
    ) -> None:
        self.theta = theta
        self.G = theta.size

        self.gene_id = np.arange(self.G)
        self.pos = np.zeros(self.G)

        self.genome_length = genome_length
        self._build()

        self.sigma_raw = sigma
        self.sigma = self.sigma_raw * np.diff(self.pos).mean()
        self.len_dist = truncnorm(
            0, np.inf, loc=np.diff(self.pos).mean(), scale=self.sigma
        )

        self.expr_ub = expr_ub
        self.expr_lb = expr_lb

    def _build(
        self,
    ) -> None:
        pos = np.random.random(self.G)
        pos[0] = 0.0
        pos = np.cumsum(pos)
        pos = pos / pos[-1]
        self.pos = pos * self.genome_length

    def in_range(
        self,
        center: float,
        width: float,
    ) -> np.ndarray:
        return np.where(np.abs((self.pos - center)) < width)[0]

    def add(
        self,
        idx: np.ndarray,
    ) -> None:

        if len(idx) > 0:
            start_pos = idx[0]
            end_pos = idx[-1]

            amp_lens = self.pos[start_pos : end_pos + 1]
            delta = amp_lens[-1] - amp_lens[0]

            p1_pos = self.pos[0 : end_pos + 1]
            p2_pos = amp_lens + delta
            p3_pos = self.pos[end_pos + 1 : :] + delta
            pos = np.hstack((p1_pos, p2_pos, p3_pos))

            p1_gen = self.gene_id[0 : end_pos + 1]
            p2_gen = self.gene_id[idx]
            p3_gen = self.gene_id[end_pos + 1 : :]
            gen_id = np.hstack((p1_gen, p2_gen, p3_gen))

            self.pos = pos
            self.gene_id = gen_id
            self.genome_length = self.pos[-1]

    def delete(
        self,
        idx: np.ndarray,
    ) -> None:

        if len(idx) > 0:
            start_pos = idx[0]
            end_pos = idx[-1]

            amp_lens = self.pos[start_pos : end_pos + 1]
            delta = amp_lens[-1] - amp_lens[0]
            p1_pos = self.pos[0:start_pos]
            p2_pos = self.pos[end_pos + 1 : :] - delta
            self.pos = np.hstack((p1_pos, p2_pos))

            self.gene_id = np.hstack(
                (self.gene_id[0:start_pos], self.gene_id[end_pos + 1 : :])
            )
            self.genome_length = self.pos[-1]

    def sample_expression(
        self,
    ) -> np.ndarray:
        scale_factor = self.genome_length / self.G
        lb = self.expr_lb * scale_factor
        ub = self.expr_ub * scale_factor
        n_genes = np.random.randint(lb, ub)
        pvals = self.theta[self.gene_id]
        pvals /= pvals.sum()
        x_raw = np.random.multinomial(n_genes, pvals)
        x_cur = np.zeros(self.G)
        for k, g in enumerate(self.gene_id):
            x_cur[g] += x_raw[k]
        return x_cur

    def sample_event(
        self,
    ) -> np.ndarray:
        center = np.random.uniform(0, self.genome_length)
        width = self.len_dist.rvs()
        return (center, width)

    def get_gene_count(
        self,
    ) -> np.ndarray:
        counts = np.zeros(self.G)
        for g in self.gene_id:
            counts[g] += 1
        return counts


class Cell:
    def __init__(
        self,
        genome,
        mutation_rate: Optional[np.ndarray],
        action_rates: Optional[np.ndarray],
        counter,
        x=0,
        y=0,
        active_state=5,
        p_move=None,
        tmat=None,
        parent=-1,
    ) -> None:

        self.counter = counter
        self.id = self.counter()
        self.parent = parent
        self.children = []

        self.genome = genome

        self.N = self.genome.G

        if mutation_rate is not None:
            self.mutation_rate = mutation_rate
        else:
            self.mutation_rate = np.random.uniform(0, 1, size=2)

        if action_rates is not None:
            self.action_rates = action_rates

            self.action_rates /= self.action_rates.sum(axis=0, keepdims=True)
        else:
            self.action_rates = np.random.dirichlet(np.ones(4), size=2).T

        if tmat is None:
            self.tmat = np.array([[0.8, 0.2], [0.4, 0.6]]).T
        else:
            self.tmat = tmat

        self.active = active_state

        assert np.all(self.action_rates >= 0), "rates must be non-negative than zero"
        self.crd = np.array([x, y])
        if p_move is None:
            self.p_move = np.random.dirichlet(0.1 * np.ones(4))
        else:
            self.p_move = p_move

        self.gamma = 0.05

        self.dx = [-1, 0, 1, 0]
        self.dy = [0, -1, 0, 1]
        self.step_scale = [1, 2]

    def search_nbrhd(
        self,
        cell_map,
    ) -> Tuple[np.ndarray, np.ndarray]:

        free_space = []
        pvals = []
        scale_x = np.random.choice(self.step_scale)
        scale_y = np.random.choice(self.step_scale)
        for k, (dx, dy) in enumerate(zip(self.dx, self.dy)):
            new_x = self.crd[0] + dx * scale_x
            new_y = self.crd[1] + dy * scale_y

            if (new_x < 0) | (new_x >= cell_map.shape[0]):
                break
            elif (new_y < 0) | (new_y >= cell_map.shape[1]):
                break

            if cell_map[new_x, new_y] == 0:
                free_space.append((new_x, new_y))
                pvals.append(self.p_move[k])

        pvals = [x / sum(pvals) for x in pvals]
        return free_space, pvals

    def genome_update(
        self,
    ) -> None:
        # check if mutation should occur
        if np.random.random() < self.mutation_rate[self.state]:
            self.mutation_rate[self.state] *= 0.95
            center, length = self.genome.sample_event()
            idx = self.genome.in_range(center, length)

            if np.random.random() < 0.5:
                self.genome.add(idx)
            else:
                self.genome.delete(idx)

    def sample_expression(
        self,
    ) -> np.ndarray:
        return self.genome.sample_expression()

    @property
    def state(
        self,
    ) -> int:
        return int(self.active < 0)

    def action(
        self,
        cell_map: np.ndarray,
    ) -> Union[List[None], List[Cell]]:

        if self.active < 0:
            new_state = np.random.multinomial(1, self.tmat[:, 1]).argmax()
            if new_state == 0:
                self.active = np.random.geometric(p=self.tmat[1, 0])
        else:
            self.active -= 1

        action = np.random.multinomial(1, self.action_rates[:, self.state]).argmax()
        if action == 3 and cell_map.sum() > 50:
            cell_map[self.crd[0], self.crd[1]] = 0
            return []
        else:
            free_space, pvals = self.search_nbrhd(cell_map)
            if len(free_space) > 0:
                new_crd_idx = np.random.choice(len(free_space), p=pvals)
                new_crd = np.array(free_space[new_crd_idx])
                if action == 1:
                    self.genome_update()
                    new_genome = copy.deepcopy(self.genome)
                    new_cell = Cell(
                        genome=new_genome,
                        mutation_rate=self.mutation_rate,
                        action_rates=self.action_rates,
                        counter=self.counter,
                        x=new_crd[0],
                        y=new_crd[1],
                        active_state=self.active,
                        p_move=(self.p_move if self.active >= 0 else None),
                        tmat=self.tmat,
                        parent=self.id,
                    )
                    self.children.append(new_cell.id)

                    cell_map[new_crd[0], new_crd[1]] = 1

                    return [self, new_cell]

                elif action == 2:
                    cell_map[self.crd[0], self.crd[1]] = 0
                    self.crd = new_crd
                    cell_map[self.crd[0], self.crd[1]] = 1
                    return [self]
                else:
                    return [self]
            else:
                return [self]


def grow_tissue(
    n_initial_cells: int,
    genome: Genome,
    action_rates: Optional[np.ndarray] = None,
    mutation_rate: Optional[np.ndarray] = None,
    transitions: Optional[np.ndarray] = None,
    max_iter: int = 500,
    max_cells: Union[int, float] = 1e4,
    domain_side_size: float = 200,
    init_coordinates: Optional[np.ndarray] = None,
    init_event_size: Optional[int] = None,
) -> Tuple[List[Cell], Dict[int, int], np.ndarray]:

    n_genes = genome.G

    if action_rates is None:
        action_rates = np.array(((0.1, 0.3), (0.4, 0.45), (0.4, 0.2), (0.0, 0.00)))
    else:
        action_rates = action_rates

    if transitions is None:
        tmat = np.array([[0.8, 0.2], [0.2, 0.8]]).T
    else:
        tmat = transitions

    if mutation_rate is None:
        mutation_rate = np.array([0.0, 0.8])
    else:
        mutation_rate = mutation_rate

    if init_coordinates is None:
        init_coordinates = []
        for i in range(n_initial_cells):
            xy = np.random.uniform(
                0.1 * domain_side_size,
                0.9 * domain_side_size,
                size=2,
            )
            init_coordinates.append(xy)
    else:
        if isinstance(init_coordinates[0], int) or isinstance(
            init_coordinates[0], float
        ):
            init_coordinates = [init_coordinates]

            assert (
                len(init_coordinates) == n_initial_cells
            ), "incorrect number of initial coordinate tuples"
            assert len(init_coordinates[0] == 2), "coordinate _tuple_ is required"

            for ic in init_coordinates:
                for ii in range(len(ic)):
                    assert (
                        ic[ii] > 0 & ic[ii] < domain_side_size
                    ), "Initial coordinate pair ({},{}) is out of bounds".format(
                        ic[0], ic[1]
                    )

    if init_event_size is None:
        init_event_size = max(int(0.03 * n_genes), 1)

    cell_map = np.zeros((domain_side_size, domain_side_size))
    counter = ut.Counter()
    lineage = {}
    population = []

    for i in range(n_initial_cells):

        seed_cell = Cell(
            genome=copy.deepcopy(genome),
            counter=counter,
            mutation_rate=mutation_rate,
            action_rates=action_rates,
            tmat=tmat,
            x=init_coordinates[i][0],
            y=init_coordinates[i][1],
        )

        start_event = np.random.randint(0, n_genes - init_event_size - 1)
        init_event_genes = np.arange(start_event, start_event + init_event_size)
        if i % 2 == 0:
            seed_cell.genome.add(init_event_genes)
        else:
            seed_cell.genome.delete(init_event_genes)

        lineage.update({seed_cell.parent: [seed_cell.id]})
        population.append(seed_cell)

    print()
    for t in range(max_iter):
        print("\r Iteration : {} | Cells : {}".format(t, len(population)), end="")
        new_cells = []
        for k in range(len(population)):
            cell = population.pop(0)
            out = cell.action(cell_map)
            for c in out:
                if c.id not in lineage:
                    lineage.update({c.id: c.parent})

            new_cells += out

        population += new_cells
        if cell_map.sum() > max_cells:
            break
    print()
    lineage.pop(-1)
    return population, lineage, cell_map


def add_background(
    population: List[Cell],
    genome: Genome,
    cell_map: np.ndarray,
    fraction_drop: float = 0,
) -> Tuple[List[Cell], np.ndarray]:

    pop_size = len(population)
    if fraction_drop > 0:
        drop = np.random.choice(
            pop_size, replace=False, size=int(pop_size * fraction_drop)
        )

        new_cell_map = copy.deepcopy(cell_map)
        for i in drop:
            new_cell_map[int(population[i].crd[0]), int(population[i].crd[1])] -= 1

        population = [c for k, c in enumerate(population) if k not in drop]

    counter = ut.Counter(max([x.id for x in population]))
    normal_population = []
    row, col = np.where(new_cell_map == 0)
    for x, y in zip(row, col):
        normal_population.append(
            Cell(
                x=x,
                y=y,
                genome=copy.deepcopy(genome),
                counter=counter,
                action_rates=None,
                mutation_rate=None,
            )
        )

    population = population + normal_population

    return population, drop
