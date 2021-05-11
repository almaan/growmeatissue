# Description
This repository contains the code for a generative process that, with several
important intermediary steps, will generate _in silico_ spatial
transcriptomics data **with associated - spatially aware - clonal (genomic) information**. In short, the
process can be described as:

1. Define a tissue domain (which cell will populate)
2.  Initialize the tissue with seeding cancer cells (these will have an amplification or deletion event)
3.  Grow the tissue (during this process, cells multiply and mutation events can occur)
4.  Add "normal" cells to regions that cancer cells do not populate
5. Place a grid of spots over the tissue
6. Identify which cells that overlay respective spot in the grid
7. Sample gene expression from cells (this is dependent on both the standard
   expression levels and the number of copies of each cell), and "capture" a
   fraction of these transcripts at the associated spot
  
Below, the model is described in more detail:

# Model
The model itself consists of two main components: the genome and the cell. Each
cell has a genome which dictates its expression profile, and each cell can also
perform certain actions. Once these two components are defined, we can apply a
set of functions to generate the _in silico_ data. Below, these components and
actions are described in a sequential manner.

## Genome
We will define the initial genome structure, that will the basis for every cells
individual genome. For the sake of simplicity, we our genome will consists of a
single chromosome with no introns. The genome has certain attributes associated with it:

- length : the total length (in base pairs) of the genome
- genes : elements in the genome, each gene has a
  - start position : where in the genome the gene starts
  - end position : where in the genome the gene ends
  - expression tendency : this is referred to as $theta$-value, and determines
    the expression level of each gene.

The genes have different sizes, and their lengths are randomly initialized such
that their average length is `n_genes / genome_length`.

### Amplification and deletion events
The Genome has two actions associated with it: _deletion_ and _addition_. Both
these events occur in a similar fashion (even though the results are quite the
opposite). When either of the actions is called, we will sample a random
position on the genome which will figure as the center of the event, and then
determine the width of the event, the latter is sampled from a truncated
normal(at zero) with zero mean and the variable `event_spread` set as the
standard deviation value. The genes that are deleted/amplified are all genes
encompassed in the region `[center - width,center+width]`. This means that we
will see co-amplification/deletion of genes, if they reside near each other on
the genome.


### Expression
We will model gene expression by a Dirichlet-Multinomial process defined as follows:

![expression_model](img/expression-model.png)

Where _x_<sub>c</sub> is the expression vector for cell _c_ and _n_<sub>cg</sub> is
the number of copies of gene _g_ that cell _c_ has. _N_<sub>c</sub> is the number of
transcripts to sample from the specific cell.

This means that the propensity to express a certain gene is partially dependent
on its baseline level (`theta`) but also eventual changes to the genome.


## Cell

Each cell will have a genome (as described above associated with it), in a
population of cells all inhabiting the same domain (tissue), the genome of all
cells will originate from the same initial genome (meaning they all have the
same `theta` values), but they might have picked up individual changes along the
course of time (e.g., amplification and deletion events.) which means that their
_n_<sub>cg</sub> values and hence `theta_hat` values differ.

Every cell object has four main actions that it can perform:

- die : if a cell dies it is removed from the population.
- move : a cell can move either one or two steps at each life cycle. The
  propensity to move in a certain direction varies between each individual cell.
  A cell can only move to a free location, if no free space is found in the
  neighborhood it will remain stagnant.
- spawn : a cell can also spawn offspring, a copy of itself with the exception
  of _eventual_ mutations to the genome. Just as for the move, a cell can only
  spawn offspring to a free space in its neighborhood, if no free space is
  available the cell will not spawn. Find more on mutations below.
- rest : a cell that rests remains in place.

The tendency for each of these actions to occur is specified with a sate of
`action_rates`, which should all sum to 1. Meaning that when we
call for a cell to perform an action, we will do by: `action ~
Categorical(p_die,p_move,p_spawn,p_rest)`. Modifying this will have a
significant impact on the population's character.


### Mutations
Mutations in the genome might occur, to the daughter cell, when a cell is
spawning. Mutations do not occur during any of the other actions. When mutating
it's a 1:1 chance of an deletion vs amplification to occur. Either type of
mutation (deletion/amplification) will invoke the above described process for
the Genome object that the particular cell hosts. The tendency to mutate during
a spawning action is determined by a `mutation_rate` similar to the above
described `action_rates`

### Cell states
To mimic the "patchy" structure of cancer clones in spatial data, we actually
need to introduce some more dynamic behavior in our population. Thus, rather
than letting the cell behave the same at all times, we introduce the idea of
having two different _states_. One state (`state 0`) where the cells mainly move
and spawn, but with basically no mutation (low `mutation_rate`), this allows
certain clones to expand and populate the tissue domain. The second state
(`state 1`) the cells are less prone to be resting, but when they spawn they
have a much higher `mutation_rate`, this allows new clones to be introduced into
the population but also for the already existing ones to spread to a certain
extent. This dynamic behavior can easily be represented by a two state Markov model.

![mm](img/state.png)

The probabilities _p_<sub>ij</sub> describes how likely we are to go **from** state
_i_ **to** state _j_, these are referred to as `transition_probabilities`. The
state is inherited by a daughter cell from its parent cell, we also make sure
that the daughter and parent cell sync in their states, until the next change of
states. Meaning they will both change state at the same time, but then become
desynchronized, this is to prevent prolific introduction of mutations into the tissue.


### Lineage Tracing
Every cell will keep track of which cell that was it's parent, meaning that
lineage tracing is easily conducted. 


## Tissue Growth

Having defined how Cells and the Genome behaves, we can continue with to
actually construct our tissue _in silico_. The process is an iterative one,
which can be described as:

1. Place _N_<sub>0</sub> initial cells in the tissue domain. We will also
   immediately introduce an abberation into the genome of these cells, this can
   be either a deletion or amplification of desired size. By doing this, we are
   guaranteed to have _N_<sub>0</sub> distinctly different clonal populations in
   our tissue. These _N_<sub>0</sub> initial cells are added to our `population`
   of cells.

2. For every cell in the `population` call for an action to occur.
3. Repeat (2) until more than `max_cells` cells are present in the population or
   `max_iterations` iterations have been executed.
4. Fill in all of the empty space in the tissue domain with normal cells,
   sharing the same initial genome as the aberrant cells, but with no mutations.

This will generate a tissue akin to the one below, where each pixel represents
and individual cell and the color indicates the total number of aberrations in
the genome.

![ex_tissue](img/tissue_example.png)


We can also look at the genomic profile and what a sample of the expression
profile from these cells would look like. The images below represent these
features (genomic and expression profile) for the subset of aberrant cells.

![sc_profile](img/sc_profile.png)


### Clustering 

As we often are more interested in groups of cells than perhaps the individual
cell itself, this implementation also implements an agglomerative clustering
step - based on the genomic profile of each cell - which shows how these
clusters share a clear spatial zonation pattern:

![sc_cluster](img/sc_cluster.png)

Since we have all of the lineages traced, we can relate these clusters to each
other, and build a form of lineage tree (here represented as a graph). The root
nodes are colored in lightblue, and arrows go from parent to child.

![sc_cluster_lineage](img/lineage.png)



# Use
Synthethic data generation


# Requirements



