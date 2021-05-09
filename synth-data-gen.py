
# coding: utf-8

# In[2]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams

import copy


# In[3]:


rcParams["figure.facecolor"] = "white"


# In[4]:


class Counter:
   def __init__(self):
       self.cnt = 0 
   def __call__(self):
       self.cnt += 1
       return self.cnt 


# In[5]:


from scipy.stats import truncnorm

class Genome:
    def __init__(self,
                 theta,
                 genome_length: int,
                 expr_ub = 1000,
                 expr_lb = 100,
                 sigma: int = 2,
                ):
        
        self.theta = theta
        self.G = theta.size
        
        self.gene_id = np.arange(self.G)
        self.pos = np.zeros(self.G)
        
        
        self.genome_length = genome_length
        self._build()
        
        self.sigma_raw = sigma
        self.sigma = self.sigma_raw * np.diff(self.pos).mean()
        self.len_dist = truncnorm(0,np.inf,loc=np.diff(self.pos).mean(), scale = self.sigma)
        
        self.expr_ub = expr_ub
        self.expr_lb = expr_lb
        
    
        
    def _build(self,):
        pos = np.random.random(self.G)
        pos[0] = 0.0
        pos = np.cumsum(pos)
        pos = pos / pos[-1]
        self.pos = pos *  self.genome_length
 
    def in_range(self,
                center: float,
                 width: float,
                ):
        return np.where(np.abs((self.pos - center)) < width)[0]
    
    def add(self,
           idx: np.ndarray,
           ):
        
        
        if len(idx) > 0:
            start_pos = idx[0]
            end_pos = idx[-1]


            amp_lens = self.pos[start_pos:end_pos+1]
            
            delta = amp_lens[-1] - amp_lens[0]
            p1_pos = self.pos[0:end_pos+1]
            p2_pos = amp_lens + delta
            p3_pos = self.pos[end_pos+1::] + delta
            pos = np.hstack((p1_pos,p2_pos,p3_pos))

            p1_gen = self.gene_id[0:end_pos+1]
            p2_gen = self.gene_id[idx]
            p3_gen = self.gene_id[end_pos+1::]
            gen_id = np.hstack((p1_gen,p2_gen,p3_gen))

            self.pos = pos
            self.gene_id = gen_id
            self.genome_length = self.pos[-1]
        
    def delete(self,
              idx: np.ndarray,
              ):
        
          
        if len(idx) > 0:
            start_pos = idx[0]
            end_pos = idx[-1]


            amp_lens = self.pos[start_pos:end_pos+1]
            delta = amp_lens[-1] - amp_lens[0]
            p1_pos = self.pos[0:start_pos]
            p2_pos = self.pos[end_pos+1::] - delta
            self.pos = np.hstack((p1_pos,p2_pos))
            
            self.gene_id = np.hstack((self.gene_id[0:start_pos],self.gene_id[end_pos+1::]))
            self.genome_length = self.pos[-1]
        
        
    def sample_expression(self,):
        scale_factor = self.genome_length / self.G
        lb = self.expr_lb * scale_factor
        ub = self.expr_ub * scale_factor
        n_genes = np.random.randint(lb,ub)
        pvals = self.theta[self.gene_id]
        pvals /= pvals.sum()
        x_raw = np.random.multinomial(n_genes,pvals)
        x_cur = np.zeros(self.G)
        for k,g in enumerate(self.gene_id):
            x_cur[g] += x_raw[k]
        return x_cur
        
    def sample_event(self,)->[float,float]:
        center = np.random.uniform(0,self.genome_length)
        width = self.len_dist.rvs()
        return (center,width)
    
    def get_gene_count(self,)->np.ndarray:
        counts = np.zeros(self.G)
        for g in self.gene_id:
            counts[g] += 1
        return counts
        
        
        
        
        


# In[157]:


class Cell:
    def __init__(self,
                 genome,
                 mutation_rate,
                 action_rates,
                 counter,
                 x = 0,
                 y = 0,
                 active_state = 5,
                 p_move = None,
                 tmat = None,
                 parent = -1
                ):
          
        self.counter = counter
        self.id = self.counter()
        self.parent = parent
        self.children = []
            
        self.genome = genome
        
        self.N = self.genome.G
         
        self.mutation_rate = mutation_rate
        self.action_rates = action_rates
        
        self.action_rates /= self.action_rates.sum(axis=0,keepdims=True)
        
        if tmat is None:
            self.tmat = np.array([[0.8,0.2],[0.4,0.6]]).T
        else:
            self.tmat = tmat
        
        self.active = active_state
        
        assert np.all(self.action_rates >= 0),        "rates must be non-negative than zero"
        self.crd = np.array([x,y])
        
        
        if p_move is None:
            self.p_move = np.random.dirichlet(0.1*np.ones(4))
        else:
            self.p_move = p_move
            
        self.gamma = 0.05
                
        self.dx = [-1,0,1,0]
        self.dy = [0,-1,0,1]
        self.step_scale = [1,2]
        
    
    def search_nbrhd(self,
                     cell_map,
                    ):
        
        free_space = []
        pvals = []
        scale_x = np.random.choice(self.step_scale)
        scale_y = np.random.choice(self.step_scale)
        for k,(dx,dy) in enumerate(zip(self.dx,self.dy)):
            new_x = self.crd[0] + dx * scale_x
            new_y = self.crd[1] + dy * scale_y
            
            if (new_x < 0) | (new_x >= cell_map.shape[0]):
                break
            elif (new_y < 0) | (new_y >= cell_map.shape[1]):
                break
            
            if cell_map[new_x,new_y] == 0:
                free_space.append((new_x,new_y))
                pvals.append(self.p_move[k])
                
        pvals = [x/sum(pvals) for x in pvals]
        return free_space,pvals
        
    def genome_update(self,):
        # check if mutation should occur
        if np.random.random() < self.mutation_rate[self.state]:
            self.mutation_rate[self.state] *= 0.95
            c,l = self.genome.sample_event()
            idx = self.genome.in_range(c,l)
            
            if np.random.random() < 0.5:
                self.genome.add(idx)
            else:                
                self.genome.delete(idx)
                
    def sample_expression(self,)->np.ndarray:
        return self.genome.sample_expression()
    
    @property
    def state(self,):
        return int(self.active < 0)
            
    def action(self,
               cell_map):
        
        if self.active < 0:
            new_state = np.random.multinomial(1,self.tmat[:,1]).argmax()
            if new_state == 0:
                self.active = np.random.geometric(p = self.tmat[1,0])
        else:
            self.active -= 1
             
        action = np.random.multinomial(1,self.action_rates[:,self.state]).argmax()
        if action == 3 and cell_map.sum() > 50:
            cell_map[self.crd[0],self.crd[1]] = 0
            return []
        else:
            free_space,pvals = self.search_nbrhd(cell_map)
            if len(free_space) > 0:
                new_crd_idx = np.random.choice(len(free_space),p=pvals)
                new_crd = np.array(free_space[new_crd_idx])
                if action == 1:
                    self.genome_update()
                    new_genome = copy.deepcopy(self.genome)
                    new_cell = Cell(genome = new_genome,
                                    mutation_rate = self.mutation_rate,
                                    action_rates = self.action_rates,
                                    counter = self.counter,
                                    x = new_crd[0],
                                    y = new_crd[1],
                                    active_state =self.active,
                                    p_move = (self.p_move if self.active >= 0 else None),
                                    tmat = self.tmat,
                                    parent = self.id,
                                   )
                    self.children.append(new_cell.id)
                    
                    cell_map[new_crd[0],new_crd[1]] = 1
                    
                    return [self,new_cell]
                
                elif action == 2:
                    cell_map[self.crd[0],self.crd[1]] = 0
                    self.crd = new_crd
                    cell_map[self.crd[0],self.crd[1]] = 1
                    return [self]
                else:
                    return [self]
            else:
                    return [self]
            
            


# In[158]:


def get_coordinstes(pop):
    crd = np.zeros((len(pop),2))
    for k,cell in enumerate(pop):
        crd[k,0] = cell.crd[0]
        crd[k,1] = cell.crd[1]
    
    return crd


# In[159]:


tmat = np.array([[0.8,0.2],[0.2,0.8]]).T
mutation_rate = np.array([0.0,0.8])
action_rates = np.array(((0.1,0.3),
                         (0.4,0.45),
                         (0.4,0.2),
                         (0.0,0.00)
                        ))
action_rates = action_rates / action_rates.sum(axis=0)
action_rates


# In[160]:


cell_map = np.zeros((200,200))
N = 500
theta = np.random.dirichlet(np.ones(N))
genome = Genome(theta,1e5,sigma = 10)
counter = Counter()
max_cells = 1e4
lineage = {}
pl = {}

cell_0 = Cell(genome = copy.deepcopy(genome),
              counter = counter,
              mutation_rate = mutation_rate,
              action_rates = action_rates,
              tmat = tmat,
              x = 50,y = 50)
cell_1 = Cell(genome = copy.deepcopy(genome),
              counter = counter,
              mutation_rate = mutation_rate,
              action_rates = action_rates,
              tmat = tmat,
              x = 150,y = 150)



lineage.update({cell_0.parent:[cell_0.id]})
lineage.update({cell_1.parent:[cell_1.id]})


cell_0.genome.add([x for x in range(3,20)])
cell_1.genome.delete([x for x in range (250,265)])

population = [cell_0,cell_1]


for t in range(500):
    new_cells = [] 
    for k in range(len(population)):
        cell = population.pop(0)
        out = cell.action(cell_map)
        for c in out:
            if c.id not in lineage:
                lineage.update({c.id:c.parent})
                
        new_cells += out
        
        
    population += new_cells
    if cell_map.sum() > max_cells:
        break


# In[161]:


crds = np.array([x.crd for x in population])
abber = np.array([(x.genome.get_gene_count() != 1).sum() for x in population])


# In[162]:


plt.figure(figsize = (10,10))
plt.scatter(crds[:,0],crds[:,1],c = abber, s = 2)
plt.axis("equal")
plt.colorbar()
plt.show()


# In[164]:


abber = np.vstack([x.genome.get_gene_count()[np.newaxis,:] for x in population])
plt.figure(figsize = (10,7))
plt.imshow(abber,cmap = plt.cm.viridis)
plt.axis("auto")
plt.colorbar()
plt.xlabel("Gene",fontsize = 20)
plt.ylabel("Cell",fontsize = 20)
plt.show()


# In[166]:


expression = np.array([x.sample_expression() for x in population])
genome_profile = np.vstack([x.genome.get_gene_count()[np.newaxis,:] for x in population])


# In[167]:


plt.figure(figsize = (10,7))
plt.imshow(expression,cmap = plt.cm.viridis)
plt.axis("auto")
plt.colorbar()
plt.xlabel("Gene",fontsize = 20)
plt.ylabel("Cell",fontsize = 20)
plt.show()


# In[171]:


from sklearn.cluster import AgglomerativeClustering
n_cluster = 20
cluster = AgglomerativeClustering(n_clusters=n_cluster, affinity='euclidean', linkage='ward')
idx = cluster.fit_predict(genome_profile)
plt.figure(figsize = (10,10))
plt.scatter(crds[:,0],crds[:,1],c = idx, s = 5,cmap = plt.cm.tab20)
plt.show()


# In[196]:


cluster_genome = np.zeros((n_cluster,N))
for clu in np.unique(idx):
    pos = (idx == clu)
    cluster_genome[clu,:] = genome_profile[pos,:].mean(axis=0).round()
    
plt.imshow(cluster_genome)
plt.axis("auto")


# In[172]:


pos_to_id = {x:population[x].id for x in range(len(population))}
id_to_pos = {v:k for k,v in pos_to_id.items()}
cmat = np.zeros((n_cluster,n_cluster))

for k1 in range(n_cluster):
    c_members = np.where(idx == k1)[0]
    for m in c_members:
        i = pos_to_id[m]
        k2 = k1
        while (k2 == k1) and (i > 0) :
            i = lineage[i]
            if i in id_to_pos.keys():
                k2 = idx[id_to_pos[i]]
                
        cmat[k1,k2] +=1
        
plt.imshow(cmat)


# In[84]:


normal_population = []
row,col = np.where(cell_map == 0)
for x,y in zip(row,col):
    normal_population.append(Cell(x=x,y=y,genome=copy.deepcopy(genome),
                                  counter = counter,
                                  mutation_rate = mutation_rate,action_rates = action_rates,tmat=tmat))


# In[85]:


normal_expression = np.array([x.sample_expression() for x in normal_population])
normal_crd = np.array([x.crd for x in normal_population])


# In[175]:





# In[185]:


all_population = population + normal_population
all_cluster = np.append(idx,(max(idx)+1) * np.ones(len(normal_population)))
all_expression = np.vstack((expression,normal_expression))
all_crd = np.vstack((crds,normal_crd))
all_abber = np.array([(x.genome.get_gene_count()[np.newaxis,:] != 1).sum() for x in all_population])
all_genome = np.vstack([x.genome.get_gene_count()[np.newaxis,:] for x in all_population])
benign = np.array([int(np.all(x.genome.get_gene_count() ==1)) for x in all_population])


# In[96]:


xx = np.arange(5,cell_map.shape[0] + 5,10)
yy = np.arange(5,cell_map.shape[1] + 5,10)
xx,yy = np.meshgrid(xx,yy)
xx = xx.flatten()
yy = yy.flatten()


# In[97]:


plt.figure(figsize = (10,10))
plt.scatter(all_crd[:,0],all_crd[:,1],c = all_abber, s = 5,cmap = plt.cm.magma)
plt.scatter(xx,yy,alpha = 0.5,facecolor = "red",edgecolor = "red", s = 100,)
plt.axis("equal")
plt.show()


# In[187]:


exp_matrix = np.zeros((xx.shape[0],N))
gen_matrix = np.zeros((xx.shape[0],N))
ben_matrix = np.zeros(xx.shape[0])
cluster_matrix = np.zeros((xx.shape[0],n_cluster+1))
depth = 5000
radius = 2.5
radius_squared = 2.5**2
for spot,(x,y) in enumerate(zip(xx,yy)):
    in_spot = []
    for k in range(all_crd.shape[0]):
        delta = (all_crd[k,0] - x)**2 + (all_crd[k,1]-y)**2
        if delta < radius_squared:
            in_spot.append(k)
    in_spot = np.array(in_spot)
    joint_expr = all_expression[in_spot,:].sum(axis=0)
    av_gen = all_genome[in_spot].mean(axis=0)
    state_ben = int(benign[in_spot].mean() >= 0.9)
    gene_prob = joint_expr / joint_expr.sum()
    ns = np.min((depth,joint_expr.sum()))
    exp_matrix[spot,:] = np.random.multinomial(ns,gene_prob)
    gen_matrix[spot,:] = av_gen
    ben_matrix[spot] = state_ben
    clu,cnt = np.unique(all_cluster[in_spot],return_counts = True)
    cluster_matrix[spot,clu.astype(int)] = cnt / cnt.sum()
    


# In[188]:


plt.imshow(cluster_matrix)
plt.axis("auto")


# In[189]:


fig,ax = plt.subplots(1,2,figsize = (10,5))
ax[0].scatter(xx,yy,c = exp_matrix[:,10],s = 120)
ax[0].set_title("Gene 13",fontsize = 20)
ax[1].scatter(xx,yy,c = exp_matrix[:,255],s = 120)
ax[1].set_title("Gene 255",fontsize = 20)
plt.show()


# In[190]:


plt.imshow(gen_matrix[ben_matrix==0,:])
plt.axis("auto")
plt.show()


# In[528]:


n_spots = exp_matrix.shape[0]
n_genes = exp_matrix.shape[1]
index = ["Spot_{}".format(x) for x in range(n_spots)]
columns = ["Gene_{}".format(x) for x in range(n_genes)]
expression_data = pd.DataFrame(exp_matrix[:,0:-1],
                               index = index,
                               columns = columns[0:-1],
                              )

expression_data = expression_data.T
annotation_data = pd.DataFrame(ben_matrix,
                               index = index,
                               columns = ["annotation"],
                              )

annotation_data.annotation = annotation_data.annotation.map({1:"benign",0:"aberrant"})

start = genome.pos[0:-1]
end = genome.pos[1::]
chrom = np.array(["chr1" for x in range(start.shape[0])])
tmp = np.hstack((chrom[:,np.newaxis],start[:,np.newaxis],end[:,np.newaxis])) 
gene_data = pd.DataFrame(tmp,
                         index = columns[0:-1],
                         columns = ["chr","start","end"],
                        )


# In[529]:


OUT_DIR = "/home/alma/w-projects/andrew-prostate/res/synth"


# In[530]:


import os.path as osp
expression_data.to_csv(osp.join(OUT_DIR,"expression.tsv"),
                      sep ="\t",
                       header = True,
                       index= True,
                      )

annotation_data.to_csv(osp.join(OUT_DIR,"annotation.tsv"),
                      sep = "\t",
                       header = False,
                       index = True,
                      )

gene_data.to_csv(osp.join(OUT_DIR,"gene_data.tsv"),
                 sep = "\t",
                 header = False,
                 index = True,
                )

