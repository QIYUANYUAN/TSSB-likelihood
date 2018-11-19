treeF = None
freqF = None
n_t = 297
tree_in = "T_297.txt"
freq_in = "F_297.txt"
result = "result_297.csv"


def openfile(tree_file, freq_file):
    global treeF,freqF;
    treeF=open(tree_file,"r");
    freqF=open(freq_file,"r");

import numpy as np
import itertools
import math

def bpdf(w,b):
    if (1-w < 1e-9): w=1-1e-9;
    return b* ( (1-w)**(b-1) )

class tree(object):

    alpha_decay=0.25
    alpha=25.0
    gamma=1.0
    nodes = []
    indices = {}
    parents = {}
    #ancestors = {}
    children = {}
    B = []
    U = []  #stick len
    F = []
    #w = []  #stick len subtree
    
    def __init__(self, nodes=[], f=[]):
        self.F=f;
        self.nodes=[v for v in nodes];
        self.indices = {v:i for i,v in enumerate(nodes)}
        self.parents = {v:None for v in nodes};
        self.children = {v:[] for v in nodes};
        #self.ancestors = {v:[] for v in nodes};
        self.B=[[] for v in nodes];

    def add_arc(self,parent,child):
        self.children[parent].append(child);
        self.parents[child]=parent;

    def creat_B(self,node):
        i=self.indices[node]
        if (self.parents[node]==None):
            self.B[i]=list(np.zeros(len(self.nodes)));
            self.B[i][i]=1.0;
        else:
            self.B[i]=list(self.B[self.indices[self.parents[node]]]);
            self.B[i][i]=1.0;
        #print(i,self.B[i],sep="\n");
        for v in self.children[node]:
            self.creat_B(v);

    def creat_phi(self):
        #print(B)
        invB=np.linalg.inv(self.B);
        self.U=np.matmul(F,invB);
        self.w=np.zeros(len(self.nodes))


    #def creat_stick_len_subtree(self, node):
    #    for v in self.children[node]:
    #        self.creat_stick_len_subtree(v);
    #        self.w[self.indices[node]]+=self.w[self.indices[v]];
    #    self.w[self.indices[node]]+=self.U[self.indices[node]];
    
    def cal_log_likelihood(self,node,depth):
        #print("--------------------", self.indices[node], "-------------------")
        if (self.F[self.indices[node]]<=1e-9): return 0;
        weight=self.U[self.indices[node]]/self.F[self.indices[node]];
        likelihood_self=bpdf(weight,(self.alpha_decay**depth)*self.alpha);
        #print(weight,likelihood_self)
        log_like=math.log(likelihood_self);
        likelihood_child=0.0;
        for p in itertools.permutations(self.children[node]):
            stick=self.F[self.indices[node]]-self.U[self.indices[node]];
            p_like=1.0;
            #print("---",p,stick,"---")
            for v in p:
                if(self.F[self.indices[v]]>1e-9): p_like*=bpdf(self.F[self.indices[v]]/stick,self.gamma);
                #print(self.F[self.indices[v]]/stick,self.gamma)
                stick-=self.U[self.indices[v]];
            likelihood_child+=p_like;
        log_like+=math.log(likelihood_child);
        for v in self.children[node]:
            log_like+=self.cal_log_likelihood(v,depth+1);
        #print("log_like", log_like);
        return log_like;


    
        
        
openfile("Data/"+tree_in,"Data/"+freq_in);


nodes=[]
F=[]

for line in freqF:
    nodes.append(line.split()[0]);
    F.append(float(line.split()[1]));

ntrees=int(treeF.readline().split()[0]);

result=open("Result/"+result,"w");

def get_next_tree():
#if True:
    nedges=int(treeF.readline().split()[0]);
    Tree=tree(nodes,F)
    for i in range(6):
        parent,child=treeF.readline().split();
        #print(parent,child);
        Tree.add_arc(parent, child);
    Tree.creat_B(nodes[0]);
    Tree.creat_phi();
    ll=Tree.cal_log_likelihood(nodes[0],0);
    
    result.write("%f\n"%ll);

for _ in range(n_t):
    get_next_tree();
    pass


result.close();
