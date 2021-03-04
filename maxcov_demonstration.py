# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 14:59:52 2021

@author: Vivian Imbriotis
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation, PillowWriter
import seaborn as sns

sns.set_style("dark")


import numpy as np


np.seterr(all="ignore")

def cov(a,b):
    return ((a-a.mean()) * (b-b.mean())).sum()

class Cube:
    color = "green"
    surf_alpha = 0.3
    edge_alpha = 0.8
    def __init__(self,ax,xmin,xmax,ymin,ymax,zmin,zmax):
        r = [-1,1]
        cube = np.zeros((2,2,2,3))
        
        #Set x-coords
        cube[0,:,:,0] = xmin
        cube[1,:,:,0] = xmax
        
        #Set y-coords
        cube[:,0,:,1] = ymin
        cube[:,1,:,1] = ymax
        
        #Set z-coords
        cube[:,:,0,2] = zmin
        cube[:,:,1,2] = zmax
        self.coords = cube
        self.ax = ax
        self.plot_everything()
        
    def plot_everything(self):
        cube = self.coords
        self.artists = [
            self.ax.plot_surface(cube[0,:,:,0],cube[0,:,:,1],cube[0,:,:,2],alpha=self.surf_alpha,color=self.color),
            self.ax.plot_surface(cube[1,:,:,0],cube[1,:,:,1],cube[1,:,:,2],alpha=self.surf_alpha,color=self.color),
            
            self.ax.plot_surface(cube[:,0,:,0],cube[:,0,:,1],cube[:,0,:,2],alpha=self.surf_alpha,color=self.color),
            self.ax.plot_surface(cube[:,1,:,0],cube[:,1,:,1],cube[:,1,:,2],alpha=self.surf_alpha,color=self.color),
            
            self.ax.plot_surface(cube[:,:,0,0],cube[:,:,0,1],cube[:,:,0,2],alpha=self.surf_alpha,color=self.color),
            self.ax.plot_surface(cube[:,:,1,0],cube[:,:,1,1],cube[:,:,1,2],alpha=self.surf_alpha,color=self.color),
        
            self.ax.plot_wireframe(cube[0,:,:,0],cube[0,:,:,1],cube[0,:,:,2],alpha=self.edge_alpha,color=self.color),
            self.ax.plot_wireframe(cube[1,:,:,0],cube[1,:,:,1],cube[1,:,:,2],alpha=self.edge_alpha,color=self.color),
            
            self.ax.plot_wireframe(cube[:,0,:,0],cube[:,0,:,1],cube[:,0,:,2],alpha=self.edge_alpha,color=self.color),
            self.ax.plot_wireframe(cube[:,1,:,0],cube[:,1,:,1],cube[:,1,:,2],alpha=self.edge_alpha,color=self.color),
            
            self.ax.plot_wireframe(cube[:,:,0,0],cube[:,:,0,1],cube[:,:,0,2],alpha=self.edge_alpha,color=self.color),
            self.ax.plot_wireframe(cube[:,:,1,0],cube[:,:,1,1],cube[:,:,1,2],alpha=self.edge_alpha,color=self.color),
        ]
        
    def remove_everything(self):
        for artist in self.artists:
            artist.remove()
        self.artists = []
    def move(self,dx,dy,dz):
        self.remove_everything()
        self.coords += np.array((dx,dy,dz))
        #Update surfaces
        self.plot_everything()
        
    def contains(self,x,y,z):
        #check x dimention
        # print(f"X between {(self.coords[0,0,0,0],self.coords[1,0,0,0])}")
        # print(f"Y between {(self.coords[0,0,0,1],self.coords[0,1,0,1])}")
        # print(f"Z between {(self.coords[0,0,0,2],self.coords[0,0,1,2])}")
        return all((self.coords[0,0,0,0]<x<self.coords[1,0,0,0],
                    self.coords[0,0,0,1]<y<self.coords[0,1,0,1],
                    self.coords[0,0,0,2]<z<self.coords[0,0,1,2]))
        
    def __del__(self):
        self.remove_everything()
    

class MaxCovFigure:
    deg_per_frame=0.3
    n_genpop=1500
    n_taxon=200
    motion_per_frame = 0.1
    frames = int(3*(4+4)/motion_per_frame) - 2
    def __init__(self,taxon=True,cube=True,scatter=True,cov=True):
        np.random.seed(0)
        self.draw_taxon = taxon
        self.draw_cube = cube
        self.draw_scatter = scatter
        self.draw_covline = cov
        self.fig = plt.figure(figsize = (16,9))
        
        #2D Axes
        self.ax_xy = self.fig.add_subplot(2,2,1)
        if self.draw_scatter: self.ax_xy.set_title("2D Slice of Data")
        self.ax_xy.set_xlim((-4,4))
        self.ax_xy.set_ylim((-4,4))
        self.xy_genpop_artist, = self.ax_xy.plot([],[],'o',color="blue")
        self.xy_taxon_artist, = self.ax_xy.plot([],[],'o',color="red")
        
        self.ax_cov = self.fig.add_subplot(2,2,3)
        if self.draw_covline: self.ax_cov.set_title("Covariance within slice")
        self.ax_cov.set_xlim((-4,4))
        self.ax_cov.set_ylim((-200,200))
        self.cov_artist, = self.ax_cov.plot([],[],color="k")
        self.ax3d = self.fig.add_subplot(1,2,2,projection="3d")
        self.ax3d.set_xlim((-4,4))
        self.ax3d.set_ylim((-4,4))
        self.ax3d.set_zlim((-4,4))
        self.gen_pop = np.random.normal(size = (self.n_genpop,3))
        self.taxon   = np.random.normal(size = (self.n_taxon,3)) * 0.4 + 0.8*np.array((2,2,-3))
        self.ax3d.scatter(self.gen_pop[:,0],self.gen_pop[:,1],self.gen_pop[:,2],color="blue")
        if self.draw_taxon: self.ax3d.scatter(self.taxon[:,0],self.taxon[:,1],self.taxon[:,2],color="red")
        self.ax3d.set_xlabel("X")
        self.ax3d.set_ylabel("Y")
        if self.draw_cube: self.cube = Cube(self.ax3d,-4,-3.5,-4,4,-4,4)
        self.ani = FuncAnimation(self.fig, self.update, 
                                 interval=50,frames=self.frames)
        self.state = 0
        self.idx=0
        self.state_change = False
        self.covline = np.ones((int(8/self.motion_per_frame)+2,2)) * np.nan
        if self.draw_taxon: self.all = np.concatenate((self.gen_pop,self.taxon))
        else: self.all = self.gen_pop
        
    def rotate(self,angle):
         self.ax3d.view_init(elev = 15, azim=angle)
    def update_cube(self):
         move = [0,0,0]
         move[self.state] += 0.1
         self.cube.move(*move)
         if self.cube.coords[0,0,0,self.state] > 4:
             self.cube.remove_everything()
             self.state += 1
             self.state_change = True
             self.state = self.state % 3
             if self.state == 0:
                 self.cube = Cube(self.ax3d,-4,-3.5,-4,4,-4,4)
             elif self.state == 1:
                 self.cube = Cube(self.ax3d,-4,4,-4,-3.5,-4,4)
             elif self.state == 2:
                 self.cube = Cube(self.ax3d,-4,4,-4,4,-4,-3.5)
    def update(self,frame):
        self.rotate(angle=frame*self.deg_per_frame)
        if self.draw_cube:    self.update_cube()
        if self.draw_scatter: self.draw_xy()
        if self.draw_covline: self.draw_cov()
    
    def draw_xy(self):
        contained = [p for p in self.gen_pop if self.cube.contains(*p)]
        containedt = [p for p in self.taxon if self.cube.contains(*p)]        
        if self.state==0:
            self.xy_genpop_artist.set_data([c[1] for c in contained],
                                           [c[2] for c in contained]) 
            if self.draw_taxon: self.xy_taxon_artist.set_data([c[1] for c in containedt],
                                           [c[2] for c in containedt]) 
            self.ax_xy.set_xlabel("Y-axis")
            self.ax_xy.set_ylabel("Z-axis")
        if self.state==1:
            self.xy_genpop_artist.set_data([c[0] for c in contained],
                                           [c[2] for c in contained]) 
            if self.draw_taxon: self.xy_taxon_artist.set_data([c[0] for c in containedt],
                                           [c[2] for c in containedt]) 
            self.ax_xy.set_xlabel("X-axis")
            self.ax_xy.set_ylabel("Z-axis")
        if self.state==2:
            self.xy_genpop_artist.set_data([c[0] for c in contained],
                                           [c[1] for c in contained]) 
            if self.draw_taxon: self.xy_taxon_artist.set_data([c[0] for c in containedt],
                                           [c[1] for c in containedt]) 
            self.ax_xy.set_xlabel("X-axis")
            self.ax_xy.set_ylabel("Y-axis")

    def draw_cov(self):
        if self.state_change:
            self.state_change=False
            self.covline = np.ones((int(8/self.motion_per_frame)+2,2)) * np.nan
            self.idx = 0
        contained = np.array([p for p in self.all if self.cube.contains(*p)])
        x = self.cube.coords[0,0,0,self.state]
        if contained.size==0:
            y = 0
        elif self.state==0:
            y = cov(contained[:,1],contained[:,2])
        elif self.state==1:
            y = cov(contained[:,0],contained[:,2])
        elif self.state==2:
            y = cov(contained[:,0],contained[:,1])
        self.covline[self.idx] = (x,y)
        self.cov_artist.set_data(self.covline[:,0],self.covline[:,1])
        self.idx += 1
    def save(self,filename = "result.gif"):
        self.ani.save(filename = filename,writer = PillowWriter(fps=20))
        
            
            
     
if __name__=="__main__":
    plt.ioff()
    fig = MaxCovFigure(taxon=False,cube=False,scatter=False,cov=False).save(
        "dim_nocube_nocov.gif")
    fig = MaxCovFigure(taxon=False,cube=True,scatter=True,cov=False).save(
        "dim_cube_nocov.gif")
    fig = MaxCovFigure(taxon=False,cube=True,scatter=True,cov=True).save(
        "dim_cube_cov.gif")
    fig = MaxCovFigure(taxon=True,cube=False,scatter=False,cov=False).save(
        "tax_nocube_nocov.gif")
    fig = MaxCovFigure(taxon=True,cube=True,scatter=True,cov=False).save(
        "tax_cube_nocov.gif")
    fig = MaxCovFigure(taxon=True,cube=True,scatter=True,cov=True).save(
        "tax_cube_cov.gif")