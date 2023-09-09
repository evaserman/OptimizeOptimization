close all
clear 
clc

nbglr       = 3;        % number of bays going from left to right
nbgtb       = 3;        % number of bays going top to bottom
cScale      = 4;       % number of elements in each bay
szmv        = 4;        % allowable movement distance of each node
cover       = 5;       % number of elements around border of problem
rmin        = 2.5;      % minimum feature size
Rc          = 4;        % radius of influence of truss nodes
eta         = 3;        % penalty parameter
vmax_c      = 0.4;      % maximum continuum material fraction
vmax_t      = 0.1;      % maximum truss material fraction
max_iter    = 50;       % maximum iterations before termination
loadcase    = 1;        % 1 - Cantilever    2 - MBB
vsm         = 0;        % Variable Thickness Sheet Switch
nelx        = nbglr*cScale+cover;   %number of elements in x
nely        = nbgtb*cScale+cover*2; %number of elements in y


time_top88(nbglr,nbgtb,nelx,nely,rmin,eta,vmax_c,vmax_t,max_iter,...
    loadcase,vsm,szmv,Rc,cScale,cover)