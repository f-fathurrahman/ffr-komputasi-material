close all; clear all;

draw_spheres('HeH.xyz', 'He', [1 0 0], 0.1, 40)
draw_spheres('HeH.xyz', 'H', [0.8 0.8 0.8], 0.1, 40)
light
draw_bonds('HeH.xyz','He', [1 0 0], 'H', [0.8 0.8 0.8], 0.1, 3.2, 6)
