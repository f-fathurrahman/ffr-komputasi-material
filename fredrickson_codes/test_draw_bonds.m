clear all; close all
width = 3;
draw_bonds('HO.xyz', 'H', 'red', 'O', 'blue', 0.5, 2.5, width)
draw_spheres('HO.xyz', 'H', 'red', 0.1, 10)
draw_spheres('HO.xyz', 'O', 'blue', 0.1, 10)
