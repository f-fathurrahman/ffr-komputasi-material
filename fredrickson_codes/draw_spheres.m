function draw_spheres(filename, atom1, color, sphere_radius, sphere_resolution)
    [atomname,x,y,z] = textread(filename,'%s %f %f %f');
    natoms = size(x);
    natoms = natoms(1);
    for j = 1:natoms
        if(strcmp(atom1,atomname(j))==1)
        [x_sphere,y_sphere,z_sphere] = sphere(sphere_resolution);
        [phi, theta_matlab, r] = cart2sph(x_sphere,y_sphere,z_sphere);
        theta = pi/2 - theta_matlab;
        r_new = sphere_radius;
        [x_new, y_new, z_new] = sph2cart(phi,theta_matlab,r_new);
        x_new = x_new + x(j);
        y_new = y_new + y(j);
        z_new = z_new + z(j);
        hold on
        surf(x_new,y_new,z_new,'FaceColor',color, 'EdgeColor','none');
    end
end
axis equal;