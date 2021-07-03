function draw_bonds(filename,atom1,color1,atom2,...
    color2,d_min,d_max,width)
    [atomname,x,y,z] = textread(filename,'%s %f %f %f');
    natoms = size(x);
    natoms = natoms(1);
    for j = 1:natoms
        for k = 1:natoms
            if( strcmp(atom1,atomname(j) ) == 1 )
                if( strcmp(atom2,atomname(k)) == 1 )
                    dist1 = ( (x(j)-x(k))^2 + (y(j)-y(k))^2 + (z(j)-z(k))^2 )^0.5;
                    fprintf('dist = %f\n', dist1)
                    if( (dist1<=d_max) && (dist1>=d_min) )
                        fprintf('Should draw bonds ...\n')
                        % TIME TO DRAW A BOND
                        hold on;
                        midpnt = [(x(j)+x(k))/2 (y(j)+y(k))/2 (z(j)+z(k))/2];
                        % Draw first half of bond with color1
                        plot3([x(j) midpnt(1)],[y(j) midpnt(2)], ...
                            [z(j) midpnt(3)], 'color',color1, ...
                            'linewidth',width);
                        % Draw second half of bond with color2
                        plot3([midpnt(1) x(k)],[midpnt(2) y(k)], ...
                              [midpnt(3) z(k)], 'color', color2, ...
                              'linewidth',width);
                    end
                end
            end
        end
    end
end
