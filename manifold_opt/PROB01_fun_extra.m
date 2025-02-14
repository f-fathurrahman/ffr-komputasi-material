function data = PROB01_fun_extra(data)
    XP = data.XP;
    data.rhoX = sum(XP.^2, 2);
end