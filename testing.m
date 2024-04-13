% Example for plate of 4x8 feet made from 26 gauge galvanized steel

length = 1.2192;
width = 2.4384;
thickness = 0.00045;
n_modes = 100;
density = 7800;
youngs_mod = 210;
poisson = 0.29;


[fLong,fTrans,fBend,fEig] = calculatePlateModes(length ...
    ,width,thickness,n_modes,density,youngs_mod,poisson);


for i = 1:9
    disp(["Bending Mode length ", i ,": " ,fBend(1,i)]);
    disp(["Bending Mode width ", i ,": " ,fBend(2,i)]);
end