# JuliaMotMod
Simulation of light activated microtubule motor system

Example input Standard MotMod: 
MotorModelv03(["Filbounce.mat",50000,0,0,4000,200.0f0,0.0f0,-200.0f0,10000.0f0,100.0f0,100000.0f0,0.3f0,4.0f0,1.0f0,1.0f0,true])

Example input Torsion Spring MotMod:
i = 1;
j = 10;
MotorModelAngv01([string("LinePhase01_BeamH",i,"um_FilLen",j,"um.mat"),5000,0,0,300,200.0f0,0.0f0,-800.0f0,j*1000.0f0,100.0f0,j*1000.0f0,0.3f0,5.0f0,1.0f0,Float32(i/2)*1000.0f0,true])


Example input:
MotorModelv01("TwoAsterFormation4.mat",4000,0,0,400,800.0,0.0,-800.0,10000.0,100.0,13000.0,0.5,5.0,true);

Linking test: MotorModelv01("TwoAsterlink03_p1.mat",10000,0,0,200,800.0,0.0,-800.0,10000.0,100.0,10000.0,0.5,3.0,true);
