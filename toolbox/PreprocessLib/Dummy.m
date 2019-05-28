%% Centerline Calc

tic
disp(['Centerline calculation started']);
vec2dotU = Vx2.*fU;
vec2dotV = Vy2.*fV;
vec2dotW = Vz2.*fW;
vec3dotU = Vx3.*fU;
vec3dotV = Vy3.*fV;
vec3dotW = Vz3.*fW;
% clear fU fV fW Vx2 Vy2 Vz2 Vx3 Vy3 Vz3
vec2dotGVF = vec2dotU+vec2dotV+vec2dotW;
vec3dotGVF = vec3dotU+vec3dotV+vec3dotW;
clear vec2dotU vec2dotV vec2dotW vec3dotU vec3dotV vec3dotW;
centerlines = ((vec2dotGVF+vec3dotGVF).^2)<.0001;
stuff=((vec2dotGVF+vec3dotGVF).^2);
% clear vec2dotGVF vec3dotGVF;
% disp(['Vesselness from GVF calculation started']);
% [Iout,Vx,Vy,Vz]=NFrangiGVF3D(ux,uy,uz,vy,vz,wz,options,m,n,o);
toc