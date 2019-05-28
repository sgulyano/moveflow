function Tx = ffd_interpolate(Ox, spline )
%FFD_INTERPOLATE interpolate position of pixel by FFD Ox with control point
%index i_array and spline parameter u_index_array and spline function Bu
Tx = zeros(size(spline.i_array));
for l = 1:4
    indexO = spline.i_array+l-1;
    val = spline.Bu(l+spline.u_index_array);
    Tx = Tx + val.*Ox(indexO);
end
end

