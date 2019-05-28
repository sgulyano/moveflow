function [ err, O_grad ] = bspl_ffd(O_trans, Bu, u_flow, Spacing, Stepsize )
%BSPL_FFD Summary of this function goes here
%   Detailed explanation goes here

err = 0;

Osize = size(O_trans);
Isize = size(u_flow);

O_grad = zeros(Osize);

% for every block of size 'Spacing'
i = 1;
for y = 1:Spacing(1):Isize(1)
    j = 1;
    for x = 1:Spacing(2):Isize(2)
        u_vec = u_flow(y:y+Spacing(1)-1, x:x+Spacing(2)-1);
        
        O_current = O_trans(i:i+3, j:j+3);
        u_init = Bu*O_current*Bu';
        
        du_init = u_init - u_vec;
        err_init = sum(du_init(:).^2);
        err = err + err_init;
        % for every 4x4 control points that affect the block
        for k = 1:4
            for l = 1:4
                % find changes in error when perturb by a step
                O_step = zeros(4,4);
                O_step(k,l) = Stepsize;

                u_step = Bu*(O_current+O_step)*Bu';
                du_step = u_step - u_vec;

                err_step = sum(du_step(:).^2);
                O_grad(i+k-1,j+l-1) = O_grad(i+k-1,j+l-1) + (err_step - err_init);
            end
        end
        
        j = j+1;
    end
    i = i+1;
end

end

