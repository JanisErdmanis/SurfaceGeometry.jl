function [p,t] = meshgen(a,b,c,step) %,boxsize)
    function f = fd(p)
        f = p(:,1).^2/a^2+p(:,2).^2/b^2+p(:,3).^2/c^2-1;
    end
    
    [p,t]=distmeshsurface(@fd,@huniform,step,[-2.1,-1.1,-1.6; 2.1,1.1,1.6]);
end

function h=huniform(p,varargin)
    h=ones(size(p,1),1);
end

