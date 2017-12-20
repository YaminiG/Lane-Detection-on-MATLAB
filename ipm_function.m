function IPM = ipm_function(image)

global fu fv cu cv h c1 c2 s1 s2

T = h*[-1/fu*c2 1/fv*s1*s2 (1/fu*cu*c2)-(1/fv*cv*s1*s2)-(c1*s2) 0; 1/fu*s2 1/fv*s1*c2...
       -(1/fu*cu*s2)-(1/fv*cv*s1*c2)-(c1*c2) 0; 0 1/fv*c1 -1/fv*cv*c1+s1 0; 0 -1/(h*fv)*c1 1/(h*fv)*cv*c1-s1/h 0];

image_hor = image(339:end,:);
[m,n] = size(image_hor);
proj = zeros(m,n,2);
for i=1:382
    for j = 1:1280
        comp_proj = T*[j;i;1;1];
        proj(i,j,1) = comp_proj(1);
        proj(i,j,2) = comp_proj(2);
    end
end    
end
