function vTmp = nonlin(v)
    vTmp = v.*(v - 0.1).*(1 - v);

% vTmp = [];
%     for j=1:size(v, 2)
%         for i=1:size(v, 1)
%             vTmp(i, j) = v(i, j)*(v(i, j) - 0.1)*(1 - v(i, j));
%         end
%     end
end
