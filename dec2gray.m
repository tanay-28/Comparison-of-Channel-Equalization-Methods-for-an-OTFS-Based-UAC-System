function [grayCoded]=dec2gray(decInput)
%convert decimal to Gray code representation %example:x=[01234567]%decimal %y=dec2gray(x) %returnsy=[01326754]%Graycoded [rows,cols]=size(decInput); grayCoded=zeros(rows,cols);
[rows,cols]=size(decInput);
grayCoded=zeros(rows,cols);
for i=1:rows
 for j=1:cols
   grayCoded(i,j)=bitxor(bitshift(decInput(i,j),-1),decInput(i,j));
 end
end