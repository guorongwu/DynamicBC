function [OutMatrix,OutMod,roire,ZeroInd,ReGroupInd,Outsortinfo] = DynamicBC_RegroupDat(InMatrix,InMod,sortinfo)
% matsize = size(InMatrix,1);
% OrigOrder = 1:matsize;
% [ix,iy] = sort(InMod);
% Re1iy = iy;
% OutMod = ix;
Inmodind = unique(InMod);
roire = zeros(length(InMod)+length(Inmodind),1);
Outsortinfo = zeros(length(InMod)+length(Inmodind),1);
st = 1;
OutMod = zeros(length(InMod)+length(Inmodind),1);
for i = 1:length(Inmodind)
    tempind = find(InMod==Inmodind(i));
    sortinfotemp = sortinfo(tempind);
%     iytemp = iy(tempind);
    [ixt,iyt] = sort(sortinfotemp,'descend');
    roire(st:st+length(sortinfotemp)-1) = tempind(iyt);
    OutMod(st:st+length(sortinfotemp)-1) = Inmodind(i);
    Outsortinfo(st:st+length(sortinfotemp)-1) = ixt;
    ZeroInd(i) = st+length(sortinfotemp);
    ReGroupInd{i} = st:st+length(sortinfotemp)-1;
    st = st+length(sortinfotemp)+1;
end
OutMatrix1 = zeros(length(InMod)+length(Inmodind),length(InMod));
OutMatrix = zeros(length(InMod)+length(Inmodind),length(InMod)+length(Inmodind));
for i = 1:length(Inmodind)
    OutMatrix1(ReGroupInd{i},:) = InMatrix(roire(ReGroupInd{i}),:);
end
for i = 1:length(Inmodind)
    OutMatrix(:,ReGroupInd{i}) = OutMatrix1(:,roire(ReGroupInd{i}));
end

end
