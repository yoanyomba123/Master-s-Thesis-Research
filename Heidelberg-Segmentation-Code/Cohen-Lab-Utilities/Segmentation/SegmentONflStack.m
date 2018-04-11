function [onfl] = SegmentOnflStack(BScans,infl,ipl,icl,opl,bv,Params)


for i = 1: size(BScans, 3)
    onfl_res = SegmentOnfl(BScans(:,:,i),infl(i,:),ipl(i,:),icl(i,:), opl(i,:),bv(i,:),Params);
    onfl(i, :) = onfl_res;
end

end

