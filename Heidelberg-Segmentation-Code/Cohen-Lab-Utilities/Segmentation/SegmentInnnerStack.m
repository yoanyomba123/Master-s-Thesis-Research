function [ipl,opl, icl] = SegmentInnnerStack(BScans,rpe, infl, med,bv, Params)
for i = 1: size(BScans, 3)
    [ipl_res, opl_res, icl_res] = SegmentInner(BScans(:,:,i),rpe(i,:),infl(i,:), med(i,:),bv(i,:),Params)
    ipl(i, :) = ipl_res;
    opl(i, :) = opl_res;
    icl(i, :) = icl_res;
end
end

