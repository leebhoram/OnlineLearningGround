% by Bhoram Lee 
function [ quat ] = dcm2quat( dcm )

eulr = dcm2eulr(dcm);
quat = eulr2quat(eulr);

end

