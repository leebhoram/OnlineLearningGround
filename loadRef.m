% prepare_reference_data
ref_data = load(posefile);
N_pose = length(ref_data);
ref_path = ref_data(:,[4 8 12]);
ref_rot = ref_data(:,[1:3 5:7 9:11]);
ref_eulr = zeros(N,3);
Ca_b = [0 0 1; 1 0 0; 0 1 0];
qb_a = dcm2quat(Ca_b');
qa_b = quatconj(qb_a);

for k=1:N
    Rb = rotConv(qb_a, reshape(ref_rot(k,:),3,3)');
    ref_eulr(k,:) = dcm2eulr(Rb); 
end

clear Ca_b qb_a qb_a