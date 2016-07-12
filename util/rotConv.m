function Rb = rotConv(qb_a, Ra)

Rb = Ra;

type = 0;
if size(Ra,1) == size(Ra,2)
    type = 1; % Rotation Matrix
    qa = dcm2quat(Ra);
elseif max(size(Ra)) == 3 && min(size(Ra)) == 1
    type = 2; % Euler angles
    qa = eulr2quat(Ra);
end

if type > 0
    q2c = quatconj(qb_a);

    qb = quatprod(qb_a, quatprod(qa, q2c));
    if type == 1
        Rb = quat2dcm(qb);
    elseif type == 2
        Rb = quat2eulr(qb);
    end
end
