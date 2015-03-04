function zt = solve_zeq(A,Zv)

B = fun_loop(A);
Q = fun_cutset(A);
Z = diag(Zv);

it = -Q(1,2:end)*([B(:,2:end)*Z; Q(2:end,2:end)]\[B(:,1);Q(2:end,1)]);

zt = 1/it;
if isa(zt,'sym')
   zt = simplify(zt); 
end

end





                                                                                                                            