function I_double2single = generate_projection_operator(Bs)    
    I_double2single = zeros(2*Bs,2*Bs-1);
    for i=1:Bs-2
        I_double2single(i,i) = 1;
    end
    for i=Bs+3:2*Bs
        I_double2single(i,i-1) = 1;
    end
    I_double2single(Bs-1,Bs-3) = -3/59;
    I_double2single(Bs-1,Bs-2) = 10/59;
    I_double2single(Bs-1,Bs-1) = 48/59;
    I_double2single(Bs-1,Bs) = 4/59;
    I_double2single(Bs,Bs-3) = 2/17;
    I_double2single(Bs,Bs-2) = -5/17;
    I_double2single(Bs,Bs-1) = 2/17;
    I_double2single(Bs,Bs) = 20/17;
    I_double2single(Bs,Bs+1) = -2/17;
    I_double2single(Bs+1,Bs-1) = -2/17;
    I_double2single(Bs+1,Bs-1:Bs+3) = I_double2single(Bs,Bs+1:-1:Bs-3) ;
    I_double2single(Bs+2,Bs:Bs+3) = I_double2single(Bs-1,Bs:-1:Bs-3) ; 

