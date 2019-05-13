function Na_new=Unimolecular_Decay_Population_Only(Na,Pf)


r1=rand(Na,1);
dtemp=(r1>Pf);

n=sum(dtemp);

Na_new=n;
