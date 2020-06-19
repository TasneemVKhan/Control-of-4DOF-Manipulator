s = tf ('s');
G1 = 1/(s+1);
G2 = (10*s+12)/(s);
G3 = ((13.68*s^3)+(59.88*s))/((s^2+4.37));
H1 = (s)/((s)*(s+1));
G23 = series(G1,G2);
HG23 = feedback(G23,H1);
OutputTF = series(G1,HG23);
margin(OutputTF);