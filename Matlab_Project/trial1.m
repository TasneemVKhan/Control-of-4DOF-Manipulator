s = tf ('s');
G1 = 1/(s+1);
G2 = (10*s+12)/(s);
G3 = (20*s);
H1 = (s)/((s)*(s+1));
G23 = series(G1,G2);
HG23 = feedback(G23,H1);
OutputTF = series(G1,HG23);
OutputTF
rlocus(OutputTF);
figure(1);

%pzmap(OutputTF);
%figure(2);
%grid on
