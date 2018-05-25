lower = -vsbar;
peak = 0;
upper = vsbar;

N=10000;
pd = makedist('Triangular','a',lower,'b',peak,'c',upper);
rng('default');
r = random(pd,N,1);
hist(r)