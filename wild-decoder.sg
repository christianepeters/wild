# Wild-Goppa-code decoder
# as in "Wild McEliece" 
# by Daniel J. Bernstein, Tanja Lange, Christiane Peters
# To appear at SAC 2010.
# Article: http://www2.mat.dtu.dk/people/C.Peters/publications/2010.wild.pdf
#
# The following is an example. Note that it's from the perspective 
# of the code Gamma_q^m(a1,...,an, g^q). 
# I didn't bother looking for a codeword in the subfield Fq.
#
# By Christiane Peters, August 10, 2010.
# Questions -> c.p.peters@tue.nl
#
#
q=3; K=GF(q)
m=5; n=q**m; t=8
F.<a>=GF(q^m)
R.<x>=PolynomialRing(F)
#g=R.random_element(t); g.is_irreducible()
g=x^8 + (2*a^3 + 2*a + 2)*x^7 + (a^4 + 2*a^2 + a)*x^6 + (2*a + 1)*x^5\
+ (a^3 + a^2)*x^4 + (a^2 + a + 2)*x^3 + (2*a^4 + a^3 + a)*x^2 + (2*a^3\
+ 2*a^2 + 2*a)*x + a^4 + 2*a^3 + a
g.is_irreducible()
gq=g^q
w=floor(gq.degree()/2)
print "g^q=",gq
L=[a^i for i in range(n-1)]+[0]
h=x^n-x
dh=h.derivative() # with this choice of h the derivative's equal to -1
f=R.random_element(n-gq.degree()-1)*gq
V=VectorSpace(F,n)
c=vector([f(ai)/dh(ai) for ai in L])
e=vector([0 for i in range(n)])
for i in range(w): 
  e[randint(0,n-1)]=randint(1,q-1)

y=c+e
points = []
for i in range(n) : 
  points=points + [(L[i],y[i]*dh/gq(L[i]))]

print "---> starting Lagrange interpolation"
time phi=R.lagrange_polynomial(points)

print "---> starting continued fraction computation"
# continued fractions 
def eea(a,b,deg): 
  l0=vector([a,1,0]); l1=vector([b,0,1])
  i=1
  while l1[0].degree() >= deg :
    i=i+1
    q_=l0[0]//l1[0]
    l2=l0-q_*l1
    l0=l1; l1=l2
  return l1

time l=eea(h,phi,n-w)
l[0].degree()-(n-w)
v0=l[1]; v1=-l[2]
l[0]-(v0*h-v1*phi)
(phi-v0*h/v1)*gq-f
print f-(phi-v0*h/v1)*gq ==0

# cheating, i.e. check against the error-locator polynomial epsilon:
#epsilon=prod([(x-L[i]) for i in range(n) if e[i]!=0])
#delta=(phi-f/gq)/h*epsilon
#delta*v1-epsilon*v0
#phi-f/gq- h/epsilon*delta
#phi-f/gq-h/v1*v0
#(phi-v0*h/v1)*gq-f
