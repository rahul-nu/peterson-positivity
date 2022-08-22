# General root system setup
lie_type = 'E'
rank = 7
R = RootSystem([lie_type,rank]).root_space()
W = WeylGroup(R,prefix="s")
alpha = R.simple_roots()
Lambdacheck = R.root_system.dual.weight_space().fundamental_weights()
s = W.simple_reflections()

# Automatically generate the ring A with variables corresponding to the simple roots
varlist = ','.join('a%s' % i for i in alpha.keys())
A = PolynomialRing(ZZ,varlist)
a = {k: A.gens()[k-1] for k in alpha.keys()} # so a[1] corresponds to alpha[1], a[2] to alpha[2], etc..

# I need to be able to turn roots into polynomials in A in order to define the ring homomorphisms
def to_poly(root):
    result = 0
    for k in a.keys():
        result += Lambdacheck[k].scalar(root)*a[k]
    return result

# subs[1] prints a list of [s[1].action(alpha[1]), s[1].action(alpha[2]),...], etc..
# The relationship with your current program is: your phi1 is A.hom(subs[1]), your phi2 is A.hom(subs[2]), etc..
subs = {k: [to_poly(s[k].action(alpha[j])) for j in alpha.keys()] for k in a.keys()}
# This is the same dict as in your program, just automatically generated now
phi = {k: A.hom(subs[k]) for k in alpha.keys()}

# an example showing this works
#print('substitution: a1 |--> -a1, a2 |--> 2a1+a2: a1 |--> %s' % phi[1](a[1]+a[2]^2))
#print('substitution: a1 |--> a1+a2, a2 |--> -a2: a1^2 + a2 |--> %s' % phi[2](a[1]^2+a[2]))
