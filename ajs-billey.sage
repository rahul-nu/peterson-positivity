#Program to compute [X^v]|_w
#e7=661620960
#e7/(1*5*7*9*11*13*17)= 864
#e8=11179629901440
#e8/prod([7,11,13,17,19,23,29])=51840
#import itertools

# General root system setup
#def mkRing(rank): # Generate the ring A with variables corresponding to the simple roots
#    varlist = ','.join('a%s' % i for i in range(1,rank+1))
#    A = PolynomialRing(ZZ,varlist)
#    a = {k: A.gens()[k-1] for k in alpha.keys()} # a[i] corresponds to the simple root alpha[i],
# subs[1] prints a list of [s[1].action(alpha[1]), s[1].action(alpha[2]),...], etc..
#    subs = {k: [to_poly(s[k].action(alpha[j])) for j in alpha.keys()] for k in a.keys()}
#phi[i] is the extension to A of the action of as s[i] on the Root System
#    phi = {k: A.hom(subs[k],A) for k in alpha.keys()}
#    return(A)


def localize(lie_type,rank,v,w='l'):
    R = RootSystem([lie_type,rank]).root_space()
    W = WeylGroup(R,prefix="s") #,implementation='coxeter3')
    s = W.simple_reflections()
    alpha = R.simple_roots()
    Lambdacheck = R.root_system.dual.weight_space().fundamental_weights()
    if (w=='l'): w=W.long_element()       #point to localize at.
#    v=prod(s[i] for i in v)               #Schubert variety X^v to be localized 
    lenv=v.length()
    varlist = ','.join('a%s' % i for i in range(1,rank+1))
    A = PolynomialRing(ZZ,varlist)
    a = {k: A.gens()[k-1] for k in alpha.keys()} # a[i] corresponds to the simple root alpha[i],

# I need to be able to turn roots into polynomials in A in order to define the ring homomorphisms
    def to_poly(root):
        result = 0
        for k in a.keys():
            result += Lambdacheck[k].scalar(root)*a[k]
        return result

    wWord=w.reduced_word()
    reducedWords=[tuple(x) for x in v.reduced_words()]

    m={}
    subwords=[]
    for vi in reducedWords:
        outerList=[[]]*lenv
        positionsOf={k:[] for k in range(1,rank+1)}
        for pos in range(len(wWord)): positionsOf[wWord[pos]].append(pos)
        for i in range(lenv): outerList[i]=positionsOf[vi[i]]

        n=[[x] for x in outerList[0]]
        for innerList in outerList[1:]:
            n=[ni+[b] for ni in n for b in innerList if ni[-1]<b]
        subwords.extend(n)

    q=0
    for subword in subwords:
        oi=0
        u=[x for x in W.elements_of_length(0)][0]
        p=1    
        for i in subword:
            if (i!=0): u=u*prod(s[wWord[k]] for k in range(oi,i))
            p*=to_poly(u.action(alpha[wWord[i]]))
            oi=i
        q+=p

    return(q)

def NvI(lie_type,rank,v=()):
    if (v==()): v=[1..rank]
    varlist = ','.join('a%s' % i for i in range(1,rank+1))
    A = PolynomialRing(ZZ,varlist)

    if (lie_type=='A'): exponents=[1..rank]
    if (lie_type in ['B','C']): exponents=[1,3,..2*rank-1]
    if (lie_type=='D'): exponents=[1,3,..2*rank-3]+[rank-1]
    if (lie_type=='F'): exponents=[1,5,7,11]
    if (lie_type=='G'): exponents=[1,5]
    if (lie_type=='E' and rank==6): exponents=[1,4,5,7,8,11]
    if (lie_type=='E' and rank==7): exponents=[1,5,7,9,11,13,17]
    if (lie_type=='E' and rank==8): exponents=[1,7,11,13,17,19,23,29]

    q=localize(lie_type,rank,v,'l') #print(q)
#   T=A.hom(['t']*rank,PolynomialRing(ZZ,'t')) #    print(T(q))
    T=A.hom([1]*rank,ZZ)
#    print ('N('+repr(v)+') = '+repr(T(q)/prod(exponents)))
    return (T(q)/prod(exponents))


#for v in vall: print (NvI('E',7,v.reduced_word())/len(v.reduced_words()))
