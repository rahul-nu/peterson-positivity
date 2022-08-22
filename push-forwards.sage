from itertools import chain, combinations

lie_type='B'
rank=3

R = RootSystem([lie_type,rank]).root_space()
W = WeylGroup(R,prefix="s")           #,implementation='coxeter3')
s = W.simple_reflections()
alpha = R.simple_roots()
Lambdacheck = R.root_system.dual.weight_space().fundamental_weights()
varlist = ','.join('a%s' % i for i in range(1,rank+1))
A = PolynomialRing(ZZ,varlist)
a = {k: A.gens()[k-1] for k in alpha.keys()} 
e=s[1]*s[1]

def to_poly(root):
    result = 0
    for k in a.keys():
        result += Lambdacheck[k].scalar(root)*a[k]
    return result

def localize(v,w):
    if (v==s[1]*s[1]): return 1
    lenv=v.length()
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
#    T=A.hom(['t']*rank,ZZ['t'])
    T=A.hom([1]*rank,ZZ)
    return(T(q))
#    return(q)


coxAll=list(chain.from_iterable(combinations([1..rank], r+1) for r in range(rank)))
coxAll=[e]+[prod(s[i] for i in x) for x in coxAll]

wAll=coxAll.copy() #w for w in W.elements_of_length(0)]
#for i in range(rank+1): wAll.extend([w for w in W.elements_of_length(i) if not (len(w.support())==w.length())])
wAll.extend([w for w in W if w not in coxAll])

#maxAll=[e,s[2],s[3],s[4],s[2]*s[3]*s[2],s[2]*s[4],s[4]*s[3]*s[4]*s[3], s[4]*s[3]*s[4]*s[2]*s[3]*s[4]*s[2]*s[3]*s[2],s[1],s[1]*s[2]*s[1],s[1]*s[3],s[1]*s[4],s[1]*s[2]*s[3]*s[1]*s[2]*s[1],s[1]*s[2]*s[1]*s[4],s[1]*s[4]*s[3]*s[4]*s[3],W.long_element()]
maxAll=[e,s[1],s[2],s[3],s[1]*s[2]*s[1],s[1]*s[3],s[3]*s[2]*s[3]*s[2], s[3]*s[2]*s[3]*s[1]*s[2]*s[3]*s[1]*s[2]*s[1]]
M=[]
for i in range(len(wAll)):
    N=[]
    for j in range(len(maxAll)): 
#        print(wAll[i],maxAll[j])
        N.append(localize(wAll[i],maxAll[j]))
#        print(N)
    M.append(N)

B=Matrix([M[i] for i in range(2^rank)]) #[[X^v]]_{v\in Cox}=B*[[w_I]_S]_{I\subset\Delta}
C=Matrix([M[i] for i in range(2^rank,len(wAll))]) #[[X^v]]_{v\in W}=C*[[w_I]_S]_{I\subset\Delta}
print(wAll[2^rank:])
D=C*B.inverse() #[[p_v]]_{v\in wAll\Cox}=D[[p_v]]_{v\in Cox}
print(D)
print(coxAll)

    
