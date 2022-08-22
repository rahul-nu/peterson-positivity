import itertools

lie_type='B'
rank=4

R = RootSystem([lie_type,rank]).root_space()
W = WeylGroup(R,prefix="s") 
s = W.simple_reflections()
alpha = R.simple_roots()
Lambdacheck = R.root_system.dual.weight_space().fundamental_weights()
w=W.long_element()       
#wTry=prod(s[i] for i in [1..rank]*rank)
#if (wTry==w): w=wTry
#else: print("Fail")
#wWord=w.reduced_word()

def to_number(root): return (sum(Lambdacheck[k].scalar(root) for k in range(1,rank+1)))

#wordClass=[ [1, 2, 3, 2, 1, 2],
# [3, 2, 1, 2, 3, 2],
# [2, 1, 2, 3, 2, 1],
# [2, 1, 3, 2, 1, 3],
# [1, 2, 1, 3, 2, 1],
# [2, 3, 2, 1, 2, 3],
# [3, 2, 1, 3, 2, 3],
# [1, 3, 2, 1, 3, 2]]

def wt(word):
    wt=1
    for i in range(len(word)-1): 
        if ({word[i],word[i+1]}=={1,3}): wt*=2
    return wt

perms=itertools.permutations([1..rank])
allCox={prod(s[i] for i in perm) for perm in perms}

myWord=[1,2,3,4,3,2,1,2,3,4,3,2,3,4,3,4]
myWord.reverse()

#for wWord in [w.reduced_word()]: #[[1..rank]*rank]:
#for wWord in [[1..rank]*rank]:
for wWord in [myWord]:
#    print(wWord)
    ans=True
    for v in allCox: #[(1,2,4,3),(1,4,2,3),(4,1,2,3)]:
        nv=0
        rv=len(v.reduced_words())
        for vi in v.reduced_words():
    #        tot=0
        #    for wWord in w.reduced_words():
            outerList=[[]]*rank
            positionsOf={k:[] for k in range(1,rank+1)}
            for pos in range(len(wWord)): positionsOf[wWord[pos]].append(pos)
            for i in range(rank): outerList[i]=positionsOf[vi[i]]

            n=[[x] for x in outerList[0]]
            for innerList in outerList[1:]:
                n=[ni+[b] for ni in n for b in innerList if ni[-1]<b]
            nv+=len(n)
        print(nv,rv)
#        ans=ans and (nv==rv)
#        if (not ans):break
#    if (ans): print(wWord,ans)
#    print(wWord,ans)

#        q=0
#        for subword in n:
#            oi=0
#            u=[x for x in W.elements_of_length(0)][0]
#            p=1    
#            for i in subword:
#                if (i!=0): u=u*prod(s[wWord[k]] for k in range(oi,i))
#                p*=to_number(u.action(alpha[wWord[i]]))
#                oi=i
#            q+=p
#        tot+=q/wt(wWord)
#    print((vi,tot))



