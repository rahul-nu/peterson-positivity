import itertools
from multiprocessing import Pool

lie_type='E'
rank=8

if (lie_type=='A'): pexponents=prod([1..rank])
if (lie_type in ['B','C']): pexponents=prod([1,3,..2*rank-1])
if (lie_type=='D'): pexponents=prod([1,3,..2*rank-3]+[rank-1])
if (lie_type=='F'): pexponents=prod([1,5,7,11])
if (lie_type=='G'): pexponents=prod([1,5])
if (lie_type=='E' and rank==6): pexponents=prod([1,4,5,7,8,11])
if (lie_type=='E' and rank==7): pexponents=prod([1,5,7,9,11,13,17])
if (lie_type=='E' and rank==8): pexponents=prod([1,7,11,13,17,19,23,29])

R = RootSystem([lie_type,rank]).root_space()
W = WeylGroup(R,prefix="s") 
s = W.simple_reflections()
alpha = R.simple_roots()
Lambdacheck = R.root_system.dual.weight_space().fundamental_weights()
w=W.long_element()       
wWord=w.reduced_word()
partialHts=[]

def to_number(root): return (sum(Lambdacheck[k].scalar(root) for k in range(1,rank+1)))

u=s[1]*s[1]
for i in wWord: 
    partialHts.append(to_number(u.action(alpha[i])))
    u*=s[i]

#    partialWords.append(partialWords[-1]*s[i])
#v=prod(s[i] for i in v) 


def NvI(v):
    lenv=rank
    reducedWords=[tuple(x) for x in v.reduced_words()]

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
        p=1    
        for i in subword:
            p*=partialHts[i] #to_number(partialWords[i].action(alpha[wWord[i]]))
        q+=p

    return(q/pexponents)

#print("Starting")
#counter=0
#perms=itertools.permutations([1..rank])
#vall=list({prod(s[i] for i in x) for x in perms})
#vall=list(vall)
#vall=[(i,vall[i]) for i in range(len(vall))]

vall=[(0, s[8]*s[7]*s[2]*s[4]*s[5]*s[6]*s[3]*s[1]),
 (1, s[8]*s[7]*s[6]*s[3]*s[2]*s[4]*s[5]*s[1]),
 (2, s[8]*s[6]*s[7]*s[2]*s[4]*s[5]*s[3]*s[1]),
 (3, s[8]*s[6]*s[7]*s[5]*s[4]*s[3]*s[2]*s[1]),
 (4, s[1]*s[3]*s[4]*s[5]*s[6]*s[7]*s[8]*s[2]),
 (5, s[7]*s[8]*s[5]*s[6]*s[4]*s[3]*s[2]*s[1]),
 (6, s[8]*s[7]*s[6]*s[3]*s[4]*s[5]*s[2]*s[1]),
 (7, s[8]*s[7]*s[6]*s[4]*s[5]*s[1]*s[3]*s[2]),
 (8, s[8]*s[1]*s[3]*s[4]*s[5]*s[6]*s[7]*s[2]),
 (9, s[6]*s[7]*s[8]*s[4]*s[5]*s[3]*s[2]*s[1]),
 (10, s[6]*s[7]*s[8]*s[5]*s[3]*s[2]*s[4]*s[1]),
 (11, s[8]*s[6]*s[7]*s[1]*s[3]*s[2]*s[4]*s[5]),
 (12, s[5]*s[6]*s[7]*s[8]*s[1]*s[3]*s[2]*s[4]),
 (13, s[3]*s[2]*s[4]*s[5]*s[6]*s[7]*s[8]*s[1]),
 (14, s[2]*s[4]*s[5]*s[6]*s[7]*s[8]*s[3]*s[1]),
 (15, s[8]*s[1]*s[3]*s[2]*s[4]*s[5]*s[6]*s[7]),
 (16, s[7]*s[8]*s[1]*s[3]*s[2]*s[4]*s[5]*s[6]),
 (17, s[8]*s[7]*s[5]*s[6]*s[3]*s[4]*s[2]*s[1]),
 (18, s[8]*s[7]*s[5]*s[6]*s[1]*s[3]*s[4]*s[2]),
 (19, s[5]*s[6]*s[7]*s[8]*s[3]*s[2]*s[4]*s[1]),
 (20, s[8]*s[4]*s[5]*s[6]*s[7]*s[1]*s[3]*s[2]),
 (21, s[5]*s[6]*s[7]*s[8]*s[4]*s[3]*s[2]*s[1]),
 (22, s[6]*s[7]*s[8]*s[3]*s[2]*s[4]*s[5]*s[1]),
 (23, s[4]*s[5]*s[6]*s[7]*s[8]*s[3]*s[2]*s[1]),
 (24, s[7]*s[8]*s[6]*s[2]*s[4]*s[5]*s[3]*s[1]),
 (25, s[6]*s[7]*s[8]*s[3]*s[4]*s[5]*s[2]*s[1]),
 (26, s[8]*s[5]*s[6]*s[7]*s[2]*s[4]*s[3]*s[1]),
 (27, s[6]*s[7]*s[8]*s[4]*s[5]*s[1]*s[3]*s[2]),
 (28, s[8]*s[7]*s[6]*s[1]*s[3]*s[4]*s[5]*s[2]),
 (29, s[8]*s[7]*s[6]*s[5]*s[2]*s[4]*s[1]*s[3]),
 (30, s[8]*s[7]*s[1]*s[3]*s[4]*s[5]*s[6]*s[2]),
 (31, s[8]*s[5]*s[6]*s[7]*s[4]*s[1]*s[3]*s[2]),
 (32, s[7]*s[8]*s[6]*s[1]*s[3]*s[2]*s[4]*s[5]),
 (33, s[8]*s[7]*s[4]*s[5]*s[6]*s[3]*s[2]*s[1]),
 (34, s[7]*s[8]*s[5]*s[6]*s[2]*s[4]*s[3]*s[1]),
 (35, s[8]*s[7]*s[5]*s[6]*s[2]*s[4]*s[1]*s[3]),
 (36, s[7]*s[8]*s[5]*s[6]*s[4]*s[1]*s[3]*s[2]),
 (37, s[7]*s[8]*s[3]*s[4]*s[5]*s[6]*s[2]*s[1]),
 (38, s[8]*s[6]*s[7]*s[5]*s[1]*s[3]*s[2]*s[4]),
 (39, s[8]*s[7]*s[2]*s[4]*s[5]*s[6]*s[1]*s[3]),
 (40, s[8]*s[6]*s[7]*s[2]*s[4]*s[5]*s[1]*s[3]),
 (41, s[4]*s[5]*s[6]*s[7]*s[8]*s[1]*s[3]*s[2]),
 (42, s[7]*s[8]*s[6]*s[5]*s[2]*s[4]*s[3]*s[1]),
 (43, s[7]*s[8]*s[6]*s[5]*s[4]*s[3]*s[2]*s[1]),
 (44, s[8]*s[2]*s[4]*s[5]*s[6]*s[7]*s[3]*s[1]),
 (45, s[6]*s[7]*s[8]*s[1]*s[3]*s[4]*s[5]*s[2]),
 (46, s[8]*s[6]*s[7]*s[4]*s[5]*s[3]*s[2]*s[1]),
 (47, s[8]*s[6]*s[7]*s[5]*s[3]*s[2]*s[4]*s[1]),
 (48, s[5]*s[6]*s[7]*s[8]*s[2]*s[4]*s[3]*s[1]),
 (49, s[8]*s[7]*s[4]*s[5]*s[6]*s[1]*s[3]*s[2]),
 (50, s[7]*s[8]*s[3]*s[2]*s[4]*s[5]*s[6]*s[1]),
 (51, s[8]*s[7]*s[5]*s[6]*s[1]*s[3]*s[2]*s[4]),
 (52, s[5]*s[6]*s[7]*s[8]*s[4]*s[1]*s[3]*s[2]),
 (53, s[7]*s[8]*s[2]*s[4]*s[5]*s[6]*s[3]*s[1]),
 (54, s[2]*s[4]*s[5]*s[6]*s[7]*s[8]*s[1]*s[3]),
 (55, s[8]*s[7]*s[6]*s[2]*s[4]*s[5]*s[3]*s[1]),
 (56, s[7]*s[8]*s[6]*s[5]*s[3]*s[4]*s[2]*s[1]),
 (57, s[7]*s[8]*s[6]*s[5]*s[1]*s[3]*s[4]*s[2]),
 (58, s[6]*s[7]*s[8]*s[5]*s[4]*s[1]*s[3]*s[2]),
 (59, s[6]*s[7]*s[8]*s[5]*s[2]*s[4]*s[3]*s[1]),
 (60, s[8]*s[7]*s[6]*s[1]*s[3]*s[2]*s[4]*s[5]),
 (61, s[8]*s[7]*s[6]*s[5]*s[4]*s[3]*s[2]*s[1]),
 (62, s[8]*s[7]*s[5]*s[6]*s[3]*s[2]*s[4]*s[1]),
 (63, s[8]*s[5]*s[6]*s[7]*s[3]*s[4]*s[2]*s[1]),
 (64, s[8]*s[5]*s[6]*s[7]*s[1]*s[3]*s[4]*s[2]),
 (65, s[6]*s[7]*s[8]*s[5]*s[3]*s[4]*s[2]*s[1]),
 (66, s[6]*s[7]*s[8]*s[5]*s[1]*s[3]*s[4]*s[2]),
 (67, s[6]*s[7]*s[8]*s[5]*s[2]*s[4]*s[1]*s[3]),
 (68, s[7]*s[8]*s[6]*s[2]*s[4]*s[5]*s[1]*s[3]),
 (69, s[8]*s[6]*s[7]*s[3]*s[2]*s[4]*s[5]*s[1]),
 (70, s[8]*s[7]*s[5]*s[6]*s[4]*s[3]*s[2]*s[1]),
 (71, s[7]*s[8]*s[6]*s[5]*s[4]*s[1]*s[3]*s[2]),
 (72, s[8]*s[6]*s[7]*s[3]*s[4]*s[5]*s[2]*s[1]),
 (73, s[8]*s[6]*s[7]*s[4]*s[5]*s[1]*s[3]*s[2]),
 (74, s[7]*s[8]*s[5]*s[6]*s[3]*s[4]*s[2]*s[1]),
 (75, s[7]*s[8]*s[5]*s[6]*s[1]*s[3]*s[4]*s[2]),
 (76, s[7]*s[8]*s[6]*s[4]*s[5]*s[3]*s[2]*s[1]),
 (77, s[1]*s[3]*s[2]*s[4]*s[5]*s[6]*s[7]*s[8]),
 (78, s[8]*s[7]*s[6]*s[5]*s[4]*s[1]*s[3]*s[2]),
 (79, s[6]*s[7]*s[8]*s[5]*s[4]*s[3]*s[2]*s[1]),
 (80, s[6]*s[7]*s[8]*s[2]*s[4]*s[5]*s[3]*s[1]),
 (81, s[8]*s[7]*s[1]*s[3]*s[2]*s[4]*s[5]*s[6]),
 (82, s[8]*s[5]*s[6]*s[7]*s[2]*s[4]*s[1]*s[3]),
 (83, s[8]*s[3]*s[2]*s[4]*s[5]*s[6]*s[7]*s[1]),
 (84, s[6]*s[7]*s[8]*s[1]*s[3]*s[2]*s[4]*s[5]),
 (85, s[7]*s[8]*s[1]*s[3]*s[4]*s[5]*s[6]*s[2]),
 (86, s[7]*s[8]*s[4]*s[5]*s[6]*s[3]*s[2]*s[1]),
 (87, s[8]*s[2]*s[4]*s[5]*s[6]*s[7]*s[1]*s[3]),
 (88, s[7]*s[8]*s[6]*s[5]*s[2]*s[4]*s[1]*s[3]),
 (89, s[7]*s[8]*s[5]*s[6]*s[2]*s[4]*s[1]*s[3]),
 (90, s[7]*s[8]*s[6]*s[3]*s[2]*s[4]*s[5]*s[1]),
 (91, s[5]*s[6]*s[7]*s[8]*s[3]*s[4]*s[2]*s[1]),
 (92, s[5]*s[6]*s[7]*s[8]*s[1]*s[3]*s[4]*s[2]),
 (93, s[7]*s[8]*s[6]*s[3]*s[4]*s[5]*s[2]*s[1]),
 (94, s[7]*s[8]*s[6]*s[5]*s[1]*s[3]*s[2]*s[4]),
 (95, s[8]*s[6]*s[7]*s[1]*s[3]*s[4]*s[5]*s[2]),
 (96, s[7]*s[8]*s[6]*s[4]*s[5]*s[1]*s[3]*s[2]),
 (97, s[3]*s[4]*s[5]*s[6]*s[7]*s[8]*s[2]*s[1]),
 (98, s[7]*s[8]*s[2]*s[4]*s[5]*s[6]*s[1]*s[3]),
 (99, s[8]*s[7]*s[6]*s[2]*s[4]*s[5]*s[1]*s[3]),
 (100, s[8]*s[7]*s[5]*s[6]*s[2]*s[4]*s[3]*s[1]),
 (101, s[7]*s[8]*s[6]*s[5]*s[3]*s[2]*s[4]*s[1]),
 (102, s[8]*s[7]*s[5]*s[6]*s[4]*s[1]*s[3]*s[2]),
 (103, s[8]*s[4]*s[5]*s[6]*s[7]*s[3]*s[2]*s[1]),
 (104, s[8]*s[7]*s[6]*s[5]*s[1]*s[3]*s[2]*s[4]),
 (105, s[8]*s[7]*s[3]*s[4]*s[5]*s[6]*s[2]*s[1]),
 (106, s[8]*s[5]*s[6]*s[7]*s[1]*s[3]*s[2]*s[4]),
 (107, s[8]*s[7]*s[6]*s[5]*s[3]*s[4]*s[2]*s[1]),
 (108, s[8]*s[7]*s[6]*s[4]*s[5]*s[3]*s[2]*s[1]),
 (109, s[8]*s[7]*s[6]*s[5]*s[1]*s[3]*s[4]*s[2]),
 (110, s[8]*s[6]*s[7]*s[5]*s[4]*s[1]*s[3]*s[2]),
 (111, s[8]*s[6]*s[7]*s[5]*s[2]*s[4]*s[3]*s[1]),
 (112, s[7]*s[8]*s[4]*s[5]*s[6]*s[1]*s[3]*s[2]),
 (113, s[5]*s[6]*s[7]*s[8]*s[2]*s[4]*s[1]*s[3]),
 (114, s[8]*s[6]*s[7]*s[5]*s[3]*s[4]*s[2]*s[1]),
 (115, s[8]*s[6]*s[7]*s[5]*s[1]*s[3]*s[4]*s[2]),
 (116, s[8]*s[3]*s[4]*s[5]*s[6]*s[7]*s[2]*s[1]),
 (117, s[8]*s[6]*s[7]*s[5]*s[2]*s[4]*s[1]*s[3]),
 (118, s[7]*s[8]*s[5]*s[6]*s[1]*s[3]*s[2]*s[4]),
 (119, s[8]*s[7]*s[6]*s[5]*s[3]*s[2]*s[4]*s[1]),
 (120, s[8]*s[5]*s[6]*s[7]*s[3]*s[2]*s[4]*s[1]),
 (121, s[8]*s[7]*s[6]*s[5]*s[2]*s[4]*s[3]*s[1]),
 (122, s[8]*s[5]*s[6]*s[7]*s[4]*s[3]*s[2]*s[1]),
 (123, s[7]*s[8]*s[6]*s[1]*s[3]*s[4]*s[5]*s[2]),
 (124, s[6]*s[7]*s[8]*s[5]*s[1]*s[3]*s[2]*s[4]),
 (125, s[8]*s[7]*s[3]*s[2]*s[4]*s[5]*s[6]*s[1]),
 (126, s[7]*s[8]*s[5]*s[6]*s[3]*s[2]*s[4]*s[1]),
 (127, s[6]*s[7]*s[8]*s[2]*s[4]*s[5]*s[1]*s[3])]

#print(len(vall))
#for v in vall: 
#    counter+=1
#    print ((counter,NvI(v)/len(v.reduced_words())))

#def myFunc(iv):
#    (i,v)=iv
#    print ((i,NvI(v)/len(v.reduced_words())))
#
#pool=Pool(4)
#pool.map(myFunc,vall)

#
#
#u=s[1]*s[1]
#for i in range(1,rank+1):
#    for j in range(1,rank+1):
#        print(to_number(u.action(alpha[j])),end=' ')   
#        u*=s[j]
#    print('')
    
#u=u*u
#(u*s[1]*s[2]).action(alpha[3])

