from itertools import combinations

W=WeylGroup(['D',4],prefix="s")
w=W.long_element()
ref=W.reflections()
[s1,s2,s3,s4]=W.simple_reflections()
v1=s2*s1*s3*s4
v2=s1*s2*s3*s4
v3=s1*s3*s2*s4
v4=s1*s3*s4*s2
v=[s2*s1*s3*s4,s1*s2*s3*s4,s1*s3*s2*s4,s1*s3*s4*s2]

[[s[i].bruhat_le(w*v[j]) for i in range(4)] for j in range(4)]

R=[(2,1,4,3),(2,1,3,4),(2,3,1,4),(2,3,4,1),(2,4,1,3),(2,4,3,1)]
full=w.reduced_word()

