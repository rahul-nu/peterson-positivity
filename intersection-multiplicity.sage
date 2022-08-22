R.<a,b,c>=PolynomialRing(QQ,3)
R.<e,d,c,b,a>=PolynomialRing(QQ,5)

t1=a
t2=a^2/2
t3=a^3/6+b
t4=a^4/24+a*b
t5=a^5/120+a^2*b/2+c
t6=a^6/720+a^3*b/3+a*c+b^2/2
t7=a^7/(7*720)+a^4*c/24+a^2*c/2+a*b^2/2+d
t8=a^8/(56*720)+a^5/120*b+a^3/6*c+a^2*b^2/4+a*d+b*c
t9=a^9/(504*720)+a^6*b/720+a^4*c/24+a^3*b^2/2+a^2*d/2+a*b*c+b^3/6+e

I=Ideal(t4,t5,t6,t7)
I.vector_space_dimension()

#t6-t5*t1-1/2*t3^2+1/3*t1^2*t4+1/144*t1^6=0
#t4==t1*t3-1/8*t1^4

