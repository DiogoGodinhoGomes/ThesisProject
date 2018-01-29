openr,1,'ref.in'
h=strarr(1,31)
readf,1,h
openw,2,'atm.in'
printf,2,h
atm=fltarr(7,5)
readf,1,atm
close,1
i=0.
j=1.
k=0.
l=0.
read,i
read,j
read,k
read,l
read,m
read,n

tshift=[-500,0,500]
pres0=[3.,1.]
pres1=[0.3,0.1,0.01]
pfac2=[0.1,0.01]


h2ovmr=10.^((float(l)*0.75) -6.)

atm[3,1]=pres0[fix(i)-1]
atm[3,2]=pres1[fix(j)-1]
atm[3,3]=atm[3,2]*pfac2[fix(k)-1]

atm[4,2]=atm[4,2]+tshift[fix(m)-1]
atm[4,3]=atm[4,3]+tshift[fix(n)-1]
atm[4,4]=atm[4,3]

atm[6,*]=h2ovmr
;atm[7,*]=ch4vmr




printf,2,atm,format='(i2,d6.3,i2,e12.5,d7.1,2e12.4)'

close,2

end
