import matplotlib.pyplot as plt

x = []
y = []
n = []

with open("25.txt",'r') as reader:
    for line in reader:
        c,a,b = map(int, line.split(" "))
        x.append(a)
        y.append(b)
        n.append(c+1)

x.append(x[0])
y.append(y[0])

plt.figure(figsize=(10, 8))
plt.xlabel('x')
plt.ylabel('y')
#plt.plot(x,y, lw=1 , zorder=1)
plt.scatter(x,y,s=30, zorder=2)
#plt.scatter(x[0],y[0],s=60,c='red',zorder=3)
'''
for i, txt in enumerate(n):
    if txt == 5 or txt == 6:
        plt.text(x[i]-10,y[i]-60, txt, fontsize=15)
    else:
        plt.text(x[i]+10,y[i], txt, fontsize=15)
'''       
plt.show()