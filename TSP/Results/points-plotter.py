from cmath import sqrt
import matplotlib.pyplot as plt
import os

cityMap = {}
cities = "210.txt"
dir = os.path.dirname(os.path.abspath(__file__))

with open(dir + "\\" + cities, "r") as file:
    for city in file:
        num, posx, posy = city.split(" ")
        num = int(num)
        posx = int(posx)
        posy = int(posy[:-1])
        cityMap[num] = [posx, posy]

x = []
y = []
n = []

results = "obrazek-para.txt"
path = 0
with open(dir + "\\" + results, "r") as file:
    for block in file:
        numbers = block.split(" ")
        numbers = [x for x in numbers if x != '']

        for num in numbers:
            num = int(num)
            posx, posy = cityMap[num]
            if len(x) > 0 and len(y) > 0:
                path += sqrt(pow(posx - x[-1], 2) + pow(posy - y[-1], 2))

            x.append(posx)
            y.append(posy)
            n.append(num + 1)

x.append(x[0])
y.append(y[0])
path += sqrt(pow(x[0] - x[-1], 2) + pow(y[0] - y[-1], 2))
print(path)

plt.figure(figsize=(10, 8))
plt.xlabel('x')
plt.ylabel('y')
plt.plot(x,y, lw=1 , zorder=1)
plt.scatter(x,y,s=30, zorder=2)
plt.scatter(x[0],y[0],s=60,c='red',zorder=3)    

'''
for i, txt in enumerate(n):
    if txt == 5 or txt == 6:
        plt.text(x[i]-10,y[i]-60, txt, fontsize=15)
    else:
        plt.text(x[i]+10,y[i], txt, fontsize=15)
'''

plt.show()