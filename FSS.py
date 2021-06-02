#Implementation of the FSS algorithm proposed in
# C. J. A. Bastos Filho, F. B. de Lima Neto, A. J. C. C. Lins, A. I. S. Nascimento and M. P. Lima,
# "A novel search algorithm based on fish school behavior,"
# 2008 IEEE International Conference on Systems, Man and Cybernetics, 2008, pp. 2646-2651,
# doi: 10.1109/ICSMC.2008.4811695.


from random import random, randrange, uniform
import matplotlib.pyplot as plt
import numpy
from matplotlib.animation import FuncAnimation



def f(v):  #evaluation function
    global fitnessFunction
    if fitnessFunction == "Sphere":
        return 1 / (1 + Sphere(v))
    elif fitnessFunction == "Beale":
        return 1/(1 + Beale(v))
    elif fitnessFunction == "Booth":
        return 1/(1 + Booth(v))
    elif fitnessFunction == "Ackley":
        return 1/(1+Ackley(v))
    elif fitnessFunction == "Rosenbrock":
        return 1/(1+Rosenbrock(v))



#============== fitness functions ============#

def Sphere(v):
    res = 0
    for i in v:
        res += i ** 2
    return res

def Beale(v):
    return ((1.5 - v[0] + v[0] * v[1]) ** 2 + (2.25 - v[0] + v[0] * v[1] ** 2) ** 2 + (
                2.625 - v[0] + v[0] * v[1] ** 3) ** 2)

def Booth(v):
    return (v[0] + 2* v[1] -7)**2 + (2*v[0] + v[1] -5)**2

def Ackley(v):
    return -20*numpy.exp(-0.2*(0.5*(v[0]**2 +v[1]**2))**0.5)- numpy.exp(0.5*(numpy.cos(2*numpy.pi*v[0]) + numpy.cos(2* numpy.pi*v[1])))+ numpy.exp(1) + 20

def Rosenbrock(v):
        return 100*(v[1] - v[0]**2)**2 + (1-v[0])**2




#============== operators ============#

def feedingOperator():
    max = abs(f(x[0][1]) - f(x[0][0]))
    for fish in x:
        if abs(f(fish[1]) - f(fish[0])) > max:
            max = abs(f(fish[1]) - f(fish[0]))
    for i in range(len(x)):
        if max > 0:
            weight = W[i][1] + ((f(x[i][1]) - f(x[i][0])) / max)
            if weight >= weight_scale:
                W[i].append(weight_scale)
            elif weight <= 1:
                W[i].append(1)
            else:
                W[i].append(weight)
            del W[i][0]
            if W[i][1] >= breedingThreshold and not i in breedingCandidates:
                breedingCandidates.append(i)



def individualMovement():
    for fish in x:
        temp = []
        for i in range(2):
            temp.append(fish[1][i] + stepind * random() *(-1)**randrange(2))
        if f(temp) > f(fish[1]) and isInBoundaries(temp):
            fish.append(temp)
            del fish[0]
        else:
            fish.append(fish[1])
            del fish[0]


def collectiveInstinctiveMovement():
    sum2 = 0
    sum1 = [0, 0]
    for fish in x:
        sum2 += f(fish[1]) - f(fish[0])
        for i in range(2):
            sum1[i] += (fish[1][i] - fish[0][i]) * (f(fish[1]) - f(fish[0]))

    for fish in x:
        if sum2 != 0:
            x1 = fish[1][0] + sum1[0] / sum2
            x2 = fish[1][1] + sum1[1] / sum2
            if isInBoundaries([x1, x2]):
                fish.append([x1, x2])
                del fish[0]


def collectiveVolitiveMovement():
    sum1 = [0, 0]
    totalWeight = 0
    prevTotalWeight = 0
    for fish in range(len(x)):
        for i in range(2):
            sum1[i] += x[fish][1][i] * W[fish][1]
            totalWeight += W[fish][1]
            prevTotalWeight += W[fish][0]
    Bari = [sum1[0] / totalWeight, sum1[1] / totalWeight]   #calculate baricenter
    for fish in x:
        x1 = []
        for i in range(2):
            if totalWeight > prevTotalWeight:
                x1.append(fish[1][i] - stepvol * random() * (fish[1][i] - Bari[i]))
            else:
                x1.append(fish[1][i] + stepvol * random() * (fish[1][i] - Bari[i]))
        if isInBoundaries(x1):
            fish.append(x1)
            del fish[0]


def breedingOperator():
    if len(breedingCandidates) > 1:
        for i in breedingCandidates:
            ratio = 0
            index = 0
            for j in breedingCandidates:    #searches for best match to breed with each breeding candidate
                if i != j:
                    dis = distance(x[j][1], x[i][1])
                    if dis > 0:
                        temp = W[j][1] / dis        #calculate weight/distance ratio
                        if temp > ratio:
                            ratio = temp
                            index = j
            weight = (W[i][1] + W[index][1]) / 2    #weight of new fish
            position = [(x[i][1][0] + x[index][1][0]) / 2, (x[i][1][1] + x[index][1][1]) / 2]   #position of new fish
            W.append([weight, weight])
            x.append([position, position])
            minWeight = W[0][1]
            minIndex = 0
            for fish in range(len(W)):      #searches for lightest fish
                if W[fish][1] < minWeight:
                    minWeight = W[fish][1]
                    minIndex = fish
            del W[minIndex]     #deletes lightest fish
            del x[minIndex]
            del breedingCandidates[breedingCandidates.index(i)] #removes the fishes that just bred with wach ither from the breeding candidates list
            del breedingCandidates[breedingCandidates.index(index)]

#==============================================================#


def distance(x1, x2):   #calculates eucledian distance between two points
    return ((x1[0] - x2[0]) ** 2 + (x1[1] - x2[1]) ** 2) ** 0.5


def isInBoundaries(p):  #checks if a point is within the aquarium boundaries
    if aquarium_boundary >= p[0] >= -aquarium_boundary and aquarium_boundary >= p[1] >= -aquarium_boundary:
        return True
    return False


def firstMovement():    #makes the first individual movement of the fish school
    for fish in x:
        while (len(fish) < 2):
            x1 = []
            for i in range(2):
                x1.append(fish[0][i] + stepind * random() * (-1) ** randrange(2))
            if isInBoundaries(x1):
                fish.append(x1)



def parameterInitialization(fitnessFunc):
    global aquarium_boundary
    global stepvol
    global stepind
    global stepvolfin
    global stepindfin

    if fitnessFunc == "Sphere":
        aquarium_boundary = 200
        stepvol = aquarium_boundary * 0.02
        stepind = aquarium_boundary * 0.2
        stepvolfin = 0.002
        stepindfin = 0.002

    elif fitnessFunc == "Beale":
        aquarium_boundary = 4.5
        stepvol = aquarium_boundary * 0.2
        stepind = aquarium_boundary * 0.02
        stepvolfin = 0.00002
        stepindfin = 0.0001

    elif fitnessFunc == "Booth":
        aquarium_boundary = 10
        stepvol = aquarium_boundary * 0.02
        stepind = aquarium_boundary * 0.1
        stepvolfin = 0.00001
        stepindfin = 0.0001

    elif fitnessFunc == "Ackley":
        aquarium_boundary = 5
        stepvol = aquarium_boundary * 0.2
        stepind = aquarium_boundary * 0.2
        stepvolfin = 0.0002
        stepindfin = 0.002

    elif fitnessFunc == "Rosenbrock":
        aquarium_boundary = 10
        stepvol = aquarium_boundary * 0.2
        stepind = aquarium_boundary * 0.02
        stepvolfin = 0.0001
        stepindfin = 0.001



def nextOperator():     #executes the next operator
    global counter
    global cycles
    if counter > 4:
        print(cycles)
        cycles += 1
        if cycles == 5000:
            exit()
        counter = 0
        updateSteps()
    if counter == 0:
        individualMovement()
    elif counter == 1:
        feedingOperator()
    elif counter == 2:
        collectiveInstinctiveMovement()
    elif counter == 3:
        collectiveVolitiveMovement()
    elif counter == 4:
        breedingOperator()
    counter += 1


def iteration(frame):
    nextOperator()
    x1 = []
    x2 = []
    for fish in x:
        x1.append(fish[1][0])
        x2.append(fish[1][1])
    ln.set_data(x1, x2)
    return ln,


def updateSteps():

    global stepind
    global stepvol
    global stepindfin
    global stepvolfin
    if stepind > aquarium_boundary * stepindfin:
        stepind *= 0.98
    if stepvol > aquarium_boundary * stepvolfin:
        stepvol *= 0.9


fitnessFunction = "Booth"           #initialize which fitness function will be used
parameterInitialization(fitnessFunction)


noOfFishes = 30
W = []
x = []
breedingCandidates = []
breedingThreshold = 400
weight_scale = 500
counter = 0
cycles = 0

for i in range(noOfFishes):
    x.append([])
    x[i].append([uniform(-aquarium_boundary, aquarium_boundary), uniform(-aquarium_boundary, aquarium_boundary)])   #initialize random positions of fishes
    W.append([0, weight_scale / 2])
firstMovement()
feedingOperator()

fig, ax = plt.subplots()
xdata, ydata = [], []
ln, = plt.plot([], [], 'o')
ax.axhline(y=0, color='k')
ax.axvline(x=0, color='k')


def init():
   ax.set_xlim(-aquarium_boundary, aquarium_boundary)
   ax.set_ylim(-aquarium_boundary, aquarium_boundary)
   return ln,


ani = FuncAnimation(fig, iteration, interval = 1, init_func=init, blit=True)
plt.show()