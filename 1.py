import numpy as np
from scipy.optimize import root,fsolve
iters = 2000
angle = 40
y_axis = 0
z_axis = 10000
x_axis = 0

numGrp = 2*np.arange(2,24,1)

for num in numGrp:
    calu = np.zeros(3)
    for i in range(iters):
        myint = np.random.randint(0,10000)
        np.random.seed(myint)

        Paxis = np.zeros((num,3))
        Daxis = np.zeros((num,3))

        xyz0 = np.array((x_axis,y_axis,z_axis)).squeeze()

        dX = 0.03*100
        dY = 0.03*100
        dZ = 0.05*100
        dL = 0.02*100
        # print(x_axis)
        # print(y_axis)
        # print(z_axis)
        Daxis[:,0] = np.random.normal(0,dX,(num,))
        Daxis[:,1] = np.random.normal(0,dY,(num,))
        Daxis[:,2] = np.random.normal(0,dZ,(num,))
        DL = np.random.normal(0,dL,(num,))
        #print(DL)

        dx0 = z_axis*np.tan(angle*np.pi/180)
        dx1 = dx0 + 2500;
        #print(dx0, dx1)

        Paxis[:,0] = np.random.uniform(dx0, dx1,size =(num,))
        #print(Paxis)
        Paxis[:,1] = 0 #np.random.uniform(0, 0,size =(num,))
        Paxis[:,2] = 0 #np.random.uniform(0, 0,size =(num,))

        distL0 = np.sqrt( (x_axis-Paxis[:,0])**2 + (y_axis-Paxis[:,1])**2+ (z_axis-Paxis[:,2])**2 )
        #print(distL0)
        distL1 = distL0 + DL
        distL1 = distL1.reshape(num,1)

        Paxis += Daxis

        results = []


        def f5(x):
            m = np.array((x[0],x[1],x[2]),dtype = np.float64).reshape(1,3)
            minus = (m - Paxis)
        #    print('###################')
        #    print(Paxis)
        #    print('#####################')
            k = np.concatenate((minus**2, -distL1**2),axis = 1)
            l = np.sum(k,keepdims = True,axis =1).reshape(-1,1)
            return np.sum(l * minus,axis = 0)

        #Sxyz_tmp = root(f5,[x_axis,y_axis,z_axis])

        Sxyz = fsolve(f5,[x_axis,y_axis,z_axis])


        Dxyz = (Sxyz - xyz0)**2
        calu += Dxyz

    calu /= iters

    print(num, np.sqrt(calu), np.sqrt(np.sum(calu)),np.sqrt(np.sum(calu[0:2])))