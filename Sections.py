from shapely import Polygon, LineString, intersection
from shapely.ops import split
from shapely.affinity import rotate
import numpy as np

import numpy as np

def rotate_points(points, line_point, angle):
    """
    Rotate points around a given line by a specified angle.
    
    :param points: List of points [(x1, y1), (x2, y2), ...]
    :param line_point: A point on the line (x, y) that serves as the pivot for rotation
    :param angle: Angle of rotation in degrees
    :return: List of rotated points [(x1', y1'), (x2', y2'), ...]
    """
    # Convert angle to radians
    theta = np.radians(angle)
    
    # Define rotation matrix
    rotation_matrix = np.array([
        [np.cos(theta), -np.sin(theta)],
        [np.sin(theta), np.cos(theta)]
    ])
    
    # Translate points to align with the line
    translated_points = np.array(points) - np.array(line_point)
    
    # Rotate the points
    rotated_points = np.dot(translated_points, rotation_matrix.T)
    
    # Translate back to original position
    final_points = rotated_points + np.array(line_point)
    
    return final_points

class Sections():
    
    def __init__(self,Fc,Fy,E_s,A_steel,d_steel,Poly,Ref):
        self.fc = Fc                    ##f'c
        self.fy = Fy                    ##fy
        self.Es = E_s                   ##Modulo de elasticidad acero
        self.As = A_steel               #
        self.ds = d_steel
        self.n_s = len(self.As)
        self.poly = Poly
        self.ref = Ref
        self.c = None
        self.Cc = None
        self.d = None
        self.et = None
        self.es = None
        self.fs = None
        self.Fs = None
        self.b1 = None
        self.a = None
        self.dt = None
        self.Pn = None
        self.Mn = None
        self.fig1 = None
        self.fig2 = None
        self.Ac = None
        section = Polygon(self.poly)
        self.section = section
        self.Ac = self.section.area

    def Rotated_properties(self,alpha=0):
        
        
        self.section_rot = rotate(self.section,alpha,(0,0))
        self.ds_rot = rotate_points(self.ds,(0,0),alpha)
        
        self.cen_c = abs(self.section.centroid.y)

    
    def Cent_plas_rec(self,fc,fy,b,h,n_s,A_s,d_s):
        h_c = h/2
        
        F = 0.85*fc*(b)*(h)
        M = F*h_c

        for i in range(n_s):
            F_o =  A_s[i]*fy
            M_o =  A_s[i]*fy*d_s[i]
            F = F + F_o
            M = M + M_o

        return M/F

    def Cent_plas(self): #Ac y cen_c son Area de concreto y centroide de area de concreto
        F = 0.85*self.fc*self.Ac
        M = F*self.cen_c

        for i in range(self.n_s):
            F_o =  self.As[i]*self.fy
            M_o =  self.As[i]*self.fy*self.ds[i]
            F = F + F_o
            M = M + M_o

        self.CP = M/F
        return 0

    def Comp_defo(self,c,alpha=0):

        #theta = np.radians(alpha)

        ## Hallando el bloque de Whitney
        #-----------------------
        if(175<=self.fc<=280):
            β1 = 0.85
        elif(280<self.fc<550):
            β1 = 0.85-0.05*(self.fc-280)/70
        elif(self.fc>=550):
            β1 = 0.65

        self.b1 = β1
        a = β1*c
        self.a = a

        x_s = []
        y_s = []
        self.poly_rot = rotate_points(self.poly,(0,0),alpha)
        for coord in self.poly_rot:
            x_s.append(coord[0])
            y_s.append(coord[1])
        x_min = min(x_s)
        x_max = max(x_s)
        y_comp = max(y_s)
        self.ycomp = y_comp

        self.line = LineString([(x_min,y_comp-a),(x_max,y_comp-a)])

        y_witmax = 0

        for geom in split(self.section_rot,self.line).geoms:
            if max(geom.exterior.xy[1])>=y_witmax:
                self.witney_sec_rot = geom
                y_witmax = max(geom.exterior.xy[1])
        
        #self.witney_sec_rot = split(self.section_rot,self.line).geoms[0]

        self.line_c_rot = intersection(self.section_rot,self.line)
        #-----------------------

        #Rotate again to get the right coordinates (Xc,Yc)
        #-----------------------
        self.witney_sec = rotate(self.witney_sec_rot,-alpha,(0,0))

        Aw = self.witney_sec.area
        self.Aw = Aw
        self.Xc = self.witney_sec.centroid.x
        self.Yc = self.witney_sec.centroid.y

        self.line_c = rotate(self.line_c_rot,-alpha,(0,0))

        #-----------------------   

        #Finding distances of steels from extreme fiber in compression
        #-----------------------
        
        dis_st = []
        for i in range(self.n_s):
            dis_st.append(y_comp-self.ds_rot[i][1])
        self.dis_st = dis_st
        #-----------------------

        #Some params
        #-----------------------
        self.dt = max(self.dis_st)

        et = abs(0.003*(c-self.dt)/c)
        ey = self.fy/self.Es
        self.ey = ey
        #-----------------------

        #Calcs
        #-----------------------

        Cc = 0.85*self.fc*Aw/1000
        self.Cc = Cc

        Pn = Cc
        Mnx = Cc*(self.Yc)
        Mny = Cc*(self.Xc)

        es_l = []
        fs_l = []
        Fs_l = []
        for i in range(self.n_s):
            es = round(0.003*(c-self.dis_st[i])/c,4)
            es_l.append(es)
            if(abs(es)<ey):
                fs = self.Es*es
            elif(abs(es)>=ey):
                fs = np.sign(es)*self.fy
            fs_l.append(fs)
            Fs = fs*self.As[i]/1000
            Fs_l.append(Fs)

            Pn = Pn + Fs
            Mnx = Mnx + Fs*(self.ds[i][1])
            Mny = Mny + Fs*(self.ds[i][0])
            
        self.es = es_l
        self.fs = fs_l
        self.Fs = Fs_l
        self.c = c
        self.et = -et
        self.Pn = Pn
        self.Mnx = Mnx/100
        self.Mny = Mny/100

        #Reduction factor
        #-----------------------
        if(self.ref=="Espirales"):
            a = 0.75
            b = 0.15
        elif(self.ref=="Otro"):
            a = 0.65
            b = 0.25

        if(et<ey):
            fi = a
        elif(ey<et<0.005):
            fi =  a + b*(et-ey)/(0.005-ey)
        elif(et>=0.005):
            fi = 0.9
        
        self.fi = fi
        #-----------------------
        return 0
    
    def d_dis(self,c):
        d = 0
        at = 0
        for i in range(self.n_s):
            if(self.ds[i]>c):
                d = d + self.ds[i]*self.As[i]
                at = at + self.As[i]
        d = d/at
        self.d = d
        return 0
    
    def Plot_comp_defo(self):
        import matplotlib.pyplot as plt
        from shapely import Polygon, LineString, intersection
        from shapely.ops import split
        import numpy as np

        #self.section = rotate(self.section,-angle,(0,0))
        #self.witney_sec = rotate(self.witney_sec,-angle,(0,0))
        #self.line_c = rotate(self.line_c,-angle,(0,0))
        #self.ds = rotate_points(self.ds,(0,0),-angle)



        x, y = self.section.exterior.xy
        fig1,ax1 = plt.subplots()
        ax1.plot(x, y, color='blue')  # Plot the exterior of the polygon
        ax1.fill(x, y, color='lightblue', alpha=0.5)

        x,y = self.witney_sec.exterior.xy
        ax1.plot(x,y,color='red')
        ax1.fill(x, y, color='red', alpha=0.5)

        x,y = self.line_c.xy
        ax1.plot(x,y,color='red')

        x_steels = []
        y_steels = []
        for i in range(self.n_s):
            x_steels.append(self.ds[i][0])
            y_steels.append(self.ds[i][1])
        ax1.scatter(x_steels,y_steels)
        ax1.set_aspect('equal')

        #for i,a_s in enumerate(self.As):
        #    ax1.annotate(f"{str(round(a_s,2))}cm2", (zeros[i],ds_array[i]))

        plt.axis('off')

        #fig2,ax2 = plt.subplots()

        #y_s = []
        #for coord in self.poly:
        #    y_s.append(coord[1])
        #y_min = min(y_s)
        #y_max = max(y_s)
        #depth = y_max-y_min

        #depth = abs(self.ds[-1])

        #ax2.plot((0,0,0.003,self.et,0),(0,depth,depth,0,0))
        #ax2.annotate("ec=0.003",(0.0008,1*depth))
        #for i in range(self.n_s):
        #    ax2.plot((0,self.es[i]),(depth-self.ds[i],depth-self.ds[i]))
        #    ax2.annotate(f"{str(self.es[i])}",(self.es[i]*0.75,(depth-self.ds[i])))

        plt.axis('off')

        self.fig1 = fig1
        #self.fig2 = fig2
        plt.show()

    

    def Save_comp_defo(self,path):
        import matplotlib.pyplot as plt
        try:
            self.fig1.savefig(f"{path}/Beam sect")
            self.fig2.savefig(f"{path}/Beam deformation compatibility")
        except:
            print("No figure specified")
    
if __name__=="__main__":
    fc = 175 #kgf/cm2
    fy = 5000 #kgf/cm2
    Es = 2*10**6 #kgf/cm2
    As = [2,2,2,2] #cm2
    ds = [(0,24),(0,-24),(8,0),(-8,0)]
    b = 25
    h = 60
    poly = [(b/2,h/2),(b/2,-h/2),(-b/2,-h/2),(-b/2,h/2)]
    ref = "Otro"
    alpha = 15
    c = 9.5

    sec1 = Sections(fc,fy,Es,As,ds,poly,ref)
    sec1.Rotated_properties(alpha)
    #CP = sec1.Cent_plas()
    sec1.Comp_defo(c,alpha)
    print(sec1.ycomp)
    sec1.Plot_comp_defo()
    print(f"Pn:{sec1.Pn}")
    print(f"Mnx:{sec1.Mnx}")
    print(f"Mny:{sec1.Mny}")
    print(f"Cc:{sec1.Cc}")