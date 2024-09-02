from Packages.Sections import Sections
import numpy as np
import copy
import matplotlib.pyplot as plt

As_steels = {"0":0,"3/8":0.71, "1/2":1.29, "5/8":2, "3/4":2.84, "7/8": 3.87, "1":5.1, "1 1/8":6.45}
diam_steels = {"0":0,"3/8":0.9525, "1/2":1.27, "5/8":1.5875, "3/4":1.905, "7/8": 2.2225, "1":2.54, "1 1/8":2.8575}

class Column():
    def __init__(self,L,fc,fy,Es):
        self.L = L
        self.fc = fc
        self.fy = fy
        self.Es = Es
        #self.b = None
        #self.h = None
        #self.R = None
        #self.poly = None
        #self.name = None
    
    def Rec_col(self,b,h):
        poly_x = [(b/2,0),(b/2,-h),(-b/2,-h),(-b/2,0)]
        poly_y = [(h/2,0),(h/2,-b),(-h/2,-b),(-h/2,0)]
        self.h = h
        self.b = b
        self.poly_x = poly_x
        self.poly_y = poly_y
        self.name = "Columna Rectangular"
        pass

    def Circ_col(self,r):
        pass

    def Demand(self):
        pass
        #intention is to be able to find the flexure moment and axial force and the shear force



class Column_section:
    def __init__(self,Column,d_stirrup,recover=4):
        self.col = Column
        self.strp = d_stirrup
        self.rec = recover

    def Steel_distribution_rec(self,n_x,n_y,d_corner,d_long):
        self.d_long = d_long
        self.d_corner = d_corner
        self.n_x = n_x
        self.n_y = n_y

        sx = (self.col.b-2*self.rec-2*diam_steels[self.strp]-diam_steels[d_corner])/(n_x-1)
        sy = (self.col.h-2*self.rec-2*diam_steels[self.strp]-diam_steels[d_long])/(n_y-1)
        self.sx = sx
        self.sy = sy

        ###------------------------------------------------------------------------------------------
        ###RESPECTO A X
        ###------------------------------------------------------------------------------------------

        diams_i = [d_long for i in range(n_x-2)]
        diams_i.append(d_corner)
        diams_i.insert(0,d_corner)
        d_i = self.rec + diam_steels[self.strp] + diam_steels[d_corner]/2
        diams_x = [diams_i]
        ds_x = [d_i]
        
        if n_y>2:
            d_temp = d_i
            for i in range(n_y-2):
                diams_x.append([d_long,d_long])
                d_temp = d_temp + sy
                ds_x.append(d_temp)
            
        
        diams_j = [d_long for i in range(n_x-2)]
        diams_j.append(d_corner)
        diams_j.insert(0,d_corner)
        diams_x.append(diams_j)

        self.diams_x = diams_x

        d_j = self.col.h - self.rec - diam_steels[self.strp] - diam_steels[d_corner]/2
        ds_x.append(d_j)
        self.ds_x = ds_x


        ###------------------------------------------------------------------------------------------
        ###RESPECTO A Y
        ###------------------------------------------------------------------------------------------

        diams_i = [d_long for i in range(n_y-2)]
        diams_i.append(d_corner)
        diams_i.insert(0,d_corner)
        d_i = self.rec + diam_steels[self.strp] + diam_steels[d_corner]/2
        diams_y = [diams_i]
        ds_y = [d_i]

        if n_x>2:
            d_temp = d_i
            for i in range(n_x-2):
                diams_y.append([d_long,d_long])
                d_temp = d_temp + sx
                ds_y.append(d_temp)

        

        diams_j = [d_long for i in range(n_y-2)]
        diams_j.append(d_corner)
        diams_j.insert(0,d_corner)
        diams_y.append(diams_j)

        self.diams_y = diams_y

        d_j = self.col.b - self.rec - diam_steels[self.strp] - diam_steels[d_corner]/2
        ds_y.append(d_j)
        self.ds_y = ds_y


    def Create_section(self):
        poly_x = self.col.poly_x
        poly_y = self.col.poly_y
        
        As_x = []
        As_y = []

        for diams in self.diams_x:
            a_s = 0
            for diam in diams:
                a_s = a_s + As_steels[diam]
            As_x.append(a_s)

        for diams in self.diams_y:
            a_s = 0
            for diam in diams:
                a_s = a_s + As_steels[diam]
            As_y.append(a_s)

        self.As_x = As_x
        self.As_y = As_y

        col_section_x = Sections(self.col.fc,self.col.fy,self.col.Es,self.As_x,self.ds_x,poly_x,"Otro")
        col_section_y = Sections(self.col.fc,self.col.fy,self.col.Es,self.As_y,self.ds_y,poly_y,"Otro")
        self.section_x = col_section_x
        self.section_y = col_section_y

    def Diag_inter(self,path,Pux,Puy,Mux,Muy):

        ###------------------------------------------------------------------------------------------
        ###RESPECTO A X
        ###------------------------------------------------------------------------------------------

        Pn_x = []
        Mn_x = []
        fiPn_x = []
        fiMn_x = []

        for c in np.arange(0.1,2*self.col.h,0.1):
            self.section_x.Comp_defo(c)
            pn = self.section_x.Pn
            mn = self.section_x.Mn
            fi = self.section_x.fi
            Pn_x.append(pn)
            Mn_x.append(mn)
            fiPn_x.append(fi*pn)
            fiMn_x.append(fi*mn)
        
        self.Pn_x = Pn_x
        self.Mn_x = Mn_x
        self.fiPn_x = fiPn_x
        self.fiMn_x = fiMn_x

        #Limites
        Ast = sum(self.section_x.As)
        Pon = 0.8*(0.85*self.col.fc*(self.section_x.Ac-Ast)+self.col.fy*Ast)/1000
        Tn = Ast*self.col.fy


        x_list = []

        for i in range(len(fiPn_x)):
            if 0.75*Pon <fiPn_x[i]:
                x_list.append(fiMn_x[i])
        
        Lim = [0.75*Pon for i in x_list]
        #print(Pon)
        fig,ax = plt.subplots()
        ax.plot(Mn_x,Pn_x)
        ax.plot(fiMn_x,fiPn_x)
        ax.plot(x_list,Lim)
        ax.plot(Mux,Pux,'ro')
        ax.grid()
        ax.set_xlabel("Mn (tonf-m)")
        ax.set_ylabel("Pn (tonf)")
        #ax.plot()

        self.fig_diag_x = fig  
        try:
            plt.savefig(f"{path}/Interaction diagram X")
        except:
            print("No figure is defined")


        ###------------------------------------------------------------------------------------------
        ###RESPECTO A Y
        ###------------------------------------------------------------------------------------------

        Pn_y = []
        Mn_y = []
        fiPn_y = []
        fiMn_y = []

        for c in np.arange(0.1,2*self.col.b,0.1):
            self.section_y.Comp_defo(c)
            pn = self.section_y.Pn
            mn = self.section_y.Mn
            fi = self.section_y.fi
            Pn_y.append(pn)
            Mn_y.append(mn)
            fiPn_y.append(fi*pn)
            fiMn_y.append(fi*mn)
        
        self.Pn_y = Pn_y
        self.Mn_y = Mn_y
        self.fiPn_y = fiPn_y
        self.fiMn_y = fiMn_y

        #Limites
        Ast = sum(self.section_y.As)
        Pon = 0.8*(0.85*self.col.fc*(self.section_y.Ac-Ast)+self.col.fy*Ast)/1000
        Tn = Ast*self.col.fy


        x_list = []

        for i in range(len(fiPn_y)):
            if 0.75*Pon <fiPn_y[i]:
                x_list.append(fiMn_y[i])
        
        Lim = [0.75*Pon for i in x_list]
        #print(Pon)
        fig,ax = plt.subplots()
        ax.plot(Mn_y,Pn_y)
        ax.plot(fiMn_y,fiPn_y)
        ax.plot(x_list,Lim)
        ax.plot(Muy,Puy,'ro')
        ax.grid()
        ax.set_xlabel("Mn (tonf-m)")
        ax.set_ylabel("Pn (tonf)")
        #ax.plot()

        self.fig_diag_y = fig  
        try:
            plt.savefig(f"{path}/Interaction diagram Y")
        except:
            print("No figure is defined")


    def Diag_inter_y(self,path):
        Pn = []
        Mn = []
        fiPn = []
        fiMn = []


    def Plot_rec_col(self,path):
        import matplotlib.pyplot as plt
        from shapely import Polygon, LineString, intersection
        from shapely.ops import split
        import numpy as np
        import matplotlib.patches as ptch
        import matplotlib

        section = Polygon(self.section_x.poly)
        x, y = section.exterior.xy
        fig,ax = plt.subplots()
        ax.plot(x, y, color='#a7ada0')  # Plot the exterior of the polygon
        ax.fill(x, y, color='#8b8f86', alpha=0.5)
        ds_array = -1*np.array(self.section_x.ds)
        


        for i in range(len(self.diams_x)):
            if i==0 or i==len(self.diams_x)-1:
                xo = -self.col.b/2 + self.rec + diam_steels[self.strp] + diam_steels[self.d_corner]/2
                sx = self.sx
            else:
                xo = -self.col.b/2 + self.rec + diam_steels[self.strp] + diam_steels[self.d_long]/2
                sx = (self.col.b-2*self.rec-2*diam_steels[self.strp]-diam_steels[self.d_long])
            for j in range(len(self.diams_x[i])):
                xpos = xo + j*sx
                bar = ptch.Circle((xpos, ds_array[i]),radius=diam_steels[self.diams_x[i][j]]/2,color="#444740")
                ax.add_patch(bar)


        d_strp = diam_steels[self.strp]
        stirrup = matplotlib.lines.Line2D([self.col.b/2-self.rec-d_strp,self.col.b/2-self.rec-d_strp,-self.col.b/2+self.rec+d_strp,-self.col.b/2+self.rec+d_strp,self.col.b/2-self.rec-d_strp],[-self.rec-d_strp,-self.col.h+self.rec+d_strp,-self.col.h+self.rec+d_strp,-self.rec-d_strp,-self.rec-d_strp],linewidth = d_strp,color="#444740")

        ax.add_line(stirrup)

        ax.set_aspect('equal')
        plt.axis('off')

        #plt.show()

        self.fig_section = fig  
        try:
            plt.savefig(f"{path}/Column Section")
        except:
            print("No figure is defined")

    def As_lim(self):
        Asmin = 0.01*self.section_x.Ac
        Asmax = 0.06*self.section_x.Ac
        As_tot = sum(self.As_x)

        self.Asmin = Asmin
        self.Asmax = Asmax
        self.As_tot = As_tot


if __name__ == '__main__':
    col = Column(2.3,210,4200,2*10**6)
    col.Rec_col(30,60)
    col_sec = Column_section(col,"3/8")
    col_sec.Steel_distribution_rec(4,3,"1","1")
    print(col_sec.diams)
    col_sec.Create_section()
    #col_sec.section.Comp_defo(40)
    #print(col_sec.section.poly)
    #print(col_sec.section.Mn)
    #print(col_sec.section.Pn)
    col_sec.Diag_inter()
    col_sec.Plot_rec_col()
    print(col_sec.diams)
    #col_sec.section.Comp_defo(10)
    #col_sec.section.Plot_comp_defo()
    #print(col_sec.As)
    #print(col_sec.ds)
    plt.show()
    #section_1 = Sections(210,4200,2*10**6,[5.1*4,5.1*2,5.1*4],[],Poly,Ref)

