from Sections import Sections
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as ptch
import pandas as pd
import plotly.graph_objects as go
from scipy.interpolate import griddata

path = "./Images"

As_steels = {"0":0,"3/8":0.71, "1/2":1.29, "5/8":2, "3/4":2.84, "7/8": 3.87, "1":5.1, "1 1/8":6.45}
diam_steels = {"0":0,"3/8":0.9525, "1/2":1.27, "5/8":1.5875, "3/4":1.905, "7/8": 2.2225, "1":2.54, "1 1/8":2.8575}

def calculate_offset_segment(p1, p2, width):
    delta = p2 - p1
    length = np.linalg.norm(delta)
    perp = np.array([-delta[1], delta[0]]) / length  # Perpendicular vector
    return p1 + (width / 2) * perp, p1 - (width / 2) * perp

def go_surface(x,y,z):
    X,Y = np.meshgrid(x,y)
    Z = griddata((x,y),z,(X,Y),method='linear')
    return go.Surface(x=X,y=Y,z=Z)

class Column():
    def __init__(self,L,fc,fy,Es):
        self.L = L
        self.fc = fc
        self.fy = fy
        self.Es = Es

    
    def Rec_col(self,b,h):
        poly = [(b/2,h/2),(b/2,-h/2),(-b/2,-h/2),(-b/2,h/2)]
        #poly_y = [(h/2,0),(h/2,-b),(-h/2,-b),(-h/2,0)]
        self.h = h
        self.b = b
        self.poly = poly
        #self.poly_y = poly_y
        self.name = "Columna Rectangular"
        pass

    def Circ_col(self,r):
        pass

    def Demand(self):
        pass
        #intention is to be able to find the flexure moment and axial force and the shear force

    class Col_section:
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

            diams_i = [d_long for i in range(n_x-2)]  ##diams i stands for the first row of diameters
            diams_i.append(d_corner)
            diams_i.insert(0,d_corner)
            d_i = self.rec + diam_steels[self.strp] + diam_steels[d_corner]/2  ##distancia de extremo superior a primer barra longitudinal
            diams_x = [diams_i]
            ds_x = [d_i] ##this is the distances from the compression fiber
            
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
       
        def Steel_coords_rec(self,n_x,n_y,d_corner,d_long):
            self.d_long = d_long
            self.d_corner = d_corner
            self.n_x = n_x
            self.n_y = n_y

            sx = (self.col.b-2*self.rec-2*diam_steels[self.strp]-diam_steels[d_corner])/(n_x-1)
            sy = (self.col.h-2*self.rec-2*diam_steels[self.strp]-diam_steels[d_long])/(n_y-1)
            self.sx = sx
            self.sy = sy

            ## Creating self.ds

            diams_i = [d_long for i in range(n_x-2)]  ##diams i stands for the first row of diameters
            diams_i.append(d_corner)
            diams_i.insert(0,d_corner)
            diams = diams_i

            d_ic = self.rec + diam_steels[self.strp] + diam_steels[d_corner]/2  ##distancia de extremo superior a primer barra (barra corner)
            d_il = self.rec + diam_steels[self.strp] + diam_steels[d_long]/2 ##distancia de extremo superior a primer barra (barra long)

            y_ic = self.col.h/2 - d_ic
            x_ic = -self.col.b/2 + d_ic

            y_il = self.col.h/2 - d_il
            x_il = -self.col.b/2 + d_il

            ds = [] ##this is the coords for all the rebars (starting with the first row)

            ds.append((x_ic,y_ic))
            for i in range(n_x-2): 
                ds.append((x_ic + (i+1)*sx, y_il))
            ds.append((-x_ic,y_ic))
            
            if n_y>2:
                y_temp = y_ic
                for i in range(n_y-2):
                    diams.append(d_long)
                    diams.append(d_long)
                    y_temp = y_temp - sy
                    ds.append((x_il,y_temp))
                    ds.append((-x_il,y_temp))
            
            diams.append(d_corner)
            for i in range(n_x-2):
                diams.append(d_long)
            diams.append(d_corner)

            #diams_j = [d_long for i in range(n_x-2)]
            
            
            #diams.append(diams_j)

            self.diams_rebar = diams

            ds.append((x_ic,-y_ic))
            for i in range(n_x-2): 
                ds.append((x_ic + (i+1)*sx, -y_il))
            ds.append((-x_ic,-y_ic))


            #d_j = self.col.h - self.rec - diam_steels[self.strp] - diam_steels[d_corner]/2
            #ds_x.append(d_j)
            self.ds_rebar = ds

       
        def Create_section_beta(self):
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
            self.section_x.Geometric_properties()
            self.section_x.Cent_plas()
            self.section_y = col_section_y
            self.section_y.Geometric_properties()
            self.section_y.Cent_plas()

        def Create_section(self,alpha=0):

            A_steels = []

            print(self.diams_rebar)
            print(self.ds_rebar)

            for diam in self.diams_rebar:
                a_s = As_steels[diam]
                A_steels.append(a_s)

            self.A_steels = A_steels

            col_section = Sections(self.col.fc,self.col.fy,self.col.Es,self.A_steels,self.ds_rebar,self.col.poly,"Otro")
            self.section = col_section
            

        def Diag_inter_beta(self,path,Pux,Puy,Mux,Muy):

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
                if 0.65*Pon <fiPn_x[i]:
                    x_list.append(fiMn_x[i])
            
            Lim = [0.65*Pon for i in x_list]
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
                if 0.65*Pon <fiPn_y[i]:
                    x_list.append(fiMn_y[i])
            
            Lim = [0.65*Pon for i in x_list]
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
     
        def Diag_inter(self,alphas):
            

            #Limites
            Ast = sum(self.section.As)
            Pon = 0.8*(0.85*self.col.fc*(self.section.Ac-Ast)+self.col.fy*Ast)/1000
            Tn = Ast*self.col.fy

            fiPn_list = []
            fiMn_x_list = []
            fiMn_y_list = []

            df_list = []
            dffi_list = []

            for alpha in alphas:
                self.section.Rotated_properties(alpha)

                Pn = []
                Mn_x = []
                Mn_y = []
                fiPn = []
                fiMn_x = []
                fiMn_y = []

                for c in np.arange(0,2*self.col.h,1):
                    #self.section.Geometric_properties(alpha)
                    
                    self.section.Comp_defo(c,alpha)
                    pn = self.section.Pn
                    mnx = self.section.Mnx
                    mny = self.section.Mny
                    fi = self.section.fi
                    Pn.append(pn)
                    Mn_x.append(mnx)
                    Mn_y.append(mny)
                    fiPn.append(fi*pn)
                    fiMn_x.append(fi*mnx)
                    fiMn_y.append(fi*mny)

                x_list = []

                for i in range(len(fiPn)):
                    if 0.65*Pon <fiPn[i]:
                        x_list.append(fiMn_x[i])
                        fiPn[i] = 0.65*Pon

                data = {'Pn':Pn,'Mn_x':Mn_x,'Mn_y':Mn_y}
                data_fi = {'fiPn':fiPn,'fiMn_x':fiMn_x,'fiMn_y':fiMn_y}

                df = pd.DataFrame(data).round(decimals=2)
                df_fi = pd.DataFrame(data_fi).round(decimals=2)
                df_list.append(df)
                dffi_list.append(df_fi)

                fiPn_list.append(fiPn)
                fiMn_x_list.append(fiMn_x)
                fiMn_y_list.append(fiMn_y)
            
            self.test_df = df_list
            self.test_dffi = dffi_list
            
            Lim = [0.65*Pon for i in x_list]
            fig,ax = plt.subplots()
            ax.plot(Mn_x,Pn)
            ax.plot(fiMn_x,fiPn)
            ax.plot(x_list,Lim)
            #ax.plot(Mux,Pux,'ro')
            ax.grid()
            ax.set_xlabel("Mn (tonf-m)")
            ax.set_ylabel("Pn (tonf)")
            #ax.plot()

            self.fig_diag = fig  
            try:
                plt.savefig(f"{path}/Interaction diagram X")
            except:
                print("No figure is defined")

            fiPn_plot = []
            fiMnx_plot = []
            fiMny_plot = []

            for i in range(len(fiPn_list)):
                for j in range(len(fiPn_list[0])):
                    fiPn_plot.append(fiPn_list[i][j])
                    fiMnx_plot.append(fiMn_x_list[i][j])
                    fiMny_plot.append(fiMn_y_list[i][j])

            self.fiPnplot = fiPn_plot
            self.fiMnx_plot = fiMnx_plot
            self.fiMny_plot = fiMny_plot
            surface_points = np.array(list(zip(fiMnx_plot, fiMny_plot, fiPn_plot)))
            self.surface_points = surface_points
            

            surface = go.Surface(x=fiMn_x_list,y=fiMn_y_list,z=fiPn_list,opacity=0.55,colorscale='blues',name='Interaction Surface')
            #surface = go_surface(fiMnx_plot,fiMny_plot,fiPn_plot)
            #print(fiPn_list)
            #print(fiMn_x_list)
            #print(fiMn_y_list)
            scatter = go.Scatter3d(x=np.array(fiMn_x_list).flatten(),y=np.array(fiMn_y_list).flatten(),z=np.array(fiPn_list).flatten(),mode='markers',marker=dict(size=3,color='black'))

            fig_go = go.Figure(data=[surface])

            fig_go.update_layout(title=dict(text='Interaction Surface'), autosize=True,
                  width=750, height=600,)
            
            #fig_go.show()
            self.go_fig = fig_go


        def Plot_rec_col(self,path):
            import matplotlib.pyplot as plt
            from shapely import Polygon, LineString, intersection
            from shapely.ops import split
            import numpy as np
            import matplotlib.patches as ptch
            import matplotlib

            

            #section = Polygon(self.section_x.poly)
            x, y = self.section.section.exterior.xy
            fig,ax = plt.subplots()
            ax.plot(x, y, color='#a7ada0')  # Plot the exterior of the polygon
            ax.fill(x, y, color='#8b8f86', alpha=0.5)

            d_strp = diam_steels[self.strp]
            
            b = self.col.b
            h = self.col.h
            rec = self.rec
            xb = b/2-rec-d_strp/2
            yb = h/2-rec-d_strp/2

            xb1 = xb-d_strp/2
            yb1 = yb-d_strp/2
            xb2 = xb+d_strp/2
            yb2 = yb+d_strp/2
            poly1 = np.array([[xb1,yb1],[xb1,-yb1],[-xb1,-yb1],[-xb1,yb1],[xb1,yb1]]) ##Stirrup
            poly2 = np.array([[xb2,yb2],[xb2,-yb2],[-xb2,-yb2],[-xb2,yb2],[xb2,yb2]])

            stp_patch2 = ptch.Polygon(poly2,color='#444740')

            stp_patch1 = ptch.Polygon(poly1,color='#C5C7C2')
            
            ax.add_patch(stp_patch2)
            ax.add_patch(stp_patch1)
            
            for i in range(len(self.diams_rebar)):
                bar = ptch.Circle((self.ds_rebar[i][0], self.ds_rebar[i][1]),radius=diam_steels[self.diams_rebar[i]]/2,color="#444740")
                ax.add_patch(bar)

            

            coords = [[xb, yb], [xb, -yb], [-xb, -yb], [-xb, yb], [xb, yb]]


            r_cor = diam_steels[self.d_corner]
            
           
            
            #ax.add_patch(box_patch)

            #stirrup = matplotlib.lines.Line2D([xb,xb,-xb,-xb,xb],[yb,-yb,-yb,yb,yb],linewidth = d_strp,color="#444740")

            #ax.add_line(stirrup)
            

            ax.set_aspect('equal', adjustable='box')
            plt.axis('off')

            #plt.show()

            self.fig_section = fig  
            try:
                plt.savefig(f"{path}/Column Section X")
            except:
                print("No figure is defined")

            

if __name__ == '__main__':
    L = 2.3 #m
    fc = 210
    fy = 4200
    Es = 2*10**6
    col1 = Column(L,fc,fy,Es)
    b = 30
    h = 60
    col1.Rec_col(b,h)
    col1_sec = col1.Col_section(col1,"5/8")
    col1_sec.Steel_coords_rec(4,3,"1","1")
    col1_sec.Create_section()
    alphas = []
    step = 15
    for alpha in np.arange(0,90+step,step):
        alphas.append(alpha)
    col1_sec.Diag_inter(alphas)
    #col1_sec.Plot_rec_col(path)
    #plt.show()