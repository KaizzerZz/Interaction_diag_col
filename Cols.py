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
