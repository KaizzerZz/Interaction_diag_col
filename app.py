#from Packages.Sections import Sections
from Packages.Cols import Column,Column_section
import streamlit as st
import numpy as np
#import subprocess
#from jinja2 import Environment, FileSystemLoader
#import os

path = "./Images"

def Define_column(fc,fy,Es,b,h,strr,n_x,n_y,d_corner,d_long,path,Pux,Puy,Mux,Muy):
    L = 1#m
    
    col = Column(L,fc,fy,Es)
    col.Rec_col(b,h)
    col_sec = Column_section(col,strr)
    col_sec.Steel_distribution_rec(n_x,n_y,d_corner,d_long)
    col_sec.Create_section()
    col_sec.Diag_inter(path,Pux,Puy,Mux,Muy)
    col_sec.Plot_rec_col(path)
    col_sec.As_lim()

    return col_sec

def Generate_report_web(column):

    st.title("Diagrama de interacción en dirección X")
    st.image('./Images/Interaction diagram X.png',use_column_width="auto",caption="Diagrama de interacción X")

    st.title("Diagrama de interacción en dirección Y")
    st.image('./Images/Interaction diagram Y.png',use_column_width="auto",caption="Diagrama de interacción Y")


with st.sidebar:
    st.title("1. Propiedades geométricas")

    b = st.number_input("Base (cm)",value=30)

    h = st.number_input("Peralte (cm)",value=50)

    st.title("2. Propiedades de los materiales")

    fc = st.number_input("Resistencia a la compresión del concreto (kg/cm2)",value=210)

    fy = st.number_input("Resistencia a la fluencia del acero (kg/cm2)",value=4200)

    Es = st.number_input("Modulo de elasticidad del acero (kg/cm2)",value=2*10**6)

    st.title("3. Acero transversal")

    strr = st.selectbox("Diámetro de estribo (pulg)",("0","3/8","1/2","5/8","3/4","7/8","1","1 1/8"),index=1)

    st.title("4. Acero longitudinal")

    n_x = st.number_input("Número de barras en X",value=3)

    n_y = st.number_input("Número de barras en Y",value=4)

    d_corner = st.selectbox("Diámetro de barra en esquinas (pulg)",("0","3/8","1/2","5/8","3/4","7/8","1","1 1/8"),index=3)

    d_long = st.selectbox("Diámetro de barra en en filas (pulg)",("0","3/8","1/2","5/8","3/4","7/8","1","1 1/8"),index=3)

    st.title("5. Momento y carga axial ultima")

    st.subheader("Respecto a eje x")

    Pux = st.number_input("Pu (tonf)",value=90)
    Mux = st.number_input("Mu (tonf-m)",value=9)

    st.subheader("Respecto a eje y")

    Puy = st.number_input("Pu (tonf)",value=100)
    Muy = st.number_input("Mu (tonf-m)",value=10)

st.title("Diseño de columna a flexocompresión")

st.write("Para definir la sección de la columna, se deberá indicar sus propiedades geométricas, del material y de la distribución de aceros.")

col = Define_column(fc,fy,Es,b,h,strr,n_x,n_y,d_corner,d_long,path,Pux,Puy,Mux,Muy)

st.image('./Images/Column Section.png',use_column_width="auto",caption="Sección de la columna")

#st.title("Consideraciones de diseño")

st.header("Consideraciones de diseño",divider=True)

st.subheader("Acero mínimo y acero máximo")

st.latex(r"As_{total}" + f"= {str(round(col.As_tot,1))} cm^2")

st.latex(r"As_{min}" + f"= {str(round(col.Asmin,1))} cm^2" )

if col.Asmin<=col.As_tot:
    st.success("Cumple acero mínimo")
else:
    st.error("No cumple acero mínimo")

st.latex(r"As_{max}" + f"= {str(round(col.Asmax,1))} cm^2")

if col.Asmax>=col.As_tot:
    st.success("Cumple acero máximo")
else:
    st.error("No cumple acero máximo")

#st.markdown(r"$$As_{max}$$" + f"= {str(round(col.Asmax,1))}")


button = st.button("Generate Interaction diagrams")

if button:
    try:
        Generate_report_web(col)
        st.success("Report generated!")
    except:
        st.error("Error")
    button=False

